module SecondQuantizedAlgebraQuantumToolboxExt

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumToolboxBackend, Op
import QuantumToolbox as QTB
import SciMLOperators as SO
import SymbolicUtils: BasicSymbolic
import SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals, sparse, spzeros
import LinearAlgebra: dot, mul!

const QTBDims = Union{Integer, AbstractVector{<:Integer}, Tuple{Vararg{T}} where {T <: Integer}}

struct VecSum{T, F, O <: SO.AbstractSciMLOperator} <: SO.AbstractSciMLOperator{T}
    coeffs::F
    vals::Vector{T}
    ops::Vector{O}
end

struct ConstCoeffs{T}
    c::Vector{T}
end
(f::ConstCoeffs)(u, p, t) = f.c

struct ConjCoeffs{F}
    f::F
end
(f::ConjCoeffs)(u, p, t) = conj.(f.f(u, p, t))

struct ProductCoeffs{FA, FB}
    a::FA
    b::FB
end
(f::ProductCoeffs)(u, p, t) = outer_values(f.a(u, p, t), f.b(u, p, t))

outer_values(a::Vector, b::Vector) = [x * y for x in a for y in b]

conj_coeffs(f::ConstCoeffs) = ConstCoeffs(conj.(f.c))
conj_coeffs(f) = ConjCoeffs(f)
product_coeffs(a::ConstCoeffs, b::ConstCoeffs) = ConstCoeffs(outer_values(a.c, b.c))
product_coeffs(a, b) = ProductCoeffs(a, b)

static_vecsum(c::Vector{T}, ops::Vector) where {T} = VecSum(ConstCoeffs(c), copy(c), ops)

Base.size(L::VecSum) = size(L.ops[1])
Base.size(L::VecSum, i::Int) = size(L.ops[1], i)
Base.eltype(::VecSum{T}) where {T} = T
SO.issquare(::VecSum) = true
SO.isconstant(L::VecSum) = L.coeffs isa ConstCoeffs
function SO.update_coefficients!(L::VecSum, u, p, t)
    SO.isconstant(L) || copyto!(L.vals, L.coeffs(u, p, t))
    return L
end
SO.update_coefficients(L::VecSum{T}, u, p, t) where {T} =
    SO.isconstant(L) ? L : VecSum(L.coeffs, convert(Vector{T}, L.coeffs(u, p, t)), L.ops)

function mul!(w::AbstractVecOrMat, L::VecSum, v::AbstractVecOrMat)
    mul!(w, L.ops[1], v, L.vals[1], false)
    for i in 2:length(L.ops)
        mul!(w, L.ops[i], v, L.vals[i], true)
    end
    return w
end
function mul!(w::AbstractVecOrMat, L::VecSum, v::AbstractVecOrMat, α, β)
    mul!(w, L.ops[1], v, L.vals[1] * α, β)
    for i in 2:length(L.ops)
        mul!(w, L.ops[i], v, L.vals[i] * α, true)
    end
    return w
end

function Base.:*(L::VecSum, v::AbstractVector)
    w = similar(v, promote_type(eltype(L), eltype(v)), size(L, 1))
    return mul!(w, L, v)
end
function Base.:*(L::VecSum, v::AbstractMatrix)
    w = similar(v, promote_type(eltype(L), eltype(v)), size(L, 1), size(v, 2))
    return mul!(w, L, v)
end

(L::VecSum)(w, v, u, p, t) = (SO.update_coefficients!(L, u, p, t); mul!(w, L, v))
(L::VecSum)(w, v, u, p, t, α, β) = (SO.update_coefficients!(L, u, p, t); mul!(w, L, v, α, β))

function SO.concretize(L::VecSum{T}) where {T}
    rows, cols, vals = Int[], Int[], T[]
    for i in eachindex(L.ops)

        append_scaled!(rows, cols, vals, L.vals[i], sparse(SO.concretize(L.ops[i])))
    end
    return sparse(rows, cols, vals, size(L)...)
end

function append_scaled!(rows::Vector{Int}, cols::Vector{Int}, vals::Vector{T}, c::T, m::SparseMatrixCSC) where {T}
    rv, nz = rowvals(m), nonzeros(m)
    for j in axes(m, 2), k in nzrange(m, j)
        push!(rows, rv[k]); push!(cols, j); push!(vals, c * nz[k])
    end
    return rows
end
Base.convert(::Type{AbstractMatrix}, L::VecSum) = SO.concretize(L)

SO.cache_operator(L::VecSum, ::AbstractVecOrMat) = L

Base.transpose(L::VecSum) = VecSum(L.coeffs, copy(L.vals), [transpose(o) for o in L.ops])
Base.adjoint(L::VecSum) = VecSum(conj_coeffs(L.coeffs), conj.(L.vals), [adjoint(o) for o in L.ops])

QTB._spre(A::VecSum) = VecSum(A.coeffs, copy(A.vals), [QTB._spre(o) for o in A.ops])
QTB._spost(B::VecSum) = VecSum(B.coeffs, copy(B.vals), [QTB._spost(o) for o in B.ops])

opprod(x::SO.MatrixOperator, y::SO.MatrixOperator) = SO.MatrixOperator(x.A * y.A)
opprod(x, y) = x * y

Base.:*(A::VecSum{T}, B::VecSum{T}) where {T} = VecSum(
    product_coeffs(A.coeffs, B.coeffs),
    outer_values(A.vals, B.vals),
    [opprod(x, y) for x in A.ops for y in B.ops],
)

@inline function SQA.numeric_operator(be::QuantumToolboxBackend, op::Op, N::Integer)
    k = op.kind
    k === SQA.OP_DESTROY    && return QTB.destroy(N)
    k === SQA.OP_CREATE     && return QTB.create(N)
    k === SQA.OP_POSITION   && return QTB.position(N)
    k === SQA.OP_MOMENTUM   && return QTB.momentum(N)

    k === SQA.OP_TRANSITION && return QTB.projection(N, Int(op.l1) - 1, Int(op.l2) - 1)
    k === SQA.OP_PAULI      && return qtb_pauli(op, N)
    k === SQA.OP_SPIN       && return qtb_spin(op, N)
    throw(ArgumentError("no QuantumToolbox numeric_operator for role $k on dim $N"))
end

function qtb_pauli(op::Op, N::Integer)
    N == 2 || throw(ArgumentError("Pauli operators require dim 2, got $N"))
    op.l1 == 1 && return QTB.sigmax()
    op.l1 == 2 && return QTB.sigmay()
    return QTB.sigmaz()
end

function qtb_spin(op::Op, N::Integer)
    j = (N - 1) // 2
    op.l1 == 1 && return QTB.jmat(j, Val(:x))
    op.l1 == 2 && return QTB.jmat(j, Val(:y))
    return QTB.jmat(j, Val(:z))
end

SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.FockSpace, N) = Int(N) + 1
SQA.numeric_basis(::QuantumToolboxBackend, h::SQA.NLevelSpace, _) = Int(h.n)
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.PauliSpace, _) = 2
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.SpinSpace, s) = Int(2s + 1)
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.PhaseSpace, N) = Int(N) + 1
function SQA.numeric_basis(be::QuantumToolboxBackend, h::SQA.ProductSpace, dims)
    SQA.check_product_dims(h, dims)
    return ntuple(i -> Int(SQA.numeric_basis(be, h.spaces[i], dims[i])), length(h.spaces))
end

function SQA.numeric_subbasis(::QuantumToolboxBackend, N::Integer, slot::Int)
    slot == 1 || throw(ArgumentError("simple QuantumToolbox dimension has no subsystem slot $slot"))
    return N
end
function SQA.numeric_subbasis(::QuantumToolboxBackend, ds::AbstractVector, slot::Int)
    1 <= slot <= length(ds) || throw(
        ArgumentError(
            "QuantumToolbox dimensions have $(length(ds)) subsystems, not slot $slot",
        ),
    )
    return ds[slot]
end
function SQA.numeric_subbasis(::QuantumToolboxBackend, ds::Tuple{Vararg{T}}, slot::Int) where {T <: Integer}
    1 <= slot <= length(ds) || throw(
        ArgumentError(
            "QuantumToolbox dimensions have $(length(ds)) subsystems, not slot $slot",
        ),
    )
    return ds[slot]
end

SQA.numeric_num_subsystems(::QuantumToolboxBackend, ::Integer) = 1
SQA.numeric_num_subsystems(::QuantumToolboxBackend, ds::AbstractVector) = length(ds)
SQA.numeric_num_subsystems(::QuantumToolboxBackend, ds::Tuple{Vararg{T}}) where {T <: Integer} = length(ds)

SQA.numeric_embed(::QuantumToolboxBackend, N::Integer, slot::Int, m) = m
function SQA.numeric_embed(::QuantumToolboxBackend, ds::AbstractVector, slot::Int, m)
    length(ds) == 1 && return m
    return QTB.tensor((i == slot ? m : QTB.qeye(ds[i]) for i in eachindex(ds))...)
end
function SQA.numeric_embed(::QuantumToolboxBackend, ds::Tuple{Vararg{T}}, slot::Int, m) where {T <: Integer}
    return QTB.tensor((i == slot ? m : QTB.to_sparse(QTB.qeye(ds[i])) for i in eachindex(ds))...)
end

SQA.numeric_identity(::QuantumToolboxBackend, N::Integer) = QTB.qeye(N)
function SQA.numeric_identity(::QuantumToolboxBackend, ds::AbstractVector)
    length(ds) == 1 && return QTB.qeye(ds[1])
    return QTB.tensor((QTB.qeye(d) for d in ds)...)
end
SQA.numeric_identity(::QuantumToolboxBackend, ds::Tuple{Vararg{T}}) where {T <: Integer} =
    QTB.tensor((QTB.qeye(d) for d in ds)...)

function SQA.numeric_assemble(be::QuantumToolboxBackend, basis, terms)
    refid = SQA.numeric_identity(be, basis)
    coeffs = ComplexF64[c for (c, _) in terms]
    ops = [SO.MatrixOperator(fuse(be, basis, factors).data) for (_, factors) in terms]
    return QTB.QobjEvo(static_vecsum(coeffs, ops), QTB.Operator(), refid.dimensions)
end

function SQA.to_numeric_td(terms::Vector{SQA.CoeffTerm}, ctx::SQA.NumericContext{QuantumToolboxBackend}, time_parameter)
    leaves = SQA.leaf_cache(ctx)
    factor_lists = [SQA.term_factors(ops, ctx, leaves) for (_, ops) in terms]
    isempty(factor_lists) && push!(factor_lists, SQA.term_factors(Op[], ctx, leaves))
    be, basis = ctx.backend, ctx.basis
    ops = [SO.MatrixOperator(fuse(be, basis, factors).data) for factors in factor_lists]
    coeffs = SQA.compile_td_coeffs((c for (c, _) in terms), time_parameter)
    refid = SQA.numeric_identity(be, basis)
    return QTB.QobjEvo(VecSum(coeffs, zeros(ComplexF64, length(ops)), ops), QTB.Operator(), refid.dimensions)
end

function fuse(be::QuantumToolboxBackend, basis, factors)
    isempty(factors) && return QTB.to_sparse(SQA.numeric_identity(be, basis))
    length(factors) == 1 && return factors[1]
    return foldl(*, factors)
end

SQA.numeric_materialize(::QuantumToolboxBackend, op, ::SQA.EagerOperator) = QTB.to_sparse(qtb_eager(op))
function SQA.numeric_materialize(::QuantumToolboxBackend, op, op_type)
    op_type === identity && return op
    eager = qtb_eager(op)
    op_type === SO.concretize && return eager.data
    return op_type(eager)
end

qtb_eager(op::QTB.QuantumObjectEvolution) =
    QTB.Qobj(SO.concretize(op.data), QTB.Operator(), op.dimensions)
qtb_eager(op::QTB.QuantumObject) = op

SQA.to_numeric_lazy(op::Op, ctx::SQA.NumericContext{QuantumToolboxBackend}) =
    SQA.to_numeric_lazy(SQA.single_qadd(SQA.CNUM_ONE, Op[op]), ctx)

function expect_ket(L::VecSum, data::AbstractVector)
    work = similar(data, ComplexF64)
    mul!(work, L, data)
    return ComplexF64(dot(data, work))
end

function trace_product(A::SparseMatrixCSC, ρ::AbstractMatrix)
    value = 0.0im
    rows = rowvals(A)
    vals = nonzeros(A)
    for j in axes(A, 2), p in nzrange(A, j)
        value += vals[p] * ρ[j, rows[p]]
    end
    return value
end
function trace_product(A::AbstractMatrix, ρ::AbstractMatrix)
    value = 0.0im
    for j in axes(A, 2), i in axes(A, 1)
        value += A[i, j] * ρ[j, i]
    end
    return value
end

function expect_density(L::VecSum, data::AbstractMatrix)
    acc = Ref(0.0im)
    for i in eachindex(L.ops)
        add_trace!(acc, L.vals[i], L.ops[i], data)
    end
    return acc[]
end
function add_trace!(acc::Base.RefValue{ComplexF64}, c::ComplexF64, op, data::AbstractMatrix)
    acc[] += c * trace_product(SO.concretize(op), data)
    return acc
end

QTB.expect(
    numop::QTB.QuantumObjectEvolution{QTB.Operator, D, <:VecSum{ComplexF64}},
    state::QTB.QuantumObject{QTB.Ket},
) where {D} = expect_ket(numop.data, state.data)
QTB.expect(
    numop::QTB.QuantumObjectEvolution{QTB.Operator, D, <:VecSum{ComplexF64}},
    state::QTB.QuantumObject{QTB.Bra},
) where {D} = QTB.expect(numop, adjoint(state))
QTB.expect(
    numop::QTB.QuantumObjectEvolution{QTB.Operator, D, <:VecSum{ComplexF64}},
    state::QTB.QuantumObject{QTB.Operator},
) where {D} = expect_density(numop.data, state.data)

SQA.numeric_expect(
    ::QuantumToolboxBackend,
    numop::QTB.QuantumObjectEvolution{QTB.Operator, D, <:VecSum{ComplexF64}},
    state::QTB.QuantumObject,
) where {D} = QTB.expect(numop, state)

SQA.numeric_expect(
    ::QuantumToolboxBackend,
    numop::QTB.QuantumObject{QTB.Operator},
    state::QTB.QuantumObject{QTB.Ket},
) = ComplexF64(dot(state.data, numop.data, state.data))
SQA.numeric_expect(
    be::QuantumToolboxBackend,
    numop::QTB.QuantumObject{QTB.Operator},
    state::QTB.QuantumObject{QTB.Bra},
) = SQA.numeric_expect(be, numop, adjoint(state))
SQA.numeric_expect(
    ::QuantumToolboxBackend,
    numop::QTB.QuantumObject{QTB.Operator},
    state::QTB.QuantumObject{QTB.Operator},
) = ComplexF64(trace_product(numop.data, state.data))

SQA.numeric_expect(::QuantumToolboxBackend, numop, state::QTB.QuantumObject) =
    ComplexF64(QTB.expect(numop, state))
function SQA.numeric_expect(
        be::QuantumToolboxBackend,
        numop,
        states::AbstractVector{<:QTB.QuantumObject},
    )
    return ComplexF64[SQA.numeric_expect(be, numop, state) for state in states]
end
SQA.numeric_expect(::QuantumToolboxBackend, numop, states::AbstractVector) =
    ComplexF64[ComplexF64(QTB.expect(numop, s)) for s in states]

SQA.numeric_backend(::QTB.AbstractQuantumObject) = QuantumToolboxBackend()
qtb_numeric_basis(ds::Tuple{Int}) = ds[1]
qtb_numeric_basis(ds::Tuple{Vararg{Int}}) = ds
qtb_numeric_basis(ds) = qtb_numeric_basis(Tuple(ds))
SQA.numeric_basis(o::QTB.AbstractQuantumObject{QTB.Bra}) = qtb_numeric_basis(o.dims[2])
SQA.numeric_basis(o::QTB.AbstractQuantumObject) = qtb_numeric_basis(o.dims[1])

SQA.to_numeric(op::Op, dims::QTBDims, d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.static_operator(
    SQA.single_qadd(SQA.CNUM_ONE, Op[op]),
    SQA.NumericContext(QuantumToolboxBackend(), dims, d),
    SQA.EagerOperator(),
)
SQA.to_numeric(q::SQA.QAdd, dims::QTBDims, d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.static_operator(q, SQA.NumericContext(QuantumToolboxBackend(), dims, d), SQA.EagerOperator())
SQA.to_numeric(x::Number, dims::QTBDims, ::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.numeric_materialize(
    QuantumToolboxBackend(),
    SQA.to_complex(x) * SQA.numeric_identity(QuantumToolboxBackend(), dims),
    SQA.EagerOperator(),
)

SQA.to_numeric(
    q::SQA.QAdd, dims::QTBDims,
    sites::AbstractDict{Int, Vector{Int}},
    d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS,
    scalar_subs::AbstractDict = SQA.NO_SCALAR_SUBS,
) = SQA.to_numeric_indexed(QuantumToolboxBackend(), dims, q, sites, d, scalar_subs)

function SQA.to_numeric(
        op::SQA.QField, dims::QTBDims;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
    )
    return SQA.to_numeric_kw(QuantumToolboxBackend(), op, dims; parameter, time_parameter, operators, adjoint_ops, op_type)
end
function SQA.to_numeric(
        x::Union{Number, BasicSymbolic}, dims::QTBDims;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
    )
    return SQA.to_numeric_kw(QuantumToolboxBackend(), x, dims; parameter, time_parameter, operators, adjoint_ops, op_type)
end

SQA.to_numeric(ops::AbstractVector{T}, dims::D; kwargs...) where {T, D <: QTBDims} =
    [SQA.to_numeric(op, dims; kwargs...) for op in ops]

SQA.to_numeric(op, state::QTB.QuantumObject; kwargs...) = SQA.to_numeric(op, SQA.numeric_basis(state); kwargs...)
SQA.to_numeric(op::SQA.QField, state::QTB.QuantumObject, d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.to_numeric(op, SQA.numeric_basis(state), d)
SQA.to_numeric(x::Number, state::QTB.QuantumObject) = SQA.to_numeric(x, SQA.numeric_basis(state))

SQA.numeric_average(op, state::QTB.QuantumObject, d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.numeric_average_impl(op, state, d)

sandwich_matrix(left::Vector{Op}, right::Vector{Op}, ctx::SQA.NumericContext, cache::Dict{SQA.SideKey}) =
    kron(transpose(SQA.side_matrix(right, ctx, cache)), SQA.side_matrix(left, ctx, cache))

super_dimensions(be::QuantumToolboxBackend, basis) =
    QTB.spre(SQA.numeric_identity(be, basis)).dimensions

function SQA.numeric_superoperator(be::QuantumToolboxBackend, basis, terms::Vector{<:Tuple{ComplexF64, Any, Any}})
    n = size(SQA.numeric_identity(be, basis), 1)
    rows, cols, vals = Int[], Int[], ComplexF64[]
    for (c, left, right) in terms
        append_scaled!(rows, cols, vals, c, kron(transpose(right), left))
    end
    data = sparse(rows, cols, vals, n^2, n^2)
    return QTB.Qobj(data, QTB.SuperOperator(), super_dimensions(be, basis))
end

function SQA.to_numeric_td(S::SQA.SuperOperator, ctx::SQA.NumericContext{QuantumToolboxBackend}, time_parameter)
    cache = SQA.side_cache(ctx)
    ops = [SO.MatrixOperator(sandwich_matrix(t.left, t.right, ctx, cache)) for (t, _) in S.arguments]
    if isempty(ops)
        n = size(SQA.numeric_identity(ctx.backend, ctx.basis), 1)
        push!(ops, SO.MatrixOperator(spzeros(ComplexF64, n^2, n^2)))
    end
    coeffs = SQA.compile_td_coeffs(values(S.arguments), time_parameter)
    vecsum = VecSum(coeffs, zeros(ComplexF64, length(ops)), ops)
    return QTB.QobjEvo(vecsum, QTB.SuperOperator(), super_dimensions(ctx.backend, ctx.basis))
end

SQA.numeric_matrix(::QuantumToolboxBackend, op::QTB.QuantumObject) = op.data

SQA.numeric_scaled_td(::QuantumToolboxBackend, op::QTB.QuantumObject, f) = QTB.QobjEvo(op, SQA.lift_pt(f))

const QTBMap = Union{SQA.SuperOperator, SQA.SuperExpr, SQA.Lindbladian}

function SQA.to_numeric(
        S::QTBMap, dims::QTBDims;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
    )
    return SQA.to_numeric_kw(QuantumToolboxBackend(), S, dims; parameter, time_parameter, operators, adjoint_ops, op_type)
end

SQA.to_numeric(S::SQA.NumericMap, dims::QTBDims, d::AbstractDict{<:SQA.QSym}) =
    SQA.to_numeric_translated(
    S, SQA.NumericContext(QuantumToolboxBackend(), dims, d), SQA.SubstitutionRules(),
    SQA.NO_TIME_PARAMETER, SQA.EagerOperator(),
)
SQA.to_numeric(::SQA.FormalMap, ::QTBDims, ::AbstractDict{<:SQA.QSym}) = SQA.formal_map_numeric_error()
SQA.to_numeric(S::QTBMap, state::QTB.QuantumObject, d::AbstractDict{<:SQA.QSym}) =
    SQA.to_numeric(S, SQA.numeric_basis(state), d)

SQA.to_numeric(
    S::SQA.NumericMap, dims::QTBDims,
    sites::AbstractDict{Int, Vector{Int}},
    d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS,
    scalar_subs::AbstractDict = SQA.NO_SCALAR_SUBS,
) = SQA.to_numeric_indexed(QuantumToolboxBackend(), dims, S, sites, d, scalar_subs)

# QuantumToolbox exports these names; forwarding keeps `spre(a)` symbolic when both are loaded.
QTB.spre(A::SQA.QField) = SQA.spre(A)
QTB.spost(B::SQA.QField) = SQA.spost(B)
QTB.sprepost(A::SQA.QField, B::SQA.QField) = SQA.sprepost(A, B)
QTB.sprepost(A::SQA.QField, B::Number) = SQA.sprepost(A, B)
QTB.sprepost(A::Number, B::SQA.QField) = SQA.sprepost(A, B)
QTB.liouvillian(H::SQA.QField, c_ops::AbstractVector{<:SQA.QField}) = SQA.liouvillian(H, c_ops)
QTB.liouvillian(H::SQA.QField, c_ops::Tuple{Vararg{SQA.QField}}) = SQA.liouvillian(H, collect(c_ops))

end
