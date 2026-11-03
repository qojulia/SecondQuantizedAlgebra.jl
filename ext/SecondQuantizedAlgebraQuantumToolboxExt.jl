module SecondQuantizedAlgebraQuantumToolboxExt

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumToolboxBackend, Op
import QuantumToolbox as QTB
import SciMLOperators as SO
import SymbolicUtils: BasicSymbolic
using Symbolics: Num
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
(f::ConstCoeffs)(out::Vector, u, p, t) = copyto!(out, f.c)

struct ConjCoeffs{F}
    f::F
end
(f::ConjCoeffs)(u, p, t) = conj.(f.f(u, p, t))
function (f::ConjCoeffs)(out::Vector, u, p, t)
    f.f(out, u, p, t)
    out .= conj.(out)
    return out
end

struct ProductCoeffs{FA, FB}
    a::FA
    b::FB
end
(f::ProductCoeffs)(u, p, t) = outer_values(f.a(u, p, t), f.b(u, p, t))
function (f::ProductCoeffs)(out::Vector, u, p, t)
    left, right = f.a(u, p, t), f.b(u, p, t)
    k = 0
    for a in left, b in right
        k += 1
        out[k] = a * b
    end
    return out
end

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
    SO.isconstant(L) || L.coeffs(L.vals, u, p, t)
    return L
end
function SO.update_coefficients(L::VecSum, u, p, t)
    SO.isconstant(L) && return L
    vals = similar(L.vals)
    L.coeffs(vals, u, p, t)
    return VecSum(L.coeffs, vals, L.ops)
end

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

function SQA.numeric_assemble_td(be::QuantumToolboxBackend, basis, terms)
    functions = Function[c isa Number ? ((p, t) -> c) : SQA.lift_pt(c) for (c, _) in terms]
    function coeffs(out::Vector{ComplexF64}, u, p, t)
        for k in eachindex(functions)
            out[k] = functions[k](p, t)
        end
        return out
    end
    coeffs(u, p, t) = ComplexF64[f(p, t) for f in functions]
    ops = [SO.MatrixOperator(fuse(be, basis, factors).data) for (_, factors) in terms]
    refid = SQA.numeric_identity(be, basis)
    return QTB.QobjEvo(VecSum(coeffs, zeros(ComplexF64, length(ops)), ops), QTB.Operator(), refid.dimensions)
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

function SQA.formal_td_operator(
        terms::Vector{Tuple{SQA.CNum, O}}, ctx::SQA.NumericContext{QuantumToolboxBackend}, time_parameter,
    ) where {O}
    ops = [SO.MatrixOperator(op.data) for (_, op) in terms]
    coeffs = SQA.compile_td_coeffs((c for (c, _) in terms), time_parameter)
    refid = SQA.numeric_identity(ctx.backend, ctx.basis)
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
function SQA.numeric_function(::QuantumToolboxBackend, ::typeof(SQA.expim), op::QTB.QuantumObject)
    data = SQA.numeric_expim(Matrix{ComplexF64}(op.data))
    return QTB.to_sparse(QTB.Qobj(data, QTB.Operator(), op.dimensions))
end

function SQA.numeric_function(::QuantumToolboxBackend, f::F, op::QTB.QuantumObject) where {F}
    return f(op)
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
SQA.to_numeric(q::SQA.QExpr, dims::QTBDims, d::AbstractDict{<:SQA.QSym}) =
    SQA.to_numeric(q, dims; operators = d, adjoint_ops = false)

SQA.to_numeric(q::SQA.QAdd, dims::QTBDims, d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.to_numeric(q, dims; operators = d, adjoint_ops = false)
SQA.to_numeric(x::Number, dims::QTBDims, ::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS) =
    SQA.to_numeric_kw(
    QuantumToolboxBackend(), x, dims;
    parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
    operators = SQA.NO_SUBS, adjoint_ops = false, op_type = SQA.EagerOperator(),
)

function SQA.to_numeric(
        q::SQA.QAdd, dims::QTBDims, sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:SQA.QSym} = SQA.NO_SUBS, scalar_subs::AbstractDict = SQA.NO_SCALAR_SUBS;
        default_time_dependence::Bool = true,
    )
    return SQA.to_numeric_kw(
        QuantumToolboxBackend(), q, dims; parameter = scalar_subs, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = d, adjoint_ops = false, op_type = SQA.EagerOperator(),
        default_time_dependence, sites = Dict{Int, Vector{Int}}(sites),
    )
end

function SQA.to_numeric(
        op::SQA.QField, dims::QTBDims;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
        default_time_dependence::Bool = true,
    )
    return SQA.to_numeric_kw(QuantumToolboxBackend(), op, dims; parameter, time_parameter, operators, adjoint_ops, op_type, default_time_dependence)
end
function SQA.to_numeric(
        x::Union{Number, BasicSymbolic}, dims::QTBDims;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
        default_time_dependence::Bool = true,
    )
    return SQA.to_numeric_kw(QuantumToolboxBackend(), x, dims; parameter, time_parameter, operators, adjoint_ops, op_type, default_time_dependence)
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

SQA.numeric_zero(be::QuantumToolboxBackend, basis) = QTB.to_sparse(0.0im * SQA.numeric_identity(be, basis))

function SQA.numeric_add_term!(result, coefficient::ComplexF64, ops::Vector{Op}, ctx::SQA.NumericContext{QuantumToolboxBackend})
    term = SQA.static_operator(SQA.OpTerm[(1.0 + 0.0im, ops)], ctx, SQA.EagerOperator())
    result.data .+= coefficient .* term.data
    return result
end

function parameter_variables!(variables::Vector{Num}, c::SQA.CNum)
    append!(variables, Num.(SQA.coefficient_variables(SQA.to_num(c))))
    return variables
end

function parameter_variables!(variables::Vector{Num}, q::SQA.QAdd, ctx::SQA.NumericContext, rules::SQA.SubstitutionRules)
    if isempty(q.indices)
        for c in values(q.arguments)
            parameter_variables!(variables, c)
        end
    else
        ctx = SQA.sum_context(q, ctx, rules.scalars)
        SQA.foreach_unrolled(q, ctx.sites, SQA.split_scalar_subs(rules.scalars)) do term, c, substitution
            parameter_variables!(variables, c)
        end
    end
    return variables
end

function parameter_variables!(variables::Vector{Num}, q::SQA.QExpr, ctx::SQA.NumericContext, rules::SQA.SubstitutionRules)
    parameter_variables!(variables, q.coeff)
    if q.kind == SQA.QEXPR_ADD || q.kind == SQA.QEXPR_MUL
        for arg in q.args
            parameter_variables!(variables, arg, ctx, rules)
        end
    end
    return variables
end

SQA.numeric_time_parameters(::Op, ::SQA.NumericContext{QuantumToolboxBackend}, ::SQA.SubstitutionRules, tp, ::Bool) = tp

function SQA.numeric_time_parameters(
        op::Union{SQA.QField, Number, BasicSymbolic}, ctx::SQA.NumericContext{QuantumToolboxBackend},
        rules::SQA.SubstitutionRules, time_parameter, default_time_dependence::Bool,
    )
    q = if op isa SQA.QField
        SQA.substitute_with(op isa SQA.QSym ? SQA.to_qadd(op) : op, rules)
    else
        SQA.single_qadd(SQA.to_cnum(SQA.as_cnum(SQA.apply_parameter(op, rules))), Op[])
    end
    result = isempty(time_parameter) ? Dict{Num, Any}() : copy(time_parameter)
    specified = isempty(result) ? Num[] : first(SQA.time_basis(result))
    for key in parameter_variables!(Num[], q, ctx, rules)
        (haskey(result, key) || any(v -> isequal(v, key), specified)) && continue
        result[key] = let key = key
            default_time_dependence ? ((p, t) -> p[key](t)) : ((p, t) -> p[key])
        end
    end
    return result
end

function inferred_site_count(index::SQA.Index, rules::SQA.SubstitutionRules)
    value = try
        SQA.to_complex(SQA.apply_parameter(SQA.index_range(index), rules))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("dimension inference needs static integer index ranges; supply `parameter` or explicit dimensions"))
    end
    if !iszero(imag(value)) || !isinteger(real(value)) || !(0 < real(value) <= typemax(Int))
        throw(ArgumentError("dimension inference needs positive integer index ranges; supply explicit dimensions"))
    end
    return Int(real(value))
end

function inferred_layout!(
        dimensions::Dict{Int, Int}, counts::Dict{Int, Int}, unindexed::Set{Int},
        q::SQA.QAdd, rules::SQA.SubstitutionRules,
    )
    for index in q.indices
        slot = Int(index.space_index)
        count = inferred_site_count(index, rules)
        get(counts, slot, count) == count || throw(ArgumentError("conflicting indexed ranges; supply explicit dimensions and sites"))
        counts[slot] = count
    end
    for op in SQA.get_operators(q)
        slot = Int(op.space_index)
        dimension = if SQA.is_transition(op)
            Int(op.nlev)
        elseif SQA.is_pauli(op)
            2
        else
            throw(ArgumentError("cannot infer the dimension of `$op`; supply an explicit basis, bosonic cutoff, or spin length"))
        end
        get(dimensions, slot, dimension) == dimension || throw(ArgumentError("conflicting subsystem dimensions; supply an explicit basis"))
        dimensions[slot] = dimension
        if SQA.has_index(op.index)
            if iszero(SQA.index_slot(op.index)) && !(op.index in q.indices)
                throw(ArgumentError("cannot infer a layout for a free index; bind it with Σ or supply explicit dimensions and sites"))
            end
            count = inferred_site_count(op.index, rules)
            get(counts, slot, count) == count || throw(ArgumentError("conflicting indexed ranges; supply explicit dimensions and sites"))
            counts[slot] = count
            SQA.index_slot(op.index) <= count || throw(ArgumentError("indexed site exceeds its range"))
        else
            push!(unindexed, slot)
        end
    end
    return nothing
end

function inferred_layout!(
        dimensions::Dict{Int, Int}, counts::Dict{Int, Int}, unindexed::Set{Int},
        q::SQA.QExpr, rules::SQA.SubstitutionRules,
    )
    for arg in q.args
        inferred_layout!(dimensions, counts, unindexed, arg, rules)
    end
    return nothing
end

function inferred_layout(op::SQA.QField, rules::SQA.SubstitutionRules)
    dimensions, counts, unindexed = Dict{Int, Int}(), Dict{Int, Int}(), Set{Int}()
    q = op isa SQA.QSym ? SQA.to_qadd(op) : op
    inferred_layout!(dimensions, counts, unindexed, q, rules)
    isempty(dimensions) && throw(ArgumentError("no operator dimensions remain to infer; supply an explicit basis"))
    slots = sort!(collect(keys(dimensions)))
    slots == collect(1:last(slots)) || throw(ArgumentError("missing subsystem dimensions; supply an explicit basis"))
    dims, sites = Int[], Dict{Int, Vector{Int}}()
    for slot in slots
        count = get(counts, slot, 1)
        slot in unindexed && count > 1 && throw(ArgumentError("ambiguous indexed and unindexed subsystem; supply explicit dimensions and sites"))
        sites[slot] = collect((length(dims) + 1):(length(dims) + count))
        append!(dims, fill(dimensions[slot], count))
    end
    return (length(dims) == 1 ? only(dims) : Tuple(dims)), sites
end

function SQA.to_numeric(
        op::SQA.QField, backend::QuantumToolboxBackend;
        parameter = SQA.NO_PARAMETER, time_parameter = SQA.NO_TIME_PARAMETER,
        operators = SQA.NO_SUBS, adjoint_ops = true, op_type = SQA.EagerOperator(),
        default_time_dependence::Bool = true,
    )
    rules = SQA.expand_parameter(parameter)
    dims, sites = inferred_layout(op, rules)
    return SQA.to_numeric_kw(
        backend, op, dims; parameter, time_parameter, operators, adjoint_ops, op_type,
        default_time_dependence, sites,
    )
end

end
