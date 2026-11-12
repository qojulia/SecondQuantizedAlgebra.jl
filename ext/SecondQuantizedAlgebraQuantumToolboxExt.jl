module SecondQuantizedAlgebraQuantumToolboxExt

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumToolboxBackend, Op
import QuantumToolbox as QTB
import SciMLOperators as SO
import SymbolicUtils: BasicSymbolic
using Symbolics: Num
import SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals, sparse
import LinearAlgebra: dot

const QTBDims = Union{Integer, AbstractVector{<:Integer}, Tuple{Vararg{T}} where {T <: Integer}}

function append_scaled!(rows::Vector{Int}, cols::Vector{Int}, vals::Vector{T}, c::T, m::SparseMatrixCSC) where {T}
    rv, nz = rowvals(m), nonzeros(m)
    for j in axes(m, 2), k in nzrange(m, j)
        push!(rows, rv[k]); push!(cols, j); push!(vals, c * nz[k])
    end
    return rows
end
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
    result = QTB.to_sparse(0.0im * SQA.numeric_identity(be, basis))
    for (c, factors) in terms
        result += QTB.to_sparse(c * fuse(be, basis, factors))
    end
    return QTB.QobjEvo(result)
end

native_coefficient(c::Number) = (p, t) -> c
native_coefficient(f::Function) = SQA.lift_pt(f)

function native_time_operator(coefficients, operators)
    entries = Tuple((op, native_coefficient(c)) for (c, op) in zip(coefficients, operators))
    op = first(operators)
    return QTB.QobjEvo(entries)::QTB.QuantumObjectEvolution{typeof(op.type), typeof(op.dimensions)}
end

function SQA.numeric_assemble_td(be::QuantumToolboxBackend, basis, terms)
    ops = [fuse(be, basis, factors) for (_, factors) in terms]
    if isempty(ops)
        return QTB.QobjEvo(SQA.numeric_identity(be, basis), (p, t) -> 0.0im)
    end
    return native_time_operator((c for (c, _) in terms), ops)
end

function SQA.to_numeric_td(terms::Vector{SQA.CoeffTerm}, ctx::SQA.NumericContext{QuantumToolboxBackend}, time_parameter)
    leaves = SQA.leaf_cache(ctx)
    ops = [fuse(ctx.backend, ctx.basis, SQA.term_factors(factors, ctx, leaves)) for (_, factors) in terms]
    if isempty(ops)
        return QTB.QobjEvo(SQA.numeric_identity(ctx.backend, ctx.basis), (p, t) -> 0.0im)
    end
    basevars, valuefuncs, p_aware = SQA.time_basis(time_parameter)
    coefficients = SQA.TDCoeff[]
    for (c, _) in terms
        push!(coefficients, SQA.td_coeff(SQA.to_num(c), basevars, valuefuncs, p_aware))
    end
    return native_time_operator(coefficients, ops)
end

function SQA.formal_td_operator(
        terms::Vector{Tuple{SQA.CNum, O}}, ctx::SQA.NumericContext{QuantumToolboxBackend}, time_parameter,
    ) where {O}
    basevars, valuefuncs, p_aware = SQA.time_basis(time_parameter)
    coefficients = SQA.TDCoeff[]
    for (c, _) in terms
        push!(coefficients, SQA.td_coeff(SQA.to_num(c), basevars, valuefuncs, p_aware))
    end
    return native_time_operator(coefficients, (op for (_, op) in terms))
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

SQA.numeric_expect(
    be::QuantumToolboxBackend, numop::QTB.QuantumObjectEvolution, state::QTB.QuantumObject,
) = SQA.numeric_expect(be, qtb_eager(numop), state)

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
    dimensions = super_dimensions(ctx.backend, ctx.basis)
    ops = [QTB.Qobj(sandwich_matrix(t.left, t.right, ctx, cache), QTB.SuperOperator(), dimensions) for (t, _) in S.arguments]
    if isempty(ops)
        identity_op = SQA.numeric_identity(ctx.backend, ctx.basis)
        return QTB.QobjEvo(QTB.spre(identity_op), (p, t) -> 0.0im)
    end
    basevars, valuefuncs, p_aware = SQA.time_basis(time_parameter)
    coefficients = SQA.TDCoeff[]
    for c in values(S.arguments)
        push!(coefficients, SQA.td_coeff(SQA.to_num(c), basevars, valuefuncs, p_aware))
    end
    return native_time_operator(coefficients, ops)
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
