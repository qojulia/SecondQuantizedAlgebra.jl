module SecondQuantizedAlgebraQuantumOpticsBaseExt

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumOpticsBackend, Op
import QuantumOpticsBase as QOB
import SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals, sparse, spdiagm
import LinearAlgebra: dot

const QOBState = Union{QOB.StateVector, QOB.AbstractOperator}

function SQA.numeric_operator(::QuantumOpticsBackend, op::Op, b::QOB.FockBasis)
    k = op.kind
    k === SQA.OP_DESTROY  && return QOB.destroy(b)
    k === SQA.OP_CREATE   && return QOB.create(b)
    k === SQA.OP_POSITION && return (QOB.destroy(b) + QOB.create(b)) / sqrt(2)
    k === SQA.OP_MOMENTUM && return im * (QOB.create(b) - QOB.destroy(b)) / sqrt(2)
    throw(ArgumentError("no QuantumOpticsBase numeric_operator for role $k on $(typeof(b))"))
end

function SQA.numeric_operator(::QuantumOpticsBackend, op::Op, b::QOB.NLevelBasis)
    SQA.is_transition(op) ||
        throw(ArgumentError("no QuantumOpticsBase numeric_operator for role $(op.kind) on $(typeof(b))"))
    return QOB.transition(b, Int(op.l1), Int(op.l2))
end

function SQA.numeric_operator(::QuantumOpticsBackend, op::Op, b::QOB.SpinBasis)
    k = op.kind
    k === SQA.OP_PAULI && return pauli_matrix(op, b)
    k === SQA.OP_SPIN  && return spin_matrix(op, b)
    throw(ArgumentError("no QuantumOpticsBase numeric_operator for role $k on $(typeof(b))"))
end

SQA.numeric_operator(::QuantumOpticsBackend, op::Op, b::QOB.Basis) =
    throw(ArgumentError("no QuantumOpticsBase numeric_operator for role $(op.kind) on $(typeof(b))"))

function SQA.numeric_operator(::QuantumOpticsBackend, op::Op, b::QOB.ManyBodyBasis)
    SQA.is_collective_transition(op) ||
        throw(ArgumentError("Op kind $(op.kind) does not act on a ManyBodyBasis"))
    onebody = b.onebodybasis
    onebody isa QOB.NLevelBasis || throw(
        ArgumentError(
            "CollectiveTransition requires a ManyBodyBasis with an NLevelBasis one-body basis; got $(typeof(onebody))",
        ),
    )
    onebody_op = QOB.transition(onebody, Int(op.l1), Int(op.l2))
    return QOB.manybodyoperator(b, onebody_op)
end

function pauli_matrix(op::Op, b::QOB.SpinBasis)
    b.spinnumber == 1 // 2 ||
        throw(ArgumentError("Pauli operators require SpinBasis(1//2), got SpinBasis($(b.spinnumber))"))
    op.l1 == 1 && return QOB.sigmax(b)
    op.l1 == 2 && return QOB.sigmay(b)
    return QOB.sigmaz(b)
end

function spin_matrix(op::Op, b::QOB.SpinBasis)
    op.l1 == 1 && return 0.5 * QOB.sigmax(b)
    op.l1 == 2 && return 0.5 * QOB.sigmay(b)
    return 0.5 * QOB.sigmaz(b)
end

SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.FockSpace, N) = QOB.FockBasis(Int(N))
SQA.numeric_basis(::QuantumOpticsBackend, h::SQA.NLevelSpace, _) = QOB.NLevelBasis(h.n)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.PauliSpace, _) = QOB.SpinBasis(1 // 2)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.SpinSpace, s) = QOB.SpinBasis(s)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.PhaseSpace, N) = QOB.FockBasis(Int(N))
function SQA.numeric_basis(be::QuantumOpticsBackend, h::SQA.ProductSpace, dims)
    SQA.check_product_dims(h, dims)
    subs = h.spaces
    return QOB.tensor(ntuple(i -> SQA.numeric_basis(be, subs[i], dims[i]), length(subs))...)
end

function SQA.numeric_subbasis(::QuantumOpticsBackend, b::QOB.Basis, slot::Int)
    slot == 1 || throw(ArgumentError("simple QuantumOptics basis has no subsystem slot $slot"))
    return b
end
function SQA.numeric_subbasis(::QuantumOpticsBackend, b::QOB.CompositeBasis, slot::Int)
    1 <= slot <= length(b.bases) || throw(
        ArgumentError(
            "QuantumOptics composite basis has $(length(b.bases)) subsystems, not slot $slot",
        ),
    )
    return b.bases[slot]
end

SQA.numeric_num_subsystems(::QuantumOpticsBackend, ::QOB.Basis) = 1
SQA.numeric_num_subsystems(::QuantumOpticsBackend, b::QOB.CompositeBasis) = length(b.bases)

SQA.numeric_embed(::QuantumOpticsBackend, b::QOB.Basis, slot::Int, m) = m
SQA.numeric_embed(::QuantumOpticsBackend, b::QOB.CompositeBasis, slot::Int, m) =
    QOB.LazyTensor(b, [slot], (m,))

SQA.numeric_identity(::QuantumOpticsBackend, b::QOB.Basis) = one(b)
SQA.numeric_identity(::QuantumOpticsBackend, b::QOB.CompositeBasis) =
    QOB.LazyTensor(b, collect(1:length(b.bases)), map(one, b.bases))

function SQA.numeric_assemble(be::QuantumOpticsBackend, b, terms)
    coeffs = ComplexF64[c for (c, _) in terms]
    return QOB.LazySum(ComplexF64, b, b, coeffs, term_operators(be, b, terms))
end

term_operators(be::QuantumOpticsBackend, b::QOB.Basis, terms) =
    [fuse(be, b, factors) for (_, factors) in terms]
term_operators(be::QuantumOpticsBackend, b::QOB.CompositeBasis, terms) =
    QOB.AbstractOperator[fuse(be, b, factors) for (_, factors) in terms]

function SQA.numeric_assemble_td(be::QuantumOpticsBackend, b, td_terms)
    coeffs = Function[td_function(c) for (c, _) in td_terms]
    ops = [QOB.sparse(fuse(be, b, factors)) for (_, factors) in td_terms]
    return QOB.TimeDependentSum(ComplexF64, b, b, coeffs, ops)
end

td_function(c::ComplexF64) = t -> c
function td_function(@nospecialize(c::Function))
    SQA.tp_reads_p(c) && throw(
        ArgumentError(
            "solver-parameter (p, t) coefficients require the QuantumToolbox backend; " *
                "QuantumOptics time-dependent coefficients are functions of t only",
        ),
    )
    return c
end

function fuse(be::QuantumOpticsBackend, b, factors)
    isempty(factors) && return term_identity(be, b)
    length(factors) == 1 && return factors[1]
    return foldl(*, factors)
end

term_identity(::QuantumOpticsBackend, b::QOB.Basis) = QOB.sparse(one(b))
term_identity(be::QuantumOpticsBackend, b::QOB.CompositeBasis) = SQA.numeric_identity(be, b)

function SQA.numeric_materialize(
        ::QuantumOpticsBackend,
        op::QOB.LazySum{BL, BR, Vector{ComplexF64}, Vector{O}},
        ::SQA.EagerOperator,
    ) where {BL, BR, O}
    rows, cols, vals = Int[], Int[], ComplexF64[]
    for (factor, term) in zip(op.factors, op.operators)
        append_scaled!(rows, cols, vals, factor, QOB.sparse(term).data)
    end
    data = sparse(rows, cols, vals, length(op.basis_l), length(op.basis_r))
    return QOB.Operator(op.basis_l, op.basis_r, data)
end

function append_scaled!(
        rows::Vector{Int}, cols::Vector{Int}, vals::Vector{ComplexF64}, c::ComplexF64,
        m::SparseMatrixCSC,
    )
    rv, nz = rowvals(m), nonzeros(m)
    for j in axes(m, 2), k in nzrange(m, j)
        push!(rows, rv[k]); push!(cols, j); push!(vals, c * nz[k])
    end
    return rows
end
SQA.numeric_materialize(::QuantumOpticsBackend, op, ::SQA.EagerOperator) = QOB.sparse(op)

function SQA.static_operator(
        terms::Vector{SQA.OpTerm}, ctx::SQA.NumericContext{QuantumOpticsBackend, B},
        ::SQA.EagerOperator,
    ) where {B <: QOB.CompositeBasis}
    b = ctx.basis
    factors = Dict{Tuple{Op, Int}, SparseMatrixCSC{ComplexF64, Int}}()
    rows, cols, vals = Int[], Int[], ComplexF64[]
    for (c, ops) in terms
        append_scaled!(rows, cols, vals, c, eager_term(ops, ctx, factors))
    end
    return QOB.Operator(b, b, sparse(rows, cols, vals, length(b), length(b)))
end

function eager_term(ops::Vector{Op}, ctx::SQA.NumericContext, factors::Dict{Tuple{Op, Int}, SparseMatrixCSC{ComplexF64, Int}})
    isempty(ops) && return sparse_identity(length(ctx.basis))
    m = eager_factor(ops[1], ctx, factors)
    for k in 2:length(ops)
        m = m * eager_factor(ops[k], ctx, factors)
    end
    return m
end

function eager_factor(op::Op, ctx::SQA.NumericContext, factors::Dict{Tuple{Op, Int}, SparseMatrixCSC{ComplexF64, Int}})
    haskey(ctx.op_subs, op) && return QOB.sparse(ctx.op_subs[op]).data
    slot = SQA.resolve_slot(op, ctx.sites)
    return get!(() -> embed_matrix(ctx.basis, slot, slot_matrix(op, ctx.basis.bases, slot, 1)), factors, (op, slot))
end

slot_matrix(::Op, ::Tuple{}, slot::Int, k::Int) =
    throw(ArgumentError("QuantumOptics composite basis has $(k - 1) subsystems, not slot $slot"))
function slot_matrix(op::Op, bases::Tuple, slot::Int, k::Int)
    k == slot && return QOB.sparse(SQA.numeric_operator(QuantumOpticsBackend(), op, first(bases))).data
    return slot_matrix(op, Base.tail(bases), slot, k + 1)
end

sparse_identity(n::Int) = spdiagm(n, n, 0 => ones(ComplexF64, n))

function embed_matrix(b::QOB.CompositeBasis, slot::Int, m::SparseMatrixCSC{ComplexF64, Int})
    before = prod(b.shape[1:(slot - 1)]; init = 1)
    after = prod(b.shape[(slot + 1):end]; init = 1)
    return kron(sparse_identity(after), kron(m, sparse_identity(before)))
end
SQA.numeric_materialize(::QuantumOpticsBackend, op, op_type) = op_type(op)

function SQA.numeric_average_impl(q::Union{Op, SQA.QAdd}, state::QOB.LazyKet, d::Dict{SQA.QSym, Union{}})
    ctx = SQA.NumericContext(QuantumOpticsBackend(), state.basis, d)
    kets = Vector{ComplexF64}[k.data for k in state.kets]
    acc = 0.0im
    for (c, ops) in average_terms(q)
        acc += c * product_expect(ops, ctx, kets)
    end
    return acc
end

average_terms(op::Op) = SQA.OpTerm[(1.0 + 0.0im, Op[op])]
average_terms(q::SQA.QAdd) = SQA.op_terms(q)

function product_expect(ops::Vector{Op}, ctx::SQA.NumericContext, kets::Vector{Vector{ComplexF64}})
    result = 1.0 + 0.0im
    for (slot, ψ) in enumerate(kets)
        v = ψ
        for op in Iterators.reverse(ops)
            SQA.resolve_slot(op, ctx.sites) == slot || continue
            v = slot_matrix(op, ctx.basis.bases, slot, 1) * v
        end
        result *= dot(ψ, v)
    end
    return result
end

SQA.operator_adjoint(::QuantumOpticsBackend, v::QOB.Operator) =
    QOB.Operator(v.basis_r, v.basis_l, copy(adjoint(v.data)))

SQA.numeric_backend(::QOBState) = QuantumOpticsBackend()

SQA.numeric_expect(::QuantumOpticsBackend, numop, state::QOBState) =
    ComplexF64(QOB.expect(numop, state))
SQA.numeric_expect(::QuantumOpticsBackend, numop, states::AbstractVector) =
    QOB.expect(QOB.sparse(numop), states)

function SQA.numeric_superoperator(::QuantumOpticsBackend, ctx::SQA.NumericContext, terms::Vector{SQA.SandwichTerm})
    b = ctx.basis
    n = length(b)
    cache = Dict{SQA.SideKey, SparseMatrixCSC{ComplexF64, Int}}()
    rows, cols, vals = Int[], Int[], ComplexF64[]
    for (c, left, right) in terms
        iszero(c) && continue
        m = kron(transpose(SQA.side_matrix(right, ctx, cache)), SQA.side_matrix(left, ctx, cache))
        append_scaled!(rows, cols, vals, c, m)
    end
    return QOB.SuperOperator((b, b), (b, b), sparse(rows, cols, vals, n^2, n^2))
end

SQA.numeric_matrix(::QuantumOpticsBackend, op::QOB.AbstractOperator) = QOB.sparse(op).data

# QuantumOpticsBase exports these names; forwarding keeps `spre(a)` symbolic when both are loaded.
QOB.spre(A::SQA.QField) = SQA.spre(A)
QOB.spost(B::SQA.QField) = SQA.spost(B)
QOB.sprepost(A::SQA.QField, B::SQA.QField) = SQA.sprepost(A, B)
QOB.sprepost(A::SQA.QField, B::Number) = SQA.sprepost(A, B)
QOB.sprepost(A::Number, B::SQA.QField) = SQA.sprepost(A, B)
QOB.liouvillian(H::SQA.QField, J::AbstractVector{<:SQA.QField}; rates = nothing) =
    SQA.liouvillian(H, J; rates)

end
