module SecondQuantizedAlgebraQuantumOpticsBaseExt

# QuantumOpticsBase backend for `to_numeric`/`numeric_average`. Provides the operator-to-
# matrix map, basis construction, lazy embedding, and vector-backed lazy assembly. Loaded
# automatically by `using QuantumOpticsBase`.

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumOpticsBackend, Op
import QuantumOpticsBase as QOB

const QOBState = Union{QOB.StateVector, QOB.AbstractOperator}

# --- operator -> matrix (closed value-branch ladder + open Val fallthrough) -------------

@inline function SQA.numeric_operator(be::QuantumOpticsBackend, op::Op, b::QOB.Basis)
    k = op.kind
    k === SQA.OP_DESTROY    && return QOB.destroy(b)
    k === SQA.OP_CREATE     && return QOB.create(b)
    k === SQA.OP_POSITION   && return (QOB.destroy(b) + QOB.create(b)) / sqrt(2)
    k === SQA.OP_MOMENTUM   && return im * (QOB.create(b) - QOB.destroy(b)) / sqrt(2)
    k === SQA.OP_TRANSITION && return QOB.transition(b, Int(op.l1), Int(op.l2))
    k === SQA.OP_PAULI      && return _pauli_matrix(op, b)
    k === SQA.OP_SPIN       && return _spin_matrix(op, b)
    return SQA.numeric_operator(be, Val(k), op, b)
end

# Separate, throwing extension point for custom kinds (`::Val` is `::Val{K} where K`, so it
# must NOT share a signature with the main `op::Op` method above).
SQA.numeric_operator(::QuantumOpticsBackend, ::Val{K}, op::Op, b::QOB.Basis) where {K} =
    throw(ArgumentError("no QuantumOpticsBase numeric_operator for role $K on $(typeof(b))"))

function _pauli_matrix(op::Op, b::QOB.SpinBasis)
    b.spinnumber == 1 // 2 ||
        throw(ArgumentError("Pauli operators require SpinBasis(1//2), got SpinBasis($(b.spinnumber))"))
    op.l1 == 1 && return QOB.sigmax(b)
    op.l1 == 2 && return QOB.sigmay(b)
    return QOB.sigmaz(b)
end

function _spin_matrix(op::Op, b::QOB.SpinBasis)
    op.l1 == 1 && return 0.5 * QOB.sigmax(b)
    op.l1 == 2 && return 0.5 * QOB.sigmay(b)
    return 0.5 * QOB.sigmaz(b)
end

# --- basis construction (HilbertSpace form) --------------------------------------------

SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.FockSpace, N) = QOB.FockBasis(Int(N))
SQA.numeric_basis(::QuantumOpticsBackend, h::SQA.NLevelSpace, _) = QOB.NLevelBasis(h.n)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.PauliSpace, _) = QOB.SpinBasis(1 // 2)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.SpinSpace, s) = QOB.SpinBasis(s)
SQA.numeric_basis(::QuantumOpticsBackend, ::SQA.PhaseSpace, N) = QOB.FockBasis(Int(N))
function SQA.numeric_basis(be::QuantumOpticsBackend, h::SQA.ProductSpace, dims)
    subs = h.spaces
    return QOB.tensor(ntuple(i -> SQA.numeric_basis(be, subs[i], dims[i]), length(subs))...)
end

# --- subsystem basis / embedding / identity --------------------------------------------

SQA.numeric_subbasis(::QuantumOpticsBackend, b::QOB.Basis, slot::Int) = b
SQA.numeric_subbasis(::QuantumOpticsBackend, b::QOB.CompositeBasis, slot::Int) = b.bases[slot]

SQA.numeric_embed(::QuantumOpticsBackend, b::QOB.Basis, slot::Int, m) = m
SQA.numeric_embed(::QuantumOpticsBackend, b::QOB.CompositeBasis, slot::Int, m) =
    QOB.LazyTensor(b, [slot], (m,))

SQA.numeric_identity(::QuantumOpticsBackend, b::QOB.Basis) = one(b)
SQA.numeric_identity(::QuantumOpticsBackend, b::QOB.CompositeBasis) =
    QOB.LazyTensor(b, collect(1:length(b.bases)), Tuple(one(bi) for bi in b.bases))

# --- vector-backed lazy assembly -------------------------------------------------------

# Static: one concrete `LazySum` type for any term count. Pass the concrete basis to the
# 5-arg constructor so `BL`/`BR` are pinned and the result is `@inferred`-stable (the 2-arg
# `LazySum(coeffs, ops)` is only runtime-n-stable, not inference-stable).
function SQA.numeric_assemble(be::QuantumOpticsBackend, b, terms)
    coeffs = ComplexF64[]
    ops = QOB.AbstractOperator[]
    for (c, factors) in terms
        push!(coeffs, c)
        push!(ops, _fuse(be, b, factors))
    end
    if isempty(ops)
        push!(coeffs, 0.0im)
        push!(ops, SQA.numeric_identity(be, b))
    end
    return QOB.LazySum(ComplexF64, b, b, coeffs, ops)
end

# Time-dependent: native TimeDependentSum, vector-backed (same n-stability story).
function SQA.numeric_assemble_td(be::QuantumOpticsBackend, b, td_terms)
    coeffs = Function[]
    ops = QOB.AbstractOperator[]
    for (c, factors) in td_terms
        push!(
            coeffs, c isa Function ? c : (
                    let c = c
                        t -> c
                end
                )
        )
        push!(ops, _fuse(be, b, factors))
    end
    return QOB.TimeDependentSum(ComplexF64, b, b, coeffs, ops)
end

# Per-term operator: identity for an empty product, the bare leaf for one factor, a lazy
# product otherwise.
function _fuse(be::QuantumOpticsBackend, b, factors)
    isempty(factors) && return SQA.numeric_identity(be, b)
    length(factors) == 1 && return factors[1]
    return QOB.LazyProduct(factors...)
end

# --- expectation -----------------------------------------------------------------------

SQA.numeric_expect(::QuantumOpticsBackend, numop, state::QOBState) =
    ComplexF64(QOB.expect(numop, state))
SQA.numeric_expect(::QuantumOpticsBackend, numop, states::AbstractVector) =
    QOB.expect(QOB.sparse(numop), states)

end # module
