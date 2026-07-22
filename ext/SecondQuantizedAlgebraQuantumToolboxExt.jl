module SecondQuantizedAlgebraQuantumToolboxExt

# QuantumToolbox backend for `to_numeric`/`numeric_average`. Per-term operators are eager
# concrete `QuantumObject`s; the laziness lives in a vector-backed `VecSum` (a small custom
# `SciMLOperators.AbstractSciMLOperator`) wrapped in a `QobjEvo`, so the result is one
# concrete type for any term count and for static or time-dependent coefficients. Loaded by
# `using QuantumToolbox` (which also loads SciMLOperators, the second extension trigger).

import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra: QuantumToolboxBackend, Op
import QuantumToolbox as QTB
import SciMLOperators as SO
import SymbolicUtils: BasicSymbolic
# `mul!` is imported from its owner `LinearAlgebra` (a weakdep trigger of this extension);
# `transpose`/`adjoint` come through Base.
import LinearAlgebra: mul!

# QuantumToolbox represents a "basis" as integer dims: an `Int` for a simple space, a
# `Vector{Int}` for a composite one.
const QTBDims = Union{Integer, AbstractVector{<:Integer}}

# =======================================================================================
# Vector-backed lazy sum
#
# SciMLOperators' `AddedOperator` is `Tuple`-locked (type-unstable for a runtime term count).
# `VecSum` uses `Vector`s instead, so `VecSum{ComplexF64}` is ONE concrete type for any term
# count (one dynamic dispatch per term, like QuantumOptics' `LazySum`). Abstract field eltypes
# are deliberate: a `ScalarOperator`'s type encodes its update function.
# =======================================================================================

struct VecSum{T} <: SO.AbstractSciMLOperator{T}
    coeffs::Vector{SO.ScalarOperator}
    ops::Vector{SO.AbstractSciMLOperator}
end
VecSum(coeffs, ops) = VecSum{ComplexF64}(coeffs, ops)

_cv(c) = convert(Number, c)

Base.size(L::VecSum) = size(L.ops[1])
Base.size(L::VecSum, i::Int) = size(L.ops[1], i)
Base.eltype(::VecSum{T}) where {T} = T
SO.issquare(::VecSum) = true
SO.isconstant(L::VecSum) = all(SO.isconstant, L.coeffs)
SO.update_coefficients!(L::VecSum, u, p, t) =
    (foreach(c -> SO.update_coefficients!(c, u, p, t), L.coeffs); L)

function mul!(w::AbstractVecOrMat, L::VecSum, v::AbstractVecOrMat)
    mul!(w, L.ops[1], v, _cv(L.coeffs[1]), false)
    for i in 2:length(L.ops)
        mul!(w, L.ops[i], v, _cv(L.coeffs[i]), true)
    end
    return w
end
function mul!(w::AbstractVecOrMat, L::VecSum, v::AbstractVecOrMat, α, β)
    mul!(w, L.ops[1], v, _cv(L.coeffs[1]) * α, β)
    for i in 2:length(L.ops)
        mul!(w, L.ops[i], v, _cv(L.coeffs[i]) * α, true)
    end
    return w
end

# Allocating apply (needed by `expect`, which routes through `dot(ψ, L, ψ) = dot(ψ, L*ψ)`).
function Base.:*(L::VecSum, v::AbstractVector)
    w = similar(v, promote_type(eltype(L), eltype(v)), size(L, 1))
    return mul!(w, L, v)
end
function Base.:*(L::VecSum, v::AbstractMatrix)
    w = similar(v, promote_type(eltype(L), eltype(v)), size(L, 1), size(v, 2))
    return mul!(w, L, v)
end

# In-place call forms used by the SciML/QuantumToolbox solver stepping.
(L::VecSum)(w, v, u, p, t) = (SO.update_coefficients!(L, u, p, t); mul!(w, L, v))
(L::VecSum)(w, v, u, p, t, α, β) = (SO.update_coefficients!(L, u, p, t); mul!(w, L, v, α, β))

SO.concretize(L::VecSum) = sum(_cv(L.coeffs[i]) * SO.concretize(L.ops[i]) for i in eachindex(L.ops))
Base.convert(::Type{AbstractMatrix}, L::VecSum) = SO.concretize(L)
# `::AbstractVecOrMat` avoids an ambiguity with the SciMLOperators fallback.
SO.cache_operator(L::VecSum, ::AbstractVecOrMat) = L

# transpose/adjoint MUST return a VecSum (not be wrapped in a TransposedOperator), or the
# mesolve Liouvillian (which tensors transpose(H) with identity) fails its cache reshape.
Base.transpose(L::VecSum{T}) where {T} =
    VecSum{T}(L.coeffs, SO.AbstractSciMLOperator[transpose(o) for o in L.ops])
Base.adjoint(L::VecSum{T}) where {T} =
    VecSum{T}(
    SO.AbstractSciMLOperator[adjoint(c) for c in L.coeffs],
    SO.AbstractSciMLOperator[adjoint(o) for o in L.ops],
)

# --- superoperator construction (mesolve / liouvillian path) ---------------------------
#
# Without these, QuantumToolbox's `_spre`/`_spost` wrap a `VecSum` in a generic lazy tensor
# (a `try/catch` `@warn` that breaks Mooncake, and untraversable by AD). Distributing over the
# terms instead keeps the superoperator a `VecSum` of concrete `MatrixOperator`s: AD-traversable
# and warning-free. `_spre`/`_spost` are linear (coefficients pass through); the dissipator's
# `O'O`/sandwich terms need `VecSum * VecSum`. All three return the single concrete `VecSum{T}`.
QTB._spre(A::VecSum{T}) where {T} =
    VecSum{T}(A.coeffs, SO.AbstractSciMLOperator[QTB._spre(o) for o in A.ops])
QTB._spost(B::VecSum{T}) where {T} =
    VecSum{T}(B.coeffs, SO.AbstractSciMLOperator[QTB._spost(o) for o in B.ops])

# Coefficient product as a PLAIN ScalarOperator (SciMLOperators' `*` returns a
# ComposedScalarOperator that does not fit `coeffs::Vector{ScalarOperator}`). Two constants
# multiply eagerly to stay constant (correct `isconstant`/static `concretize`); otherwise a
# closure defers to update time. The non-mutating `update_coefficients` is AD-safe.
function _scalprod(c1::SO.ScalarOperator, c2::SO.ScalarOperator)
    (SO.isconstant(c1) && SO.isconstant(c2)) && return SO.ScalarOperator(_cv(c1) * _cv(c2))
    return SO.ScalarOperator(
        0.0im,
        (a, u, p, t) -> _cv(SO.update_coefficients(c1, u, p, t)) * _cv(SO.update_coefficients(c2, u, p, t)),
    )
end

# Operator product: materialized MatrixOperator fast path (the lazy `*` fallback does not
# arise after `_spre`/`_spost`, which yield MatrixOperators).
_opprod(x::SO.MatrixOperator, y::SO.MatrixOperator) = SO.MatrixOperator(x.A * y.A)
_opprod(x, y) = x * y

function Base.:*(A::VecSum{T}, B::VecSum{T}) where {T}
    n = length(A.ops) * length(B.ops)
    coeffs = Vector{SO.ScalarOperator}(undef, n)
    ops = Vector{SO.AbstractSciMLOperator}(undef, n)
    k = 1
    for i in eachindex(A.ops), j in eachindex(B.ops)
        coeffs[k] = _scalprod(A.coeffs[i], B.coeffs[j])
        ops[k] = _opprod(A.ops[i], B.ops[j])
        k += 1
    end
    return VecSum{T}(coeffs, ops)
end

# =======================================================================================
# operator -> QuantumObject (closed value-branch ladder + open Val fallthrough)
# =======================================================================================

@inline function SQA.numeric_operator(be::QuantumToolboxBackend, op::Op, N::Integer)
    k = op.kind
    k === SQA.OP_DESTROY    && return QTB.destroy(N)
    k === SQA.OP_CREATE     && return QTB.create(N)
    k === SQA.OP_POSITION   && return QTB.position(N)
    k === SQA.OP_MOMENTUM   && return QTB.momentum(N)
    # QuantumToolbox is 0-indexed; SQA `Transition` levels are 1-indexed (as in QuantumOptics).
    k === SQA.OP_TRANSITION && return QTB.projection(N, Int(op.l1) - 1, Int(op.l2) - 1)
    k === SQA.OP_PAULI      && return _qtb_pauli(op, N)
    k === SQA.OP_SPIN       && return QTB.jmat((N - 1) // 2, _axis(op))
    return SQA.numeric_operator(be, Val(k), op, N)
end

SQA.numeric_operator(::QuantumToolboxBackend, ::Val{K}, op::Op, N::Integer) where {K} =
    throw(ArgumentError("no QuantumToolbox numeric_operator for role $K on dim $N"))

function _qtb_pauli(op::Op, N::Integer)
    N == 2 || throw(ArgumentError("Pauli operators require dim 2, got $N"))
    op.l1 == 1 && return QTB.sigmax()
    op.l1 == 2 && return QTB.sigmay()
    return QTB.sigmaz()
end

_axis(op::Op) = op.l1 == 1 ? :x : op.l1 == 2 ? :y : :z

# =======================================================================================
# basis (integer dims) / subsystem / embedding / identity
# =======================================================================================

SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.FockSpace, N) = Int(N) + 1
SQA.numeric_basis(::QuantumToolboxBackend, h::SQA.NLevelSpace, _) = Int(h.n)
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.PauliSpace, _) = 2
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.SpinSpace, s) = Int(2s + 1)
SQA.numeric_basis(::QuantumToolboxBackend, ::SQA.PhaseSpace, N) = Int(N) + 1
function SQA.numeric_basis(be::QuantumToolboxBackend, h::SQA.ProductSpace, dims)
    SQA._check_product_dims(h, dims)
    return Int[SQA.numeric_basis(be, h.spaces[i], dims[i]) for i in eachindex(h.spaces)]
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

SQA.numeric_num_subsystems(::QuantumToolboxBackend, ::Integer) = 1
SQA.numeric_num_subsystems(::QuantumToolboxBackend, ds::AbstractVector) = length(ds)

# Within-term embedding is EAGER concrete (QuantumToolbox warns on lazy tensor); the lazy sum
# is the only lazy layer.
SQA.numeric_embed(::QuantumToolboxBackend, N::Integer, slot::Int, m) = m
function SQA.numeric_embed(::QuantumToolboxBackend, ds::AbstractVector, slot::Int, m)
    length(ds) == 1 && return m
    return QTB.tensor((i == slot ? m : QTB.qeye(ds[i]) for i in eachindex(ds))...)
end

SQA.numeric_identity(::QuantumToolboxBackend, N::Integer) = QTB.qeye(N)
function SQA.numeric_identity(::QuantumToolboxBackend, ds::AbstractVector)
    length(ds) == 1 && return QTB.qeye(ds[1])
    return QTB.tensor((QTB.qeye(d) for d in ds)...)
end

# =======================================================================================
# vector-backed lazy assembly (static + time-dependent share the one VecSum type)
# =======================================================================================

function SQA.numeric_assemble(be::QuantumToolboxBackend, basis, terms)
    refid = SQA.numeric_identity(be, basis)
    coeffs = SO.ScalarOperator[]
    ops = SO.AbstractSciMLOperator[]
    for (c, factors) in terms
        push!(coeffs, SO.ScalarOperator(ComplexF64(c)))
        push!(ops, SO.MatrixOperator(_fuse(be, basis, factors).data))
    end
    if isempty(ops)
        push!(coeffs, SO.ScalarOperator(0.0im))
        push!(ops, SO.MatrixOperator(refid.data))
    end
    return QTB.QobjEvo(VecSum(coeffs, ops), QTB.Operator(), refid.dimensions)
end

function SQA.numeric_assemble_td(be::QuantumToolboxBackend, basis, td_terms)
    refid = SQA.numeric_identity(be, basis)
    coeffs = SO.ScalarOperator[]
    ops = SO.AbstractSciMLOperator[]
    for (c, factors) in td_terms
        push!(coeffs, _td_scalar(c))
        push!(ops, SO.MatrixOperator(_fuse(be, basis, factors).data))
    end
    return QTB.QobjEvo(VecSum(coeffs, ops), QTB.Operator(), refid.dimensions)
end

# SciML update signature is the four-arg `(a, u, p, t) -> c`, not the two-arg `(p, t)`.
_td_scalar(c::Function) = SO.ScalarOperator(0.0im, (a, u, p, t) -> ComplexF64(c(t)))
_td_scalar(c) = SO.ScalarOperator(ComplexF64(c))
# p-aware: thread `p` into the update function so the QobjEvo is differentiable wrt `p`.
_td_scalar(c::SQA.PControlCoeff) = SO.ScalarOperator(0.0im, (a, u, p, t) -> ComplexF64(c.f(p, t)))

function _fuse(be::QuantumToolboxBackend, basis, factors)
    isempty(factors) && return SQA.numeric_identity(be, basis)
    length(factors) == 1 && return factors[1]
    return foldl(*, factors)
end

# =======================================================================================
# materialization (op_type applied once at the top of `to_numeric`)
# =======================================================================================

# The default (`nothing`) has the same semantics as QuantumOptics: return an eager sparse
# backend operator. `op_type = identity` is the explicit opt-in to the lazily assembled
# `QobjEvo`/`VecSum`. Other callables are applied to an eager `QuantumObject`;
# `SciMLOperators.concretize` is special-cased to return its raw matrix as advertised.
SQA.numeric_materialize(::QuantumToolboxBackend, op, ::Nothing) = QTB.to_sparse(_qtb_eager(op))
function SQA.numeric_materialize(::QuantumToolboxBackend, op, op_type)
    op_type === identity && return op
    eager = _qtb_eager(op)
    op_type === SO.concretize && return eager.data
    return op_type(eager)
end

# Reduce the assembled object to an eager `QuantumObject`: concretize the lazy `VecSum`
# behind a `QobjEvo`; a bare scalar path is already an eager `QuantumObject`.
_qtb_eager(op::QTB.QuantumObjectEvolution) = QTB.Qobj(SO.concretize(op.data); dims = op.dimensions)
_qtb_eager(op::QTB.QuantumObject) = op

# =======================================================================================
# expectation + backend resolution + state/dims convenience
# =======================================================================================

SQA.numeric_expect(::QuantumToolboxBackend, numop, state::QTB.QuantumObject) =
    ComplexF64(QTB.expect(numop, state))
SQA.numeric_expect(::QuantumToolboxBackend, numop, states::AbstractVector) =
    ComplexF64[ComplexF64(QTB.expect(numop, s)) for s in states]

SQA.numeric_backend(::QTB.AbstractQuantumObject) = QuantumToolboxBackend()
function SQA.numeric_basis(o::QTB.AbstractQuantumObject)
    ds = collect(Int, o.dims[1])
    return length(ds) == 1 ? ds[1] : ds
end

# Positional dims form (the QuantumToolbox analog of `to_numeric(op, basis, d)`).
SQA.to_numeric(op::Op, dims::QTBDims, d::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS) =
    SQA.numeric_materialize(
    QuantumToolboxBackend(),
    SQA._to_numeric_static(
        SQA._single_qadd(SQA._CNUM_ONE, Op[op]),
        SQA.NumericContext(QuantumToolboxBackend(), dims, d),
    ),
    nothing,
)
SQA.to_numeric(q::SQA.QAdd, dims::QTBDims, d::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS) =
    SQA.numeric_materialize(
    QuantumToolboxBackend(),
    SQA._to_numeric_static(q, SQA.NumericContext(QuantumToolboxBackend(), dims, d)),
    nothing,
)
SQA.to_numeric(x::Number, dims::QTBDims, ::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS) =
    SQA.numeric_materialize(
    QuantumToolboxBackend(),
    SQA._to_complex(x) * SQA.numeric_identity(QuantumToolboxBackend(), dims),
    nothing,
)

# Indexed positional dims form (the QTB analog of `to_numeric(q, basis, sites, …)`).
SQA.to_numeric(
    q::SQA.QAdd, dims::QTBDims,
    sites::AbstractDict{Int, Vector{Int}},
    d::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS,
    scalar_subs::AbstractDict = SQA._NO_SCALAR_SUBS,
) = SQA._to_numeric_indexed(QuantumToolboxBackend(), dims, q, sites, d, scalar_subs)

# Keyword dims form.
function SQA.to_numeric(
        op::SQA.QField, dims::QTBDims;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{SQA.QSym, Any}(), adjoint_ops = true, op_type = nothing,
    )
    return SQA._to_numeric_kw(QuantumToolboxBackend(), op, dims; parameter, time_parameter, operators, adjoint_ops, op_type)
end
function SQA.to_numeric(
        x::Union{Number, BasicSymbolic}, dims::QTBDims;
        parameter = Dict(), time_parameter = Dict(),
        operators = Dict{SQA.QSym, Any}(), adjoint_ops = true, op_type = nothing,
    )
    return SQA._to_numeric_kw(QuantumToolboxBackend(), x, dims; parameter, time_parameter, operators, adjoint_ops, op_type)
end

# Vector of operators on the direct QTB-dims form (mirrors the `Basis` method in api.jl).
SQA.to_numeric(ops::AbstractVector, dims::QTBDims; kwargs...) =
    [SQA.to_numeric(op, dims; kwargs...) for op in ops]

# State convenience (dims derived from the state).
SQA.to_numeric(op, state::QTB.QuantumObject; kwargs...) = SQA.to_numeric(op, SQA.numeric_basis(state); kwargs...)
SQA.to_numeric(op::SQA.QField, state::QTB.QuantumObject, d::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS) =
    SQA.to_numeric(op, SQA.numeric_basis(state), d)
SQA.to_numeric(x::Number, state::QTB.QuantumObject) = SQA.to_numeric(x, SQA.numeric_basis(state))

SQA.numeric_average(op, state::QTB.QuantumObject, d::AbstractDict{<:SQA.QSym} = SQA._NO_SUBS) =
    SQA._numeric_average(op, state, d)

end # module
