module ExactStrategyPrototype

using LinearAlgebra
import SecondQuantizedAlgebra as SQA

export ExactCertificate, ExactVerified, NoScalarLift, ScalarLift,
    DiagonalPhaseBlock, RotationBlock, SqueezeBlock, NormalizedInvolutionBlock,
    ScaledInvolutionBlock, NilpotentBlock, UserMapBlock, UnsupportedBlock,
    exact_certificate, body_velocity, split_nambu_map, site_interleaved_map

abstract type AbstractExactBlock end
abstract type ExactStrategy end

struct DiagonalPhaseStrategy <: ExactStrategy end
struct RotationStrategy <: ExactStrategy end
struct SqueezeStrategy <: ExactStrategy end
struct InvolutionStrategy{S} <: ExactStrategy end
struct NilpotentStrategy{K} <: ExactStrategy end
struct UserMapStrategy <: ExactStrategy end

struct ExactVerified
    checks::UInt8
end

struct NoScalarLift end

struct ScalarLift{T}
    value::T
end

"""
Short-lived proof object consumed by rule materialization. The coordinate basis stays in the
closed-adjoint problem and the fields below do not survive in `UnitaryTransform`.
"""
struct ExactCertificate{S <: ExactStrategy, M <: AbstractMatrix, V, L}
    strategy::S
    forward::M
    inverse::M
    verification::V
    scalar_lift::L
end

struct DiagonalPhaseBlock{V <: AbstractVector, L} <: AbstractExactBlock
    phases::V
    inverse_phases::V
    scalar_lift::L
end

struct RotationBlock{T, L} <: AbstractExactBlock
    cosine::T
    sine::T
    scalar_lift::L
end

struct SqueezeBlock{T, L} <: AbstractExactBlock
    cosineh::T
    sineh::T
    scalar_lift::L
end

struct NormalizedInvolutionBlock{S, M <: AbstractMatrix, T, L} <: AbstractExactBlock
    action::M
    cosine::T
    sine::T
    scalar_lift::L
end

struct ScaledInvolutionBlock{M <: AbstractMatrix, T} <: AbstractExactBlock
    action::M
    square::T
end

struct NilpotentBlock{K, M <: AbstractMatrix, T, L} <: AbstractExactBlock
    action::M
    parameter::T
    scalar_lift::L
end

struct UserMapBlock{M <: AbstractMatrix, L} <: AbstractExactBlock
    forward::M
    inverse::M
    scalar_lift::L
end

struct UnsupportedBlock{M <: AbstractMatrix} <: AbstractExactBlock
    action::M
end

DiagonalPhaseBlock(phases, inverse_phases) =
    DiagonalPhaseBlock(phases, inverse_phases, NoScalarLift())
RotationBlock(c, s) = RotationBlock(c, s, NoScalarLift())
SqueezeBlock(c, s) = SqueezeBlock(c, s, NoScalarLift())
NormalizedInvolutionBlock(::Val{S}, action, c, s, lift = NoScalarLift()) where {S} =
    NormalizedInvolutionBlock{S, typeof(action), typeof(c), typeof(lift)}(action, c, s, lift)
NilpotentBlock(::Val{K}, action, parameter, lift = NoScalarLift()) where {K} =
    NilpotentBlock{K, typeof(action), typeof(parameter), typeof(lift)}(
        action, parameter, lift,
    )
UserMapBlock(forward, inverse) = UserMapBlock(forward, inverse, NoScalarLift())

@noinline _refuse(message::AbstractString) = throw(ArgumentError(message))

_scalar_zero(x) = zero(x)
_scalar_one(x) = one(x)
_scalar_add(x, y) = x + y
_scalar_mul(x, y) = x * y
_scalar_neg(x) = -x
_scalar_div(x, y) = x / y

_scalar_zero(::SQA.Coeff) = SQA._CNUM_ZERO
_scalar_one(::SQA.Coeff) = SQA._CNUM_ONE
_scalar_add(x::SQA.Coeff, y::SQA.Coeff) = SQA._add_cnum(x, y)
_scalar_mul(x::SQA.Coeff, y::SQA.Coeff) = SQA._mul_cnum(x, y)
_scalar_mul(x::Integer, y::SQA.Coeff) = SQA._mul_cnum(SQA._to_cnum(x), y)
_scalar_mul(x::SQA.Coeff, y::Integer) = SQA._mul_cnum(x, SQA._to_cnum(y))
_scalar_neg(x::SQA.Coeff) = SQA._neg_cnum(x)
_scalar_div(x::SQA.Coeff, y::Integer) = x / y

function _zero_matrix(example, n::Int, m::Int)
    out = Matrix{typeof(example)}(undef, n, m)
    fill!(out, _scalar_zero(example))
    return out
end

function _identity(example, n::Int)
    out = _zero_matrix(example, n, n)
    @inbounds for i in 1:n
        out[i, i] = _scalar_one(example)
    end
    return out
end


function _matmul(left::AbstractMatrix, right::AbstractMatrix)
    size(left, 2) == size(right, 1) || throw(DimensionMismatch())
    example = _scalar_mul(left[1, 1], right[1, 1])
    out = _zero_matrix(example, size(left, 1), size(right, 2))
    @inbounds for j in axes(right, 2), k in axes(left, 2), i in axes(left, 1)
        out[i, j] = _scalar_add(out[i, j], _scalar_mul(left[i, k], right[k, j]))
    end
    return out
end

function _lincomb(a, left::AbstractMatrix, b, right::AbstractMatrix)
    size(left) == size(right) || throw(DimensionMismatch())
    example = _scalar_add(_scalar_mul(a, left[1, 1]), _scalar_mul(b, right[1, 1]))
    out = Matrix{typeof(example)}(undef, size(left))
    @inbounds for i in eachindex(out)
        out[i] = _scalar_add(_scalar_mul(a, left[i]), _scalar_mul(b, right[i]))
    end
    return out
end

function _verify_inverse(forward::AbstractMatrix{T}, inverse::AbstractMatrix{T}) where {T}
    size(forward, 1) == size(forward, 2) ||
        _refuse("exact strategy requires a square coefficient map")
    size(inverse) == size(forward) ||
        _refuse("exact inverse has size $(size(inverse)); expected $(size(forward))")
    identity = _identity(forward[1, 1], size(forward, 1))
    _matmul(forward, inverse) == identity ||
        _refuse("candidate exact map failed the forward-times-inverse identity")
    _matmul(inverse, forward) == identity ||
        _refuse("candidate exact map failed the inverse-times-forward identity")
    return ExactVerified(0x03)
end

function _certificate(strategy::S, forward::M, inverse::M, lift::L) where {
        S <: ExactStrategy, M <: AbstractMatrix, L,
    }
    verification = _verify_inverse(forward, inverse)
    return ExactCertificate(strategy, forward, inverse, verification, lift)
end

function exact_certificate(block::DiagonalPhaseBlock)
    length(block.phases) == length(block.inverse_phases) ||
        _refuse("diagonal phase map and exact inverse have different lengths")
    forward = Matrix(Diagonal(block.phases))
    inverse = Matrix(Diagonal(block.inverse_phases))
    return _certificate(DiagonalPhaseStrategy(), forward, inverse, block.scalar_lift)
end

function exact_certificate(block::RotationBlock)
    c, s = block.cosine, block.sine
    negative_s = _scalar_neg(s)
    forward = [c negative_s; s c]
    inverse = [c s; negative_s c]
    return _certificate(RotationStrategy(), forward, inverse, block.scalar_lift)
end

function exact_certificate(block::SqueezeBlock)
    c, s = block.cosineh, block.sineh
    negative_s = _scalar_neg(s)
    forward = [c s; s c]
    inverse = [c negative_s; negative_s c]
    return _certificate(SqueezeStrategy(), forward, inverse, block.scalar_lift)
end

function exact_certificate(block::NormalizedInvolutionBlock{S}) where {S}
    S === 1 || S === -1 ||
        _refuse("normalized involution square must be +I or -I, got $S")
    action = block.action
    n, m = size(action)
    n == m || _refuse("normalized involution action must be square")
    identity = _identity(action[1, 1], n)
    _matmul(action, action) == _lincomb(S, identity, 0, identity) ||
        _refuse("declared normalized involution does not square to $(S)I")
    c, s = block.cosine, block.sine
    cc = _scalar_mul(c, c)
    ss = _scalar_mul(s, s)
    relation = S === 1 ? _scalar_add(cc, _scalar_neg(ss)) : _scalar_add(cc, ss)
    relation == _scalar_one(relation) ||
        _refuse("involution coefficients do not satisfy their exact normalization")
    forward = _lincomb(c, identity, s, action)
    inverse = _lincomb(c, identity, _scalar_neg(s), action)
    return _certificate(InvolutionStrategy{S}(), forward, inverse, block.scalar_lift)
end

# This is intentionally unsupported: normalizing by sqrt(square) would require a branch.
exact_certificate(::ScaledInvolutionBlock) = nothing

function exact_certificate(block::NilpotentBlock{K}) where {K}
    K >= 2 || _refuse("nilpotent order must be at least 2, got $K")
    action = block.action
    n, m = size(action)
    n == m || _refuse("nilpotent action must be square")
    typeof(block.parameter) === eltype(action) ||
        _refuse("nilpotent parameter and action coefficients must share one concrete type")
    A = copy(action)
    identity = _identity(A[1, 1], n)
    power = copy(identity)
    forward = copy(identity)
    inverse = copy(identity)
    factorial = 1
    positive_parameter = block.parameter
    negative_parameter = _scalar_neg(positive_parameter)
    positive_power = _scalar_one(positive_parameter)
    negative_power = _scalar_one(negative_parameter)
    for degree in 1:(K - 1)
        power = _matmul(power, A)
        factorial *= degree
        positive_power = _scalar_mul(positive_power, positive_parameter)
        negative_power = _scalar_mul(negative_power, negative_parameter)
        forward = _lincomb(
            _scalar_one(positive_power), forward,
            _scalar_div(positive_power, factorial), power,
        )
        inverse = _lincomb(
            _scalar_one(negative_power), inverse,
            _scalar_div(negative_power, factorial), power,
        )
    end
    _matmul(power, A) == _zero_matrix(A[1, 1], n, n) ||
        _refuse("declared nilpotent action is not zero at order $K")
    return _certificate(NilpotentStrategy{K}(), forward, inverse, block.scalar_lift)
end

function exact_certificate(block::UserMapBlock)
    return _certificate(UserMapStrategy(), block.forward, block.inverse, block.scalar_lift)
end

exact_certificate(::UnsupportedBlock) = nothing

"""Nonscalar gauge data come from body velocities; the scalar lift is not inferred here."""
body_velocity(certificate::ExactCertificate, derivative::AbstractMatrix) =
    _matmul(certificate.inverse, derivative)

"""
Temporarily expose a split-Nambu view from canonical site-interleaved coordinates. The
certificate is converted back before materialization, so ordering metadata is not persisted.
"""
split_nambu_map(map::AbstractMatrix, permutation::AbstractVector{<:Integer}) =
    map[permutation, permutation]

function site_interleaved_map(map::AbstractMatrix, permutation::AbstractVector{<:Integer})
    inverse_permutation = invperm(permutation)
    return map[inverse_permutation, inverse_permutation]
end

end

