"""
    ExactComplex
    BigExactComplex

The two tiers of exact numbers, `GaussianRational{Int}` and `GaussianRational{BigInt}`.

Exactness is a property of the type: an exact value is a `GaussianRational`, and a
`ComplexF64` is always inexact. Every function that can produce an exact scalar takes the
exact type `E` of its tier. Arithmetic in `ExactComplex` is checked and throws
`OverflowError`; `tiered` (`coeff.jl`) catches that and redoes the whole computation in
`BigExactComplex`, so the steps must not mutate their inputs. A result returns to
`ExactComplex` whenever every scalar fits.
"""
const ExactComplex = GaussianRational{Int}
const BigExactComplex = GaussianRational{BigInt}
const ExactScalar = Union{ExactComplex, BigExactComplex}
const CoeffScalar = Union{ComplexF64, ExactComplex, BigExactComplex}
const NativeScalar = Union{ComplexF64, ExactComplex}
const RawExpression = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}

struct NativeSlot
    re::Int
    im::Int
    den::Int
end
const FLOAT_DENOMINATOR = 0

@inline holds_float(s::Union{NativeSlot, ExactScalar}) = s.den == FLOAT_DENOMINATOR
@inline NativeSlot(z::ExactComplex) = NativeSlot(z.re, z.im, z.den)
@inline NativeSlot(z::ComplexF64) =
    NativeSlot(reinterpret(Int, real(z)), reinterpret(Int, imag(z)), FLOAT_DENOMINATOR)
@inline function native_scalar(s::NativeSlot)
    holds_float(s) && return ComplexF64(reinterpret(Float64, s.re), reinterpret(Float64, s.im))
    return unsafe_gaussian(s.re, s.im, s.den)
end

@inline pack_scalar(::Type{E}, z::E) where {E <: ExactScalar} = z
@inline pack_scalar(::Type{GaussianRational{T}}, z::ComplexF64) where {T} = unsafe_gaussian(
    T(reinterpret(Int, real(z))), T(reinterpret(Int, imag(z))), T(FLOAT_DENOMINATOR),
)
@inline function unpack_scalar(s::ExactScalar)
    holds_float(s) || return s
    return ComplexF64(reinterpret(Float64, Int(s.re)), reinterpret(Float64, Int(s.im)))
end

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactScalar) = z

@inline fits_small(x::BigInt) = typemin(Int) < x <= typemax(Int)
@inline fits_small(x::Rational{BigInt}) = fits_small(numerator(x)) && fits_small(denominator(x))
@inline fits_small(z::BigExactComplex) = fits_small(z.re) && fits_small(z.im) && fits_small(z.den)

@inline as_tier(::Type{E}, z::E) where {E <: ExactScalar} = z
@inline as_tier(::Type{BigExactComplex}, z::ExactComplex) = BigExactComplex(z)
@inline as_tier(::Type{<:ExactScalar}, x) = x
@inline function as_tier(::Type{ExactComplex}, z::BigExactComplex)
    fits_small(z) || throw(OverflowError("exact scalar exceeds Int"))
    return unsafe_gaussian(Int(z.re), Int(z.im), Int(z.den))
end

@inline to_float_scalar(z::ComplexF64) = z
@inline to_float_scalar(z::ExactScalar) = ComplexF64(z)

@inline scalar_mul(a::GaussianRational{T}, b::GaussianRational{T}) where {T} = a * b
@inline scalar_mul(a::ComplexF64, b::CoeffScalar) = normalize_scalar(a * to_float_scalar(b))
@inline scalar_mul(a::ExactScalar, b::ComplexF64) = normalize_scalar(to_float_scalar(a) * b)
@inline scalar_add(a::GaussianRational{T}, b::GaussianRational{T}) where {T} = a + b
@inline scalar_add(a::ComplexF64, b::CoeffScalar) = normalize_scalar(a + to_float_scalar(b))
@inline scalar_add(a::ExactScalar, b::ComplexF64) = normalize_scalar(to_float_scalar(a) + b)
