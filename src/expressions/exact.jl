"""
Exact scalar tiers and the radical-folding primitives that keep a numeric radical atom's
exponent in `(0, 1)`. Exactness is a property of the type: an exact value is a
`GaussianRational` (`gaussian.jl`), and a `ComplexF64` is always inexact. Nothing here depends
on `Monomial`, `Poly`, or `Coeff`; those types are defined later (`monomial.jl`, `cnum.jl`) and
depend on this file. See `docs/src/devdocs.md` "Exact coefficients" for the invariants.

Every function that can produce an exact scalar takes the exact type `E` of its tier,
`ExactComplex` or `BigExactComplex`. Arithmetic in `ExactComplex` is checked and throws
`OverflowError`; `tiered` (`cnum.jl`) catches that and redoes the whole computation in
`BigExactComplex`.
"""
const ExactComplex = GaussianRational{Int}
const BigExactComplex = GaussianRational{BigInt}
const ExactScalar = Union{ExactComplex, BigExactComplex}
const CoeffScalar = Union{ComplexF64, ExactComplex, BigExactComplex}
# The inline value of a native coefficient: a small exact value or an inexact float.
const NativeScalar = Union{ComplexF64, ExactComplex}

# A `NativeScalar` stored with a concrete layout: an `ExactComplex` as its fields, or with
# `den == 0`, which no exact value has, the bit patterns of a `ComplexF64`. A union field
# would carry a type selector through every `Coeff` copy and call.
struct NativeSlot
    re::Int
    im::Int
    den::Int
end
@inline NativeSlot(z::ExactComplex) = NativeSlot(z.re, z.im, z.den)
@inline NativeSlot(z::ComplexF64) =
    NativeSlot(reinterpret(Int, real(z)), reinterpret(Int, imag(z)), 0)
@inline function native_scalar(s::NativeSlot)::NativeScalar
    iszero(s.den) && return ComplexF64(reinterpret(Float64, s.re), reinterpret(Float64, s.im))
    return unsafe_gaussian(s.re, s.im, s.den)
end

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactScalar) = z

@inline fits_small(x::BigInt) = typemin(Int) < x <= typemax(Int)
@inline fits_small(x::Rational{BigInt}) = fits_small(numerator(x)) && fits_small(denominator(x))
@inline fits_small(z::BigExactComplex) = fits_small(z.re) && fits_small(z.im) && fits_small(z.den)

# The value of `z` in the exact tier `E`. Narrowing a value that does not fit
# `ExactComplex` throws the same `OverflowError` as small-tier arithmetic, so a computation
# that meets a big operand is redone in the big tier by `tiered`.
@inline as_tier(::Type{E}, z::E) where {E <: ExactScalar} = z
@inline as_tier(::Type{<:ExactScalar}, z::ComplexF64) = z
@inline as_tier(::Type{BigExactComplex}, z::ExactComplex) = BigExactComplex(z)
# Arguments that hold no exact scalar pass through unchanged.
@inline as_tier(::Type{<:ExactScalar}, x) = x
@inline function as_tier(::Type{ExactComplex}, z::BigExactComplex)
    fits_small(z) || throw(OverflowError("exact scalar exceeds Int"))
    return unsafe_gaussian(Int(z.re), Int(z.im), Int(z.den))
end

@inline to_float_scalar(z::ComplexF64) = z
@inline to_float_scalar(z::ExactScalar) = ComplexF64(z)

# Exact with exact stays exact in the tier `E`; anything combined with a float is a float.
# Both operands of an exact operation are already in `E`, so mixing tiers is a `MethodError`.
# The exact methods are spelled per tier: under `where {E}` the analysis of the method on its
# own widens `a` and `b` to the union of tiers independently.
for E in (ExactComplex, BigExactComplex)
    @eval begin
        @inline scalar_mul(::Type{$E}, a::$E, b::$E) = a * b
        @inline scalar_add(::Type{$E}, a::$E, b::$E) = a + b
        @inline scalar_inv(::Type{$E}, z::$E) = inv(z)
    end
end
@inline scalar_mul(::Type{<:ExactScalar}, a::ComplexF64, b::ComplexF64) =
    normalize_scalar(a * b)
@inline scalar_mul(::Type{E}, a::E, b::ComplexF64) where {E <: ExactScalar} =
    normalize_scalar(to_float_scalar(a) * b)
@inline scalar_mul(::Type{E}, a::ComplexF64, b::E) where {E <: ExactScalar} =
    normalize_scalar(a * to_float_scalar(b))
@inline scalar_add(::Type{<:ExactScalar}, a::ComplexF64, b::ComplexF64) =
    normalize_scalar(a + b)
@inline scalar_add(::Type{E}, a::E, b::ComplexF64) where {E <: ExactScalar} =
    normalize_scalar(to_float_scalar(a) + b)
@inline scalar_add(::Type{E}, a::ComplexF64, b::E) where {E <: ExactScalar} =
    normalize_scalar(a + to_float_scalar(b))
@inline scalar_inv(::Type{<:ExactScalar}, z::ComplexF64) = inv(z)

# A radical atom is a `Const{SymReal}` of an `Int`; the type test first keeps the field
# access static for the abstract factor element type.
@inline function is_radical_atom(s::SymbolicUtils.BasicSymbolic)
    s isa SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal} || return false
    return SymbolicUtils.isconst(s) && s.val isa Int
end

radical_power(::Type{E}, p::Int, n::Int) where {E <: ExactScalar} = E(p, 0)^n

@inline function needs_radical_fold(
        syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}},
    )
    @inbounds for i in eachindex(syms)
        is_radical_atom(syms[i]) || continue
        e = exps[i]
        (e <= 0 || e >= 1) && return true
    end
    return false
end

const RADICAL_TRIAL_BOUND = 1 << 16

function prime_factorization!(factors::Vector{Tuple{Int, Int}}, n::Integer)::Bool
    p = 2
    while p < RADICAL_TRIAL_BOUND && p * p <= n
        k = 0
        while n % p == 0
            n = div(n, p); k += 1
        end
        k > 0 && push!(factors, (p, k))
        p += p == 2 ? 1 : 2
    end
    n == 1 && return true
    n < big(RADICAL_TRIAL_BOUND)^2 || return false
    push!(factors, (Int(n), 1))
    return true
end
