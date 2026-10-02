"""
Exact scalar tier: Gaussian-rational arithmetic, the native/exact classification of
`ComplexF64`, and the radical-folding primitives that keep a numeric radical atom's
exponent in `(0, 1)`. Nothing here depends on `Monomial`, `Poly`, or `Coeff`. Those types
are defined later (`monomial.jl`, `cnum.jl`) and depend on this file, which is why it is
included first. See `docs/src/devdocs.md` "Exact coefficients" for the invariants.

Every function that can produce an exact scalar takes the exact type `E` of its tier,
`ExactComplex` or `BigExactComplex`. Arithmetic in `ExactComplex` is Base's checked
`Rational{Int}` arithmetic and throws `OverflowError`; `tiered` (`cnum.jl`) catches that
and redoes the whole computation in `BigExactComplex`.
"""
const ExactComplex = Complex{Rational{Int}}
const BigExactComplex = Complex{Rational{BigInt}}
const ExactScalar = Union{ExactComplex, BigExactComplex}
const CoeffScalar = Union{ComplexF64, ExactComplex, BigExactComplex}

const MAX_EXACT_FLOAT = maxintfloat(Float64)

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactScalar) = z

@inline fits_small(x::Rational{BigInt}) =
    typemin(Int) <= numerator(x) <= typemax(Int) && denominator(x) <= typemax(Int)
@inline fits_small(z::BigExactComplex) = fits_small(real(z)) && fits_small(imag(z))

# The value of `z` in the exact tier `E`. Narrowing a value that does not fit
# `ExactComplex` throws the same `OverflowError` as small-tier arithmetic, so a computation
# that meets a big operand is redone in the big tier by `tiered`.
@inline as_tier(::Type{E}, z::E) where {E <: ExactScalar} = z
@inline as_tier(::Type{<:ExactScalar}, z::ComplexF64) = z
@inline as_tier(::Type{BigExactComplex}, z::ExactComplex) = BigExactComplex(z)
# Arguments that hold no exact scalar pass through unchanged.
@inline as_tier(::Type{<:ExactScalar}, x) = x
@inline function as_tier(::Type{ExactComplex}, z::BigExactComplex)
    fits_small(z) || throw(OverflowError("exact scalar exceeds Rational{Int}"))
    return ExactComplex(Rational{Int}(real(z)), Rational{Int}(imag(z)))
end

# A `ComplexF64` holds an exact value when both parts are integers within `2^53`.
@inline function is_exact_float(z::ComplexF64)
    re, im = real(z), imag(z)
    return abs(re) <= MAX_EXACT_FLOAT && abs(im) <= MAX_EXACT_FLOAT &&
        isinteger(re) && isinteger(im)
end
@inline exact_integer(::Type{E}, z::ComplexF64) where {E <: ExactScalar} =
    E(Int(real(z)), Int(imag(z)))

@inline to_float_scalar(z::ComplexF64) = z
@inline to_float_scalar(z::ExactScalar) = ComplexF64(Float64(real(z)), Float64(imag(z)))

@inline native_product_exact(a::ComplexF64, b::ComplexF64) =
    (abs(real(a)) + abs(imag(a))) * (abs(real(b)) + abs(imag(b))) <= MAX_EXACT_FLOAT
@inline native_sum_exact(s::ComplexF64) =
    abs(real(s)) < MAX_EXACT_FLOAT && abs(imag(s)) < MAX_EXACT_FLOAT

@inline function scalar_mul(::Type{E}, a::ComplexF64, b::ComplexF64) where {E}
    native_product_exact(a, b) && return normalize_scalar(a * b)
    (is_exact_float(a) && is_exact_float(b)) || return normalize_scalar(a * b)
    return exact_integer(E, a) * exact_integer(E, b)
end
@inline scalar_mul(::Type{E}, a::E, b::E) where {E <: ExactScalar} = a * b
@inline function scalar_mul(::Type{E}, a::E, b::ComplexF64) where {E <: ExactScalar}
    is_exact_float(b) && return a * exact_integer(E, b)
    return normalize_scalar(to_float_scalar(a) * b)
end
@inline scalar_mul(::Type{E}, a::ComplexF64, b::E) where {E <: ExactScalar} =
    scalar_mul(E, b, a)

@inline function scalar_add(::Type{E}, a::ComplexF64, b::ComplexF64) where {E}
    s = a + b
    native_sum_exact(s) && return normalize_scalar(s)
    (is_exact_float(a) && is_exact_float(b)) || return normalize_scalar(s)
    return exact_integer(E, a) + exact_integer(E, b)
end
@inline scalar_add(::Type{E}, a::E, b::E) where {E <: ExactScalar} = a + b
@inline function scalar_add(::Type{E}, a::E, b::ComplexF64) where {E <: ExactScalar}
    is_exact_float(b) && return a + exact_integer(E, b)
    return normalize_scalar(to_float_scalar(a) + b)
end
@inline scalar_add(::Type{E}, a::ComplexF64, b::E) where {E <: ExactScalar} =
    scalar_add(E, b, a)

@inline function scalar_inv(::Type{E}, z::ComplexF64) where {E}
    is_exact_float(z) && return inv(exact_integer(E, z))
    return inv(z)
end
@inline scalar_inv(::Type{E}, z::E) where {E <: ExactScalar} = inv(z)

function exact_pow(z::ExactScalar, n::Int)
    n == typemin(Int) && throw(OverflowError("exact power with exponent typemin(Int)"))
    return n >= 0 ? z^n : inv(z)^(-n)
end

@inline is_radical_atom(s::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.isconst(s) && s.val isa Int

radical_power(::Type{E}, p::Int, n::Int) where {E <: ExactScalar} = exact_pow(E(p, 0), n)

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
