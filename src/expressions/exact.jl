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
const NativeScalar = Union{ComplexF64, ExactComplex}

struct NativeSlot
    re::Int
    im::Int
    den::Int
end
const FLOAT_DENOMINATOR = 0

@inline holds_float(s::NativeSlot) = s.den == FLOAT_DENOMINATOR
@inline NativeSlot(z::ExactComplex) = NativeSlot(z.re, z.im, z.den)
@inline NativeSlot(z::ComplexF64) =
    NativeSlot(reinterpret(Int, real(z)), reinterpret(Int, imag(z)), FLOAT_DENOMINATOR)
@inline function native_scalar(s::NativeSlot)::NativeScalar
    holds_float(s) && return ComplexF64(reinterpret(Float64, s.re), reinterpret(Float64, s.im))
    return unsafe_gaussian(s.re, s.im, s.den)
end

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactScalar) = z

@inline fits_small(x::BigInt) = typemin(Int) < x <= typemax(Int)
@inline fits_small(x::Rational{BigInt}) = fits_small(numerator(x)) && fits_small(denominator(x))
@inline fits_small(z::BigExactComplex) = fits_small(z.re) && fits_small(z.im) && fits_small(z.den)

@inline as_tier(::Type{E}, z::E) where {E <: ExactScalar} = z
@inline as_tier(::Type{<:ExactScalar}, z::ComplexF64) = z
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
