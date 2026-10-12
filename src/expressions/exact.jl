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
const RawExpression = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}

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
@inline function native_scalar(s::NativeSlot)
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

# A radical prime is stored as an `Int` when it fits, else as a `BigInt`. `Const` hashconsing
# keeps `Const(7)` and `Const(big(7))` apart, so the width must be canonical.
const RadicalPrime = Union{Int, BigInt}
@inline radical_prime(p::Integer) = fits_small(BigInt(p)) ? Int(p) : BigInt(p)
@inline radical_prime(p::Int) = p
@inline radical_atom(p::Integer) = SymbolicUtils.Const{SymbolicUtils.SymReal}(radical_prime(p))

@inline function is_radical_atom(s::SymbolicUtils.BasicSymbolic)
    s isa SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal} || return false
    SymbolicUtils.isconst(s) || return false
    v = s.val
    return v isa Int || v isa BigInt
end

# A `BigInt` prime does not fit the small tier, so its power throws `OverflowError`.
radical_power(::Type{E}, p::RadicalPrime, n::Int) where {E <: ExactScalar} = E(p, 0)^n

@inline function needs_radical_fold(
        syms::Vector{RawExpression}, exps::Vector{Rational{Int}},
    )
    @inbounds for i in eachindex(syms)
        is_radical_atom(syms[i]) || continue
        e = exps[i]
        (e <= 0 || e >= 1) && return true
    end
    return false
end

const RADICAL_TRIAL_BOUND = 1 << 16
# Pollard-Brent steps spent splitting one composite cofactor: about 5 ms for a cofactor below
# `2^63` (native arithmetic) and 30 ms above it. That splits off prime factors up to about
# `2^36`; a cofactor it cannot split stays a symbolic leaf, and is retried at each recognition.
const RADICAL_RHO_BUDGET = 1 << 18

"""
    prime_factorization!(factors, n) -> Bool

Push the `(prime, multiplicity)` pairs of `n > 0` onto `factors`, in increasing order. Trial
division runs below `RADICAL_TRIAL_BOUND`; a leftover cofactor below its square is prime. A
larger cofactor is tested with GMP's probable-prime test and split by Pollard-Brent rho.
Returns `false` when a composite cofactor cannot be split within `RADICAL_RHO_BUDGET`.
"""
function prime_factorization!(factors::Vector{Tuple{RadicalPrime, Int}}, n::Integer)
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
    if n < big(RADICAL_TRIAL_BOUND)^2
        push!(factors, (Int(n), 1))
        return true
    end
    primes = BigInt[]
    split_cofactor!(primes, BigInt(n), Ref(RADICAL_RHO_BUDGET)) || return false
    sort!(primes)
    i = 1
    while i <= length(primes)
        j = i
        while j < length(primes) && primes[j + 1] == primes[i]
            j += 1
        end
        push!(factors, (radical_prime(primes[i]), j - i + 1))
        i = j + 1
    end
    return true
end

# Miller-Rabin rounds after GMP's Baillie-PSW test; no BPSW pseudoprime is known.
is_probable_prime(n::BigInt) =
    ccall((:__gmpz_probab_prime_p, Base.GMP.libgmp), Cint, (Ref{BigInt}, Cint), n, 25) > 0

# `n` has no prime factor below `RADICAL_TRIAL_BOUND`.
function split_cofactor!(primes::Vector{BigInt}, n::BigInt, budget::Base.RefValue{Int})
    if n < big(RADICAL_TRIAL_BOUND)^2 || is_probable_prime(n)
        push!(primes, n)
        return true
    end
    # Rho finds a repeated prime only slowly, so perfect powers are peeled off first.
    for k in 2:(ndigits(n; base = 2) ÷ 16)
        r = integer_root(n, k)
        if r^k == n
            for _ in 1:k
                split_cofactor!(primes, r, budget) || return false
            end
            return true
        end
    end
    for c in 1:4
        d = n < typemax(Int) ? BigInt(pollard_brent(UInt64(n), UInt64(c), budget)) :
            pollard_brent(n, BigInt(c), budget)
        iszero(d) && return false
        d == n && continue
        return split_cofactor!(primes, d, budget) && split_cofactor!(primes, div(n, d), budget)
    end
    return false
end

# The cycle search is written against these updates. `UInt64` updates return new values;
# `BigInt` updates write into their first argument through GMP, since an allocating `BigInt`
# step costs about ten times more.
rho_step!(y::UInt64, c::UInt64, n::UInt64, ::UInt64) = ((widemul(y, y) + c) % n) % UInt64
function rho_accumulate!(q::UInt64, x::UInt64, y::UInt64, n::UInt64, ::UInt64)
    return (widemul(q, x > y ? x - y : y - x) % n) % UInt64
end
rho_gcd!(::UInt64, a::UInt64, b::UInt64, n::UInt64, ::UInt64) = gcd(a > b ? a - b : b - a, n)
rho_assign!(::UInt64, y::UInt64) = y

const MPZ = Base.GMP.MPZ
function rho_step!(y::BigInt, c::BigInt, n::BigInt, t::BigInt)
    MPZ.mul!(t, y, y); MPZ.add!(t, t, c)
    return MPZ.fdiv_r!(y, t, n)
end
function rho_accumulate!(q::BigInt, x::BigInt, y::BigInt, n::BigInt, t::BigInt)
    MPZ.sub!(t, x, y); MPZ.mul!(q, q, t)
    return MPZ.fdiv_r!(q, q, n)
end
function rho_gcd!(g::BigInt, a::BigInt, b::BigInt, n::BigInt, t::BigInt)
    MPZ.sub!(t, a, b)
    return MPZ.gcd!(g, t, n)
end
rho_assign!(x::BigInt, y::BigInt) = MPZ.set!(x, y)

"""
    pollard_brent(n, c, budget)

A nontrivial divisor of the composite `n` from Brent's cycle search on `y ↦ y² + c mod n`,
`n` itself when this `c` fails, or `0` once `budget` steps are spent.
"""
function pollard_brent end
# One concrete method per arithmetic, since a method generic over both would be analysed by
# JET with their union.
for T in (UInt64, BigInt)
    @eval function pollard_brent(n::$T, c::$T, budget::Base.RefValue{Int})
        batch = 128
        y, x, ys, t = $T(2), $T(0), $T(0), $T(0)
        q, g = $T(1), $T(1)
        r = 1
        while isone(g)
            x = rho_assign!(x, y)
            for _ in 1:r
                y = rho_step!(y, c, n, t)
            end
            k = 0
            while k < r && isone(g)
                ys = rho_assign!(ys, y)
                for _ in 1:min(batch, r - k)
                    y = rho_step!(y, c, n, t)
                    q = rho_accumulate!(q, x, y, n, t)
                end
                g = rho_gcd!(g, q, zero($T), n, t)
                k += batch
            end
            budget[] -= 2r
            budget[] < 0 && return zero($T)
            r *= 2
        end
        g == n || return g
        # The batch overshot: replay it one step at a time.
        while true
            ys = rho_step!(ys, c, n, t)
            g = rho_gcd!(g, x, ys, n, t)
            isone(g) || return g
        end
        return
    end
end
