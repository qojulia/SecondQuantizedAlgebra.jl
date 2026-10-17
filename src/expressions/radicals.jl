const RadicalPrime = Union{Int, BigInt}
@inline radical_prime(p::Int) = p
@inline radical_prime(p::BigInt) = fits_small(p) ? Int(p) : p
@inline radical_atom(p::RadicalPrime) = SymbolicUtils.Const{SymbolicUtils.SymReal}(radical_prime(p))

@inline function is_radical_atom(s::SymbolicUtils.BasicSymbolic)
    s isa SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal} || return false
    SymbolicUtils.isconst(s) || return false
    return s.val isa RadicalPrime
end

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

"""
    fold_radicals(E, scalar, syms, exps)

Move the integer part of every radical exponent into `scalar`, keeping radicals canonical.

A radical atom is the hashconsed `SymbolicUtils.Const{SymReal}(p)` of a prime `p` with an
exponent in `(0, 1)`, so `√6`, `√2·√3` and `√24/2` are the same monomial. A rational
radicand becomes its prime atoms through `numeric_radical`: trial division below `2^16`,
then `Primes.isprime` and Pollard-Brent rho under a fixed step budget for a larger cofactor.
A radicand that cannot be split within the budget stays an unevaluated symbolic leaf.
"""
function fold_radicals(
        ::Type{E}, scalar::Union{ComplexF64, E}, syms::Vector{RawExpression},
        exps::Vector{Rational{Int}},
    ) where {E <: ExactScalar}
    osyms = RawExpression[]
    oexps = Rational{Int}[]
    sizehint!(osyms, length(syms)); sizehint!(oexps, length(exps))
    @inbounds for i in eachindex(syms)
        s = syms[i]
        e = exps[i]
        if is_radical_atom(s)
            q = fld(numerator(e), denominator(e))
            if q != 0
                scalar = scalar_mul(scalar, radical_power(E, s.val::RadicalPrime, q))
                e -= q
            end
            iszero(e) && continue
        end
        push!(osyms, s); push!(oexps, e)
    end
    return (scalar, osyms, oexps)
end

function integer_root(n::Integer, q::Int)
    target = BigInt(n)
    guess = round(BigInt, BigFloat(target)^(one(BigFloat) / q))
    for m in (guess - 1, guess, guess + 1)
        m >= 0 && m^q == target && return m
    end
    return guess
end

const RADICAL_TRIAL_BOUND = 1 << 16
const RADICAL_RHO_BUDGET = 1 << 18

"""
    prime_factorization!(factors, n) -> Bool

Push the `(prime, multiplicity)` pairs of `n > 0` onto `factors`, in increasing order. Trial
division runs below `RADICAL_TRIAL_BOUND`; a leftover cofactor below its square is prime. A
larger cofactor is tested with `Primes.isprime` and split by Pollard-Brent rho.
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
    if n < RADICAL_TRIAL_BOUND^2
        push!(factors, (Int(n), 1))
        return true
    end
    primes = BigInt[]
    split_cofactor!(primes, BigInt(n), Ref(RADICAL_RHO_BUDGET)) || return false
    for prime in unique!(sort(primes))
        push!(factors, (radical_prime(prime), count(==(prime), primes)))
    end
    return true
end

function split_cofactor!(primes::Vector{BigInt}, n::BigInt, budget::Base.RefValue{Int})
    if n < RADICAL_TRIAL_BOUND^2 || (fits_small(n) ? isprime(Int(n)) : isprime(n))
        push!(primes, n)
        return true
    end
    for k in 2:(ndigits(n; base = 2) ÷ 16)
        r = integer_root(n, k)
        if r^k == n
            for _ in 1:k
                split_cofactor!(primes, r, budget) || return false
            end
            return true
        end
    end
    c = 1
    while (d = rho_divisor(n, c, budget)) == n
        c += 1
    end
    iszero(d) && return false
    return split_cofactor!(primes, d, budget) && split_cofactor!(primes, div(n, d), budget)
end

function rho_divisor(n::BigInt, c::Int, budget::Base.RefValue{Int})
    fits_small(n) && return BigInt(pollard_brent(UInt64(n), UInt64(c), budget))
    return pollard_brent(n, BigInt(c), budget)
end

rho_step!(y::UInt64, c::UInt64, n::UInt64, ::UInt64) = ((widemul(y, y) + c) % n) % UInt64
function rho_accumulate!(q::UInt64, x::UInt64, y::UInt64, n::UInt64, ::UInt64)
    return (widemul(q, x > y ? x - y : y - x) % n) % UInt64
end
rho_gcd!(::UInt64, a::UInt64, b::UInt64, n::UInt64, ::UInt64) = gcd(a > b ? a - b : b - a, n)
rho_gcd!(::UInt64, a::UInt64, n::UInt64) = gcd(a, n)
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
rho_gcd!(g::BigInt, a::BigInt, n::BigInt) = MPZ.gcd!(g, a, n)
rho_assign!(x::BigInt, y::BigInt) = MPZ.set!(x, y)

"""
    pollard_brent(n, c, budget)

A nontrivial divisor of the composite `n` from Brent's cycle search on `y ↦ y² + c mod n`,
`n` itself when this `c` fails, or `0` once `budget` steps are spent.
"""
function pollard_brent end
for T in (UInt64, BigInt)
    @eval function pollard_brent(n::$T, c::$T, budget::Base.RefValue{Int})
        batch = 128
        y, x, ys, t = $T(2), $T(0), $T(0), $T(0)
        q, g = $T(1), $T(1)
        r = 1
        while isone(g)
            budget[] -= 2r
            budget[] < 0 && return zero($T)
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
                g = rho_gcd!(g, q, n)
                k += batch
            end
            r *= 2
        end
        g == n || return g
        g = rho_assign!(g, one($T))
        while isone(g)
            ys = rho_step!(ys, c, n, t)
            g = rho_gcd!(g, x, ys, n, t)
        end
        return g
    end
end

function exact_rational_root(value::Rational{BigInt}, r::Rational{Int})
    q = denominator(r)
    rn, rd = integer_root(numerator(value), q), integer_root(denominator(value), q)
    exact = rn^q == numerator(value) && rd^q == denominator(value) && !(iszero(value) && numerator(r) < 0)
    return exact, exact ? (rn // rd)^numerator(r) : zero(Rational{BigInt})
end

function numeric_radical(value::Rational{BigInt}, r::Rational{Int}, x)
    exact, root = exact_rational_root(value, r)
    exact && return to_cnum(root)
    value > 0 || return raw_coeff(x, true)
    numerator_factors = Tuple{RadicalPrime, Int}[]
    denominator_factors = Tuple{RadicalPrime, Int}[]
    prime_factorization!(numerator_factors, numerator(value)) || return raw_coeff(x, true)
    prime_factorization!(denominator_factors, denominator(value)) || return raw_coeff(x, true)
    syms = RawExpression[]
    exps = Rational{Int}[]
    for (factors, sign) in ((numerator_factors, 1), (denominator_factors, -1)), (p, k) in factors
        push!(syms, radical_atom(p)); push!(exps, sign * k * r)
    end
    return tiered(monomial_terms, EXACT_ONE, syms, exps)
end


function negative_radical(value::Rational{BigInt}, r::Rational{Int}, x, real_root::Bool)
    magnitude = numeric_radical(-value, r, x)
    if magnitude.tail isa RawSymbolicCoeff
        real_root && return raw_coeff(x, true)
        return native(ComplexF64(Float64(value))^r)
    end
    unit = real_root ? EXACT_NEG1 : EXACT_IM^numerator(r)
    return magnitude * native(unit)
end

function radical_form(x::RawExpression)
    SymbolicUtils.iscall(x) || return (x, 0 // 1)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if length(args) == 1
        op === sqrt && return (only(args), 1 // 2)
        op === cbrt && return (only(args), 1 // 3)
    elseif op === (^) && length(args) == 2 && is_rational_exponent(args[2])
        return (args[1], rational_exponent(args[2]))
    end
    return (x, 0 // 1)
end

function split_radicals(x::RawExpression)
    factors = is_call_of(x, *) ? SymbolicUtils.arguments(x) : (x,)
    radicals = CNUM_ONE
    rest = RawExpression[]
    for f in factors
        f isa RawExpression || return (CNUM_ONE, RawExpression[x])
        base, r = radical_form(f)
        if !iszero(r) && is_exact_number(base) && exact_number(base) >= 0
            c = numeric_radical(exact_number(base), r, f)
            if is_polynomial_coeff(c)
                radicals = radicals * c
                continue
            end
        end
        push!(rest, f)
    end
    return (radicals, rest)
end

function radical_factors(m::Monomial)
    factors = Tuple{Rational{Int}, RadicalPrime}[]
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        is_radical_atom(s) || continue
        e = m.exps[i]
        p = s.val::RadicalPrime
        k = p isa Int && p < RADICAL_TRIAL_BOUND ? mergeable_factor(factors, e) : 0
        if iszero(k)
            push!(factors, (e, p))
        else
            product, overflow = Base.mul_with_overflow(factors[k][2]::Int, p::Int)
            overflow ? push!(factors, (e, p)) : (factors[k] = (e, product))
        end
    end
    return insertion_sort!(factors, radical_factor_less)
end
function mergeable_factor(factors::Vector{Tuple{Rational{Int}, RadicalPrime}}, e::Rational{Int})
    for k in eachindex(factors)
        (factors[k][1] == e && factors[k][2] isa Int) && return k
    end
    return 0
end
@inline radical_factor_less(a::Tuple{Rational{Int}, RadicalPrime}, b::Tuple{Rational{Int}, RadicalPrime}) =
    a[1] == b[1] ? a[2] < b[2] : a[1] < b[1]

function radical_expression(n::RadicalPrime, f::Rational{Int})
    base = radical_atom(n)
    f == 1 // 2 && return sqrt(base)
    f == 1 // 3 && return cbrt(base)
    return SymbolicUtils.term(^, base, f; type = Real)
end
