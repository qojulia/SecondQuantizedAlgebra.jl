"""
Exact scalar tier: Gaussian-rational arithmetic that stays exact through the small
(`ExactComplex`) and big (`BigExactComplex`) tiers, checked small-tier arithmetic that
promotes on overflow, the native/exact classification of `ComplexF64` (I4/I7), and the
radical-folding primitives that keep a numeric radical atom's exponent in `(0, 1)`
(I2/I5). Nothing here depends on `Monomial`, `Poly`, or `Coeff` — those types are
defined later (`monomial.jl`, `cnum.jl`) and depend on this file, not the other way
around; that dependency direction, not file position within `expressions/`, is why this
file is included first. See `docs/src/devdocs.md` "Exact coefficients" for the
invariants (I1-I7) this machinery maintains.
"""
const ExactComplex = Complex{Rational{Int}}
const BigExactComplex = Complex{Rational{BigInt}}
const CoeffScalar = Union{ComplexF64, ExactComplex, BigExactComplex}
const ExactScalar = Union{ExactComplex, BigExactComplex}
const SmallScalar = Union{ComplexF64, ExactComplex}

const MAX_EXACT_FLOAT = maxintfloat(Float64)

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactComplex) = z
@inline normalize_scalar(z::BigExactComplex) = canonical_exact(z)

@inline fits_int(x::Rational{BigInt}) =
    -typemax(Int) <= numerator(x) <= typemax(Int) && denominator(x) <= typemax(Int)

function canonical_exact(z::BigExactComplex)::ExactScalar
    re, im = real(z), imag(z)
    (fits_int(re) && fits_int(im)) || return z
    return ExactComplex(
        Rational{Int}(Int(numerator(re)), Int(denominator(re))),
        Rational{Int}(Int(numerator(im)), Int(denominator(im))),
    )
end

@inline widen_exact(z::ExactComplex) = BigExactComplex(z)
@inline widen_exact(z::BigExactComplex) = z

@inline function coprime_parts(a::Int, b::Int)
    g = gcd(a, b)
    return (div(a, g), div(b, g))
end

@inline function small_rational(n::Int, d::Int, overflow::Bool)
    return (Base.unsafe_rational(n, d), overflow | (n == typemin(Int)))
end

@inline function checked_rational_mul(x::Rational{Int}, y::Rational{Int})
    xn, yd = coprime_parts(numerator(x), denominator(y))
    xd, yn = coprime_parts(denominator(x), numerator(y))
    n, o1 = Base.Checked.mul_with_overflow(xn, yn)
    d, o2 = Base.Checked.mul_with_overflow(xd, yd)
    return small_rational(n, d, o1 | o2)
end

@inline function checked_rational_add(x::Rational{Int}, y::Rational{Int})
    xd, yd = coprime_parts(denominator(x), denominator(y))
    a, o1 = Base.Checked.mul_with_overflow(numerator(x), yd)
    b, o2 = Base.Checked.mul_with_overflow(numerator(y), xd)
    n, o3 = Base.Checked.add_with_overflow(a, b)
    d, o4 = Base.Checked.mul_with_overflow(denominator(x), yd)
    (o1 | o2 | o3 | o4) && return (zero(Rational{Int}), true)
    g = gcd(n, d)
    return small_rational(div(n, g), div(d, g), false)
end

@inline function checked_exact_mul(a::ExactComplex, b::ExactComplex)
    ar, ai, br, bi = real(a), imag(a), real(b), imag(b)
    rr, o1 = checked_rational_mul(ar, br)
    ii, o2 = checked_rational_mul(ai, bi)
    ri, o3 = checked_rational_mul(ar, bi)
    ir, o4 = checked_rational_mul(ai, br)
    (o1 | o2 | o3 | o4) && return (a, true)
    re, o5 = checked_rational_add(rr, -ii)
    im, o6 = checked_rational_add(ri, ir)
    return (ExactComplex(re, im), o5 | o6)
end

@inline function checked_exact_add(a::ExactComplex, b::ExactComplex)
    re, o1 = checked_rational_add(real(a), real(b))
    im, o2 = checked_rational_add(imag(a), imag(b))
    return (ExactComplex(re, im), o1 | o2)
end

@inline function exact_mul(a::ExactScalar, b::ExactScalar)::ExactScalar
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = checked_exact_mul(a, b)
        overflow || return z
    end
    return canonical_exact(widen_exact(a) * widen_exact(b))
end

@inline function exact_add(a::ExactScalar, b::ExactScalar)::ExactScalar
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = checked_exact_add(a, b)
        overflow || return z
    end
    return canonical_exact(widen_exact(a) + widen_exact(b))
end

function exact_inv(a::ExactComplex)::ExactScalar
    z = try
        inv(a)
    catch err
        err isa OverflowError || rethrow()
        return canonical_exact(inv(widen_exact(a)))
    end
    (numerator(real(z)) == typemin(Int) || numerator(imag(z)) == typemin(Int)) &&
        return widen_exact(z)
    return z
end
exact_inv(a::BigExactComplex)::ExactScalar = canonical_exact(inv(a))

function exact_pow(z::ExactScalar, n::Int)::ExactScalar
    n == typemin(Int) && throw(OverflowError("exact power with exponent typemin(Int)"))
    e = abs(n)
    result::ExactScalar = ExactComplex(1 // 1, 0 // 1)
    base::ExactScalar = z
    while e > 0
        isodd(e) && (result = exact_mul(result, base))
        e >>= 1
        e > 0 && (base = exact_mul(base, base))
    end
    return n >= 0 ? result : exact_inv(result)
end

scalar_conj(a::ExactComplex) = conj(a)
scalar_conj(a::BigExactComplex)::ExactScalar = canonical_exact(conj(a))
scalar_conj(a::ComplexF64) = conj(a)

@inline to_float_scalar(z::ComplexF64) = z
@inline to_float_scalar(z::ExactScalar) = ComplexF64(Float64(real(z)), Float64(imag(z)))

@inline function integer_scalar(z::ComplexF64)
    re, im = real(z), imag(z)
    (abs(re) <= MAX_EXACT_FLOAT && abs(im) <= MAX_EXACT_FLOAT) || return nothing
    (isinteger(re) && isinteger(im)) || return nothing
    return ExactComplex(Int(re) // 1, Int(im) // 1)
end

@inline native_product_exact(a::ComplexF64, b::ComplexF64) =
    (abs(real(a)) + abs(imag(a))) * (abs(real(b)) + abs(imag(b))) <= MAX_EXACT_FLOAT
@inline native_sum_exact(s::ComplexF64) =
    abs(real(s)) < MAX_EXACT_FLOAT && abs(imag(s)) < MAX_EXACT_FLOAT

@noinline function native_mul_wide(a::ComplexF64, b::ComplexF64)::CoeffScalar
    ea, eb = integer_scalar(a), integer_scalar(b)
    (ea === nothing || eb === nothing) && return normalize_scalar(a * b)
    return exact_mul(ea, eb)
end

@noinline function native_add_wide(a::ComplexF64, b::ComplexF64)::CoeffScalar
    ea, eb = integer_scalar(a), integer_scalar(b)
    (ea === nothing || eb === nothing) && return normalize_scalar(a + b)
    return exact_add(ea, eb)
end

@inline function scalar_inv(z::ComplexF64)
    exact = integer_scalar(z)
    return exact === nothing ? inv(z) : exact_inv(exact)
end
@inline scalar_inv(z::ExactScalar) = exact_inv(z)

function scalar_mul(@nospecialize(a::CoeffScalar), @nospecialize(b::CoeffScalar))::CoeffScalar
    if a isa ComplexF64
        if b isa ComplexF64
            native_product_exact(a, b) && return normalize_scalar(a * b)
            return native_mul_wide(a, b)
        end
        return exact_float_mul(b, a)
    end
    b isa ComplexF64 && return exact_float_mul(a, b)
    return exact_mul(a, b)
end

@inline function exact_float_mul(a::ExactScalar, b::ComplexF64)::CoeffScalar
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(to_float_scalar(a) * b) : exact_mul(a, ib)
end

function scalar_add(@nospecialize(a::CoeffScalar), @nospecialize(b::CoeffScalar))::CoeffScalar
    if a isa ComplexF64
        if b isa ComplexF64
            s = a + b
            native_sum_exact(s) && return normalize_scalar(s)
            return native_add_wide(a, b)
        end
        return exact_float_add(b, a)
    end
    b isa ComplexF64 && return exact_float_add(a, b)
    return exact_add(a, b)
end

@inline function exact_float_add(a::ExactScalar, b::ComplexF64)::CoeffScalar
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(to_float_scalar(a) + b) : exact_add(a, ib)
end

@inline is_radical_atom(s::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.isconst(s) && s.val isa Int

radical_power(p::Int, n::Int)::ExactScalar = exact_pow(ExactComplex(p // 1, 0 // 1), n)

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
