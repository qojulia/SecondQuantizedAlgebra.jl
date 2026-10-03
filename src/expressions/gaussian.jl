"""
    GaussianRational{T}

An exact Gaussian rational `(re + im*i) / den`, an element of `ℚ(i)`. The denominator is
positive and `gcd(re, im, den) == 1`, so every value has exactly one representation and
`isequal` compares fields. A Gaussian integer has `den == 1`, and its arithmetic needs no
`gcd`.

It is a storage type, not a `Number`: methods of a new `Number` subtype on `==`, `hash` and
the arithmetic operators invalidate SymbolicUtils code that dispatches on abstract `Number`
arguments, which slowed its `simplify` by about a fifth. An exact value enters a symbolic
expression as a `Rational` or `Complex{Rational}` instead.

`T = Int` is the small tier. Its arithmetic is checked and throws `OverflowError`, like
`Rational{Int}`; the kernels `mul_checked`/`add_checked` report the overflow as a flag
instead. No field of a small-tier value is `typemin(Int)`, so negation and conjugation
never overflow. `T = BigInt` is the big tier.
"""
struct Unchecked end

struct GaussianRational{T <: Union{Int, BigInt}}
    re::T
    im::T
    den::T

    # For fields that already satisfy the invariants; `unsafe_gaussian` is the entry point.
    GaussianRational{T}(::Unchecked, re::T, im::T, den::T) where {T <: Union{Int, BigInt}} =
        new{T}(re, im, den)
end

@inline unsafe_gaussian(re::T, im::T, den::T) where {T <: Union{Int, BigInt}} =
    GaussianRational{T}(Unchecked(), re, im, den)

@inline is_edge(x::Int) = x == typemin(Int)
@inline is_edge(::BigInt) = false

# Divide out the common content. `den > 0` and no field is an edge value.
@inline function normalized(re::T, im::T, den::T) where {T}
    isone(den) && return unsafe_gaussian(re, im, den)
    g = gcd(gcd(re, im), den)
    isone(g) && return unsafe_gaussian(re, im, den)
    return unsafe_gaussian(div(re, g), div(im, g), div(den, g))
end

@inline overflow_result(z::GaussianRational, overflow::Bool) =
    overflow ? throw(OverflowError("Gaussian rational exceeds Int")) : z

# An integer as a field of tier `T`; a value outside the small tier throws `OverflowError`.
@inline function field(::Type{Int}, x::Integer)
    typemin(Int) < x <= typemax(Int) || throw(OverflowError("Gaussian rational exceeds Int"))
    return Int(x)
end
@inline field(::Type{BigInt}, x::Integer) = BigInt(x)

GaussianRational{T}(re::Integer, im::Integer) where {T} =
    unsafe_gaussian(field(T, re), field(T, im), one(T))
function GaussianRational{T}(re::Rational, im::Rational) where {T}
    a, b = field(T, denominator(re)), field(T, denominator(im))
    g = gcd(a, b)
    den, o1 = Base.mul_with_overflow(div(a, g), b)
    r, o2 = Base.mul_with_overflow(field(T, numerator(re)), div(b, g))
    i, o3 = Base.mul_with_overflow(field(T, numerator(im)), div(a, g))
    overflow = o1 | o2 | o3 | is_edge(den) | is_edge(r) | is_edge(i)
    return overflow_result(overflow ? unsafe_gaussian(r, i, den) : normalized(r, i, den), overflow)
end
GaussianRational{T}(z::Complex{<:Integer}) where {T} = GaussianRational{T}(real(z), imag(z))
GaussianRational{T}(z::Complex{<:Rational}) where {T} = GaussianRational{T}(real(z), imag(z))
GaussianRational{T}(x::Union{Integer, Rational}) where {T} = GaussianRational{T}(x, zero(x))
GaussianRational{T}(re::Integer, im::Rational) where {T} = GaussianRational{T}(re // 1, im)
GaussianRational{T}(re::Rational, im::Integer) where {T} = GaussianRational{T}(re, im // 1)
GaussianRational{BigInt}(z::GaussianRational{Int}) =
    unsafe_gaussian(BigInt(z.re), BigInt(z.im), BigInt(z.den))

Base.zero(::Type{GaussianRational{T}}) where {T} = unsafe_gaussian(zero(T), zero(T), one(T))
Base.one(::Type{GaussianRational{T}}) where {T} = unsafe_gaussian(one(T), zero(T), one(T))
Base.zero(z::GaussianRational) = zero(typeof(z))
Base.one(z::GaussianRational) = one(typeof(z))
Base.iszero(z::GaussianRational) = iszero(z.re) && iszero(z.im)
Base.isone(z::GaussianRational) = isone(z.re) && iszero(z.im) && isone(z.den)
Base.isreal(z::GaussianRational) = iszero(z.im)
Base.real(z::GaussianRational) = z.re // z.den
Base.imag(z::GaussianRational) = z.im // z.den

# Checked kernels: the result and whether an `Int` operation overflowed. Gaussian integers,
# the common case, skip the denominator arithmetic and the normalization.
@inline function mul_checked(a::GaussianRational{T}, b::GaussianRational{T}) where {T}
    (isone(a.den) && isone(b.den)) || return mul_checked_rational(a, b)
    if iszero(a.im) && iszero(b.im)   # real integers: one product
        re, o = Base.mul_with_overflow(a.re, b.re)
        overflow = o | is_edge(re)
        return (overflow ? a : unsafe_gaussian(re, zero(T), one(T)), overflow)
    end
    rr, o1 = Base.mul_with_overflow(a.re, b.re)
    ii, o2 = Base.mul_with_overflow(a.im, b.im)
    ri, o3 = Base.mul_with_overflow(a.re, b.im)
    ir, o4 = Base.mul_with_overflow(a.im, b.re)
    re, o5 = Base.sub_with_overflow(rr, ii)
    im, o6 = Base.add_with_overflow(ri, ir)
    overflow = o1 | o2 | o3 | o4 | o5 | o6 | is_edge(re) | is_edge(im)
    return (overflow ? a : unsafe_gaussian(re, im, one(T)), overflow)
end
@noinline function mul_checked_rational(
        a::GaussianRational{T}, b::GaussianRational{T},
    ) where {T}
    rr, o1 = Base.mul_with_overflow(a.re, b.re)
    ii, o2 = Base.mul_with_overflow(a.im, b.im)
    ri, o3 = Base.mul_with_overflow(a.re, b.im)
    ir, o4 = Base.mul_with_overflow(a.im, b.re)
    re, o5 = Base.sub_with_overflow(rr, ii)
    im, o6 = Base.add_with_overflow(ri, ir)
    den, o7 = Base.mul_with_overflow(a.den, b.den)
    overflow = o1 | o2 | o3 | o4 | o5 | o6 | o7 | is_edge(re) | is_edge(im) | is_edge(den)
    overflow && return (a, true)
    return (normalized(re, im, den), false)
end

@inline function add_checked(a::GaussianRational{T}, b::GaussianRational{T}) where {T}
    if a.den == b.den
        re, o1 = Base.add_with_overflow(a.re, b.re)
        im, o2 = Base.add_with_overflow(a.im, b.im)
        overflow = o1 | o2 | is_edge(re) | is_edge(im)
        overflow && return (a, true)
        return (normalized(re, im, a.den), false)
    end
    g = gcd(a.den, b.den)
    da, db = div(a.den, g), div(b.den, g)
    x1, o1 = Base.mul_with_overflow(a.re, db)
    x2, o2 = Base.mul_with_overflow(b.re, da)
    y1, o3 = Base.mul_with_overflow(a.im, db)
    y2, o4 = Base.mul_with_overflow(b.im, da)
    re, o5 = Base.add_with_overflow(x1, x2)
    im, o6 = Base.add_with_overflow(y1, y2)
    den, o7 = Base.mul_with_overflow(a.den, db)
    overflow = o1 | o2 | o3 | o4 | o5 | o6 | o7 | is_edge(re) | is_edge(im) | is_edge(den)
    overflow && return (a, true)
    return (normalized(re, im, den), false)
end

@inline Base.:*(a::GaussianRational{T}, b::GaussianRational{T}) where {T} =
    overflow_result(mul_checked(a, b)...)
@inline Base.:+(a::GaussianRational{T}, b::GaussianRational{T}) where {T} =
    overflow_result(add_checked(a, b)...)
@inline Base.:-(z::GaussianRational) = unsafe_gaussian(-z.re, -z.im, z.den)
@inline Base.:-(a::GaussianRational{T}, b::GaussianRational{T}) where {T} = a + (-b)
@inline Base.conj(z::GaussianRational) = unsafe_gaussian(z.re, -z.im, z.den)

# One method per tier, so that analysis of the method on its own sees concrete fields: a
# `where {T}` signature widens each field to `Union{Int, BigInt}` independently.
for T in (Int, BigInt)
    @eval function Base.inv(z::GaussianRational{$T})
        iszero(z) && throw(DivideError())
        r2, o1 = Base.mul_with_overflow(z.re, z.re)
        i2, o2 = Base.mul_with_overflow(z.im, z.im)
        den, o3 = Base.add_with_overflow(r2, i2)
        re, o4 = Base.mul_with_overflow(z.den, z.re)
        im, o5 = Base.mul_with_overflow(z.den, z.im)
        overflow = o1 | o2 | o3 | o4 | o5 | is_edge(den) | is_edge(re) | is_edge(im)
        return overflow_result(overflow ? z : normalized(re, -im, den), overflow)
    end
end

# Square and multiply. The base is squared only while bits remain, so a power that fits never
# overflows on a discarded square.
function Base.:^(z::GaussianRational, n::Integer)
    n == typemin(n) && throw(OverflowError("exact power with exponent typemin"))
    n < 0 && return inv(z)^(-n)
    result = isodd(n) ? z : one(z)
    base = z
    n >>= 1
    while n > 0
        base *= base
        isodd(n) && (result *= base)
        n >>= 1
    end
    return result
end

# Normalized fields make equality a field comparison. A float scalar of the same polynomial
# compares by value, and `hash` agrees with the `Complex{Rational}` of the same value.
Base.:(==)(a::GaussianRational, b::GaussianRational) =
    a.re == b.re && a.im == b.im && a.den == b.den
Base.isequal(a::GaussianRational, b::GaussianRational) = a == b
Base.:(==)(a::GaussianRational, b::ComplexF64) = Complex(real(a), imag(a)) == b
Base.:(==)(a::ComplexF64, b::GaussianRational) = b == a
Base.isequal(a::GaussianRational, b::ComplexF64) = isequal(Complex(real(a), imag(a)), b)
Base.isequal(a::ComplexF64, b::GaussianRational) = isequal(b, a)
# A Gaussian integer hashes like the `Complex{Int}` of its value, which equals the hash of the
# `Complex{Rational}`; only a proper fraction needs the reduced parts.
Base.hash(z::GaussianRational, h::UInt) =
    isone(z.den) ? hash(Complex(z.re, z.im), h) : hash(Complex(real(z), imag(z)), h)

Base.abs2(z::GaussianRational) = (big(z.re)^2 + big(z.im)^2) // big(z.den)^2
Base.abs(z::GaussianRational) = abs(ComplexF64(z))
Base.Complex{Float64}(z::GaussianRational) =
    ComplexF64(Float64(z.re // z.den), Float64(z.im // z.den))

Base.show(io::IO, z::GaussianRational) = show(io, Complex(real(z), imag(z)))
