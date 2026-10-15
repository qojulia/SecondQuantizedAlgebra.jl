struct Native end
const NATIVE = Native()

struct RawSymbolicCoeff
    expr::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    known_real::Bool

    function RawSymbolicCoeff(
            expr::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}, known_real::Bool = false,
        )
        SymbolicUtils.symtype(expr) <: Number ||
            throw(ArgumentError("a symbolic coefficient must have numeric symtype"))
        isempty(SymbolicUtils.shape(expr)) ||
            throw(ArgumentError("a symbolic coefficient must be scalar"))
        return new(expr, known_real)
    end
end

"""
    Coeff

Exact coefficient of an operator term.

Iterating an expression and [`exponential_form`](@ref) return coefficients in this form. It
holds exact numbers, such as `1//2 + im`, polynomials in the parameters, and any other
symbolic expression, and it supports the usual arithmetic with numbers, parameters and
other coefficients. [`to_num`](@ref) converts it to Symbolics. It is public but not
exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: exponential_form, Coeff

julia> @variables θ;

julia> exponential_form(cos(θ)) isa Coeff
true
```

See also [`to_num`](@ref), [`get_prefactor`](@ref).
"""
struct Coeff
    slot::NativeSlot
    tail::Union{Native, Poly, RawSymbolicCoeff}
end
"""
    CNum

Alias of [`Coeff`](@ref).
"""
const CNum = Coeff

const NUM_ZERO = Num(0)
const NUM_ONE = Num(1)
const NUM_NEG_ONE = Num(-1)
const EXACT_ONE = one(ExactComplex)
const EXACT_IM = ExactComplex(0, 1)
const EXACT_NEG1 = ExactComplex(-1, 0)

@inline native(z::ComplexF64) = Coeff(NativeSlot(normalize_scalar(z)), NATIVE)
@inline native(z::ExactComplex) = Coeff(NativeSlot(z), NATIVE)
@inline native_scalar(c::Coeff) = native_scalar(c.slot)
const ZERO_SLOT = NativeSlot(zero(ExactComplex))
@inline raw_coeff(x::RawExpression, known_real::Bool) =
    Coeff(ZERO_SLOT, RawSymbolicCoeff(x, known_real))
@inline poly_coeff(p::Poly) = Coeff(ZERO_SLOT, p)
@inline is_native(c::Coeff) = c.tail isa Native
@inline is_poly(c::Coeff) = c.tail isa Poly

const CNUM_ZERO = native(zero(ExactComplex))
const CNUM_ONE = native(EXACT_ONE)
const CNUM_NEG1 = native(EXACT_NEG1)
const CNUM_IM = native(EXACT_IM)
const CNUM_NEG_IM = native(-EXACT_IM)
const CNUM_HALF = native(ExactComplex(1 // 2, 0 // 1))

@inline scalar_coeff(x::NativeScalar) = native(x)
function scalar_coeff(x::BigExactComplex)
    fits_small(x) && return native(as_tier(ExactComplex, x))
    return poly_coeff(Poly(Monomial{BigExactComplex}[constant_term(BigExactComplex, x)]))
end

const ExactReal = Union{Integer, Rational}

@inline exact_value(::Type{E}, re::T, im::T) where {E <: ExactScalar, T <: ExactReal} = E(re, im)
exact_coeff(x::Union{ExactReal, Complex{<:ExactReal}}) =
    isinf(x) ? native(ComplexF64(x)) : tiered(exact_value, real(x), imag(x))

to_cnum(x::Coeff) = x
to_cnum(x::Num) = to_cnum(SymbolicUtils.unwrap(x))
to_cnum(x::ExactScalar) = scalar_coeff(x)
@inline to_cnum(x::Int) = is_edge(x) ? exact_coeff(x) : native(unsafe_gaussian(x, 0, 1))
to_cnum(x::ExactReal) = exact_coeff(x)
to_cnum(x::Complex{<:ExactReal}) = exact_coeff(x)
function to_cnum(x::Real)
    z = ComplexF64(x)
    return z == x ? native(z) : raw_coeff(SymbolicUtils.unwrap(Num(x)), false)
end
function to_cnum(x::Complex)
    z = ComplexF64(x)
    z == x && return native(z)
    raw_re = SymbolicUtils.unwrap(Num(real(x)))
    raw_im = SymbolicUtils.unwrap(Num(imag(x)))
    return raw_coeff(raw_re + im * raw_im, false)
end
function to_cnum(x::Complex{Num})
    re, im = real(x), imag(x)
    raw_re, raw_im = SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im)
    if is_phase(raw_re) && is_const_zero(raw_im)
        return phase_coeff(only(SymbolicUtils.arguments(raw_re)))
    end
    return cnum(re, im)
end

function mul_by_im(c::Coeff)
    tail = c.tail
    tail isa Native && return mul_native(native_scalar(c), EXACT_IM)
    tail isa Poly && return tiered(poly_scale, tail, EXACT_IM)
    return from_raw(im * tail.expr, false)
end

cnum(re::Num, im::Num) =
    add_cnum(to_cnum(SymbolicUtils.unwrap(re)), mul_by_im(to_cnum(SymbolicUtils.unwrap(im))))

function from_poly(terms::Vector{Monomial{ExactComplex}})
    isempty(terms) && return CNUM_ZERO
    length(terms) == 1 && isempty(terms[1].syms) && return scalar_coeff(monomial_scalar(terms[1]))
    return poly_coeff(Poly(terms))
end
function from_poly(terms::Vector{Monomial{BigExactComplex}})
    all(fits_small, terms) && return from_poly(as_tier(ExactComplex, terms))
    return poly_coeff(Poly(terms))
end

@inline tier_result(terms::Vector{<:Monomial}) = from_poly(terms)
@inline tier_result(z::CoeffScalar) = scalar_coeff(z)

@inline is_big(::Vector{Monomial{BigExactComplex}}) = true
@inline is_big(::Union{BigExactComplex, BigInt, Rational{BigInt}}) = true
@inline is_big(p::Poly) = p.terms isa Vector{Monomial{BigExactComplex}}
@inline is_big(c::Coeff) = c.tail isa Poly && is_big(c.tail)
@inline is_big(parts::Vector{Coeff}) = any(is_big, parts)
@inline is_big(_) = false

@inline function as_tier(::Type{ExactComplex}, p::Poly)
    terms = p.terms
    terms isa Vector{Monomial{ExactComplex}} && return terms
    throw(OverflowError("exact scalar exceeds the small tier"))
end
@inline as_tier(::Type{BigExactComplex}, p::Poly) = as_tier(BigExactComplex, p.terms)

@inline function tiered(f::F, args::Vararg{Any, N}) where {F, N}
    any(is_big, args) && return big_tier(f, args)
    try
        return tier_result(f(ExactComplex, map(a -> as_tier(ExactComplex, a), args)...))
    catch err
        err isa OverflowError || rethrow()
        return big_tier(f, args)
    end
end
@noinline big_tier(f::F, args::Tuple) where {F} =
    tier_result(f(BigExactComplex, map(a -> as_tier(BigExactComplex, a), args)...))

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = to_cnum(x)

Coeff(x::Number) = to_cnum(x)

function Base.isequal(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Native
        return tb isa Native && isequal(native_scalar(a), native_scalar(b))
    elseif ta isa Poly
        return tb isa Poly && isequal(ta, tb)
    else
        return tb isa RawSymbolicCoeff && ta.known_real == tb.known_real &&
            isequal(ta.expr, tb.expr)
    end
end
Base.:(==)(a::Coeff, b::Coeff) = isequal(a, b)
function Base.hash(c::Coeff, h::UInt)
    t = c.tail
    t isa Native && return hash(native_scalar(c), hash(:CoeffNative, h))
    t isa Poly && return hash(t, hash(:CoeffSym, h))
    return hash(t.known_real, hash(t.expr, hash(:CoeffRaw, h)))
end

for op in (:isequal, :(==))
    @eval begin
        Base.$op(a::Coeff, b::Number) = isequal(a, to_cnum(b))
        Base.$op(a::Number, b::Coeff) = isequal(to_cnum(a), b)
    end
end

Base.broadcastable(c::Coeff) = Ref(c)
Base.copy(c::Coeff) = c

Base.iszero(c::Coeff) = iszero_cnum(c)
Base.isone(c::Coeff) = isequal(c, CNUM_ONE)
Base.conj(c::Coeff) = conj_cnum(c)
qadjoint(c::Coeff) = inner_adjoint(c)
function inner_adjoint(c::Coeff)
    t = c.tail
    t isa RawSymbolicCoeff || return conj_cnum(c)
    return to_cnum(inner_adjoint(t.expr))
end
Base.adjoint(c::Coeff) = conj_cnum(c)
Base.transpose(c::Coeff) = c
Base.zero(::Type{Coeff}) = CNUM_ZERO
Base.one(::Type{Coeff}) = CNUM_ONE
Base.zero(::Coeff) = CNUM_ZERO
Base.one(::Coeff) = CNUM_ONE
Base.oneunit(::Type{Coeff}) = CNUM_ONE
Base.oneunit(::Coeff) = CNUM_ONE

function native_div(a::NativeScalar, b::NativeScalar)
    if a isa ExactComplex && b isa ExactComplex && !iszero(b)
        return tiered(exact_quotient, a, b)
    end
    return native(to_float_scalar(a) / to_float_scalar(b))
end
@inline exact_quotient(::Type{E}, a::E, b::E) where {E} = a * inv(b)

function inverse_terms(::Type{E}, terms::Vector{Monomial{E}}) where {E}
    monomial = only(terms)
    return monomial_terms(E, inv(monomial_scalar(monomial)), monomial.syms, -monomial.exps)
end
@inline divided_terms(::Type{E}, terms::Vector{Monomial{E}}, z::CoeffScalar) where {E} =
    poly_scale(terms, inv(as_tier(E, z)))

function Base.inv(c::Coeff)
    tail = c.tail
    tail isa Native && return native_div(EXACT_ONE, native_scalar(c))
    if tail isa Poly
        length(tail.terms) == 1 && return tiered(inverse_terms, tail)
        return div_cnum(CNUM_ONE, c)
    end
    return normalized_leaf(inv(tail.expr), tail.known_real)
end

Base.:-(c::Coeff) = neg_cnum(c)
Base.:+(a::Coeff, b::Coeff) = add_cnum(a, b)
Base.:-(a::Coeff, b::Coeff) = add_cnum(a, neg_cnum(b))
Base.:*(a::Coeff, b::Coeff) = mul_cnum(a, b)
function div_cnum(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return native_div(native_scalar(a), native_scalar(b))
    ta, tb = a.tail, b.tail
    tb isa Native && ta isa Poly && return tiered(divided_terms, ta, native_scalar(b))
    if tb isa Poly && length(tb.terms) == 1
        inverse = inv(b)
        inverse_tail = inverse.tail
        inverse_tail isa Poly || return mul_cnum(a, inverse)
        ta isa Native && return tiered(poly_scale, inverse_tail, native_scalar(a))
        ta isa Poly && return tiered(poly_mul, ta, inverse_tail)
        any(is_radical_atom, only(inverse_tail.terms).syms) && return mul_cnum(a, inverse)
    end
    return normalized_leaf(raw_binary(/, a, b), cnum_is_real(a) && cnum_is_real(b))
end
Base.:/(a::Coeff, b::Coeff) = div_cnum(a, b)

const CoeffOperand = Union{Number, SymbolicUtils.BasicSymbolic}
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::Coeff, b::CoeffOperand) = $op(a, to_cnum(b))
        Base.$op(a::CoeffOperand, b::Coeff) = $op(to_cnum(a), b)
    end
end

function raw_conj(x)
    x isa Number && return conj(x)
    SymbolicUtils.symtype(x) <: Real && return x
    is_phase(x) && return expim_expanded(-only(SymbolicUtils.arguments(x)))
    is_call_of(x, conj) && return only(SymbolicUtils.arguments(x))
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === complex
            return raw_complex(raw_conj(args[1]), -raw_conj(args[2]))
        end
        op === (+) && return foldl(+, map(raw_conj, args))
        op === (*) && return foldl(*, map(raw_conj, args))
        if op === (/)
            return raw_conj(args[1]) / raw_conj(args[2])
        end
    end
    return conj(x)
end

@inline function conj_atom(s::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.symtype(s) <: Real && return s
    is_call_of(s, conj) && return SymbolicUtils.arguments(s)[1]
    return SymbolicUtils.unwrap(conj(s))
end

conj_poly(p::Poly) = tiered(conj_terms, p)

function conj_terms(::Type{E}, source::Vector{Monomial{E}}) where {E}
    terms = Vector{Monomial{E}}(undef, length(source))
    @inbounds for k in eachindex(source)
        m = source[k]
        n = length(m.syms)
        if n == 0
            terms[k] = conj_monomial(m, m.syms, m.exps)
            continue
        end
        nsyms = Vector{RawExpression}(undef, n)
        nexps = copy(m.exps)
        for i in 1:n
            s = m.syms[i]
            if is_phase(s)
                nsyms[i] = s
                nexps[i] = -nexps[i]
            else
                nsyms[i] = conj_atom(s)
            end
        end
        sort_factors!(nsyms, nexps)
        terms[k] = conj_monomial(m, nsyms, nexps)
    end
    return canonical_terms!(terms)
end

@inline function conj_cnum(c::Coeff)
    t = c.tail
    t isa Native && return native(conj(native_scalar(c)))
    t isa Poly && return conj_poly(t)
    return normalized_leaf(raw_conj(t.expr), t.known_real)
end

@inline function iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, NUM_ZERO)
end

@inline function iszero_cnum(c::Coeff)
    is_native(c) && return iszero(native_scalar(c))
    c.tail isa Poly && return false
    return is_const_zero(c.tail.expr)
end

@inline function raw_is_negative_of(
        positive::SymbolicUtils.BasicSymbolic, negative::SymbolicUtils.BasicSymbolic,
    )
    SymbolicUtils.isconst(positive) && SymbolicUtils.isconst(negative) &&
        return is_const_zero(positive + negative)
    SymbolicUtils.iscall(negative) || return false
    op = SymbolicUtils.operation(negative)
    args = SymbolicUtils.arguments(negative)
    if op === (*) && length(args) >= 2 && is_const_minus_one(args[1])
        if SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (*)
            positive_args = SymbolicUtils.arguments(positive)
            length(args) - 1 == length(positive_args) || return false
            @inbounds for i in eachindex(positive_args)
                isequal(args[i + 1], positive_args[i]) || return false
            end
            return true
        end
        return length(args) == 2 && isequal(args[2], positive)
    end
    if op === (+) && SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (+)
        positive_args = SymbolicUtils.arguments(positive)
        length(args) == length(positive_args) || return false
        @inbounds for i in eachindex(positive_args)
            raw_is_negative_of(positive_args[i], args[i]) || return false
        end
        return true
    end
    if op === (/) && length(args) == 2 &&
            SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (/)
        positive_args = SymbolicUtils.arguments(positive)
        length(positive_args) == 2 || return false
        return isequal(args[2], positive_args[2]) &&
            (
            raw_is_negative_of(positive_args[1], args[1]) ||
                raw_is_negative_of(args[1], positive_args[1])
        )
    end
    return false
end

@inline function isneg_cnum(a::RawSymbolicCoeff, b::RawSymbolicCoeff)
    return a.known_real == b.known_real &&
        (
        raw_is_negative_of(a.expr, b.expr) ||
            raw_is_negative_of(b.expr, a.expr)
    )
end

@noinline wide_product(a::ExactComplex, b::ExactComplex) =
    scalar_coeff(BigExactComplex(a) * BigExactComplex(b))
@noinline wide_sum(a::ExactComplex, b::ExactComplex) =
    scalar_coeff(BigExactComplex(a) + BigExactComplex(b))

@inline function mul_native(a::NativeScalar, b::NativeScalar)
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = mul_checked(a, b)
        return overflow ? wide_product(a, b) : native(z)
    end
    return native(to_float_scalar(a) * to_float_scalar(b))
end
@inline function add_native(a::NativeScalar, b::NativeScalar)
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = add_checked(a, b)
        return overflow ? wide_sum(a, b) : native(z)
    end
    return native(to_float_scalar(a) + to_float_scalar(b))
end

@inline function mul_cnum(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return mul_native(native_scalar(a), native_scalar(b))
    return mul_cnum_slow(a, b)
end

@noinline function mul_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return tiered(poly_mul, ta, tb)
    elseif ta isa Poly && tb isa Native
        return tiered(poly_scale, ta, native_scalar(b))
    elseif tb isa Poly && ta isa Native
        return tiered(poly_scale, tb, native_scalar(a))
    end
    raw = if is_native(a) || is_native(b)
        raw_binary(*, a, b)
    else
        raw_product(raw_tail(a), raw_tail(b))
    end
    return from_raw(raw, cnum_is_real(a) && cnum_is_real(b))
end

@inline is_polynomial_coeff(c::Coeff) = !(c.tail isa RawSymbolicCoeff)

@inline negated_terms(terms::Vector{Monomial{E}}) where {E} =
    Monomial{E}[with_scalar(m, normalize_scalar(-monomial_scalar(m))) for m in terms]

@inline function neg_cnum(a::Coeff)
    t = a.tail
    t isa Native && return native(-native_scalar(a))
    t isa Poly && return poly_coeff(Poly(negated_terms(t.terms)))
    return from_raw(-t.expr, t.known_real)
end

@inline function add_cnum(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return add_native(native_scalar(a), native_scalar(b))
    is_native(a) && iszero(native_scalar(a)) && return b
    is_native(b) && iszero(native_scalar(b)) && return a
    return add_cnum_slow(a, b)
end

@noinline function add_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return tiered(poly_add, ta, tb)
    elseif ta isa Poly && tb isa Native
        return tiered(constant_added_terms, ta, native_scalar(b))
    elseif tb isa Poly && ta isa Native
        return tiered(constant_added_terms, tb, native_scalar(a))
    end
    (ta isa RawSymbolicCoeff && tb isa RawSymbolicCoeff && isneg_cnum(ta, tb)) &&
        return CNUM_ZERO
    rational_sum_is_zero(a, b) && return CNUM_ZERO
    raw = raw_binary(+, a, b)
    return from_raw(raw, cnum_is_real(a) && cnum_is_real(b))
end

@inline constant_added_terms(::Type{E}, terms::Vector{Monomial{E}}, z::CoeffScalar) where {E} =
    poly_add(terms, Monomial{E}[constant_term(E, z)])

@inline function pow_cnum_nonnegative(base::Coeff, n::Int)
    result = CNUM_ONE
    while n > 0
        isodd(n) && (result = mul_cnum(result, base))
        n >>= 1
        n > 0 && (base = mul_cnum(base, base))
    end
    return result
end

function pow_cnum_integer(base::Coeff, n::Int)
    n >= 0 && return pow_cnum_nonnegative(base, n)
    n == typemin(Int) &&
        return div_cnum(CNUM_ONE, mul_cnum(pow_cnum_nonnegative(base, typemax(Int)), base))
    return div_cnum(CNUM_ONE, pow_cnum_nonnegative(base, -n))
end

Base.:^(base::Coeff, n::Integer) = pow_cnum_integer(base, Int(n))

poly_is_real(p::Poly) = terms_are_real(p.terms)
function terms_are_real(terms::Vector{<:Monomial})
    for monomial in terms
        scalar_is_real(monomial) || return false
        for symbol in monomial.syms
            (is_phase(symbol) || is_imaginary_unit(symbol)) && return false
        end
    end
    return true
end

@inline raw_is_real(t::RawSymbolicCoeff) =
    t.known_real || SymbolicUtils.symtype(t.expr) <: Real

@inline function cnum_is_real(c::Coeff)
    t = c.tail
    t isa Native && return iszero(imag(native_scalar(c)))
    t isa Poly && return poly_is_real(t)
    return raw_is_real(t)
end

LinearAlgebra.dot(a::Coeff, b::Coeff) = mul_cnum(conj_cnum(a), b)
LinearAlgebra.dot(a::Coeff, b::Number) = LinearAlgebra.dot(a, to_cnum(b))
LinearAlgebra.dot(a::Number, b::Coeff) = LinearAlgebra.dot(to_cnum(a), b)

function LinearAlgebra.norm(c::Coeff, p::Real = 2)
    p == 0 && return iszero(c) ? NUM_ZERO : NUM_ONE
    return abs(c)
end

LinearAlgebra.symmetric(c::Coeff, ::Symbol = :U) = c
LinearAlgebra.symmetric_type(::Type{Coeff}) = Coeff
LinearAlgebra.hermitian(c::Coeff, ::Symbol = :U) = to_cnum(real(c))
LinearAlgebra.hermitian_type(::Type{Coeff}) = Coeff

LinearAlgebra.rmul!(A::AbstractArray{Coeff}, b::Coeff) = (A .= A .* b)
LinearAlgebra.lmul!(a::Coeff, B::AbstractArray{Coeff}) = (B .= a .* B)

function LinearAlgebra.det(A::AbstractMatrix{Coeff})
    LinearAlgebra.checksquare(A)
    rows, cols = axes(A)
    isempty(rows) && return CNUM_ONE
    length(rows) == 1 && return A[first(rows), first(cols)]
    if istriu(A) || istril(A)
        return foldl((acc, ij) -> acc * A[ij[1], ij[2]], zip(rows, cols); init = CNUM_ONE)
    end
    top = first(rows)
    rest = rows[(begin + 1):end]
    acc = CNUM_ZERO
    for (k, j) in enumerate(cols)
        term = A[top, j] * LinearAlgebra.det(A[rest, filter(!=(j), cols)])
        acc = isodd(k) ? acc + term : acc - term
    end
    return acc
end
