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

const Coefficient = Union{Number, SymbolicUtils.BasicSymbolic, Coeff}

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
@inline poly_coeff(p::Poly) = Coeff(ZERO_SLOT, p)
@inline symbolic_coeff(x::RawExpression, known_real::Bool) = Coeff(ZERO_SLOT, RawSymbolicCoeff(x, known_real))
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
    return z == x ? native(z) : symbolic_coeff(SymbolicUtils.unwrap(Num(x)), false)
end
function to_cnum(x::Complex)
    z = ComplexF64(x)
    z == x && return native(z)
    raw_re = SymbolicUtils.unwrap(Num(real(x)))
    raw_im = SymbolicUtils.unwrap(Num(imag(x)))
    return symbolic_coeff(raw_re + im * raw_im, false)
end
function to_cnum(x::Complex{Num})
    re, im = real(x), imag(x)
    raw_re, raw_im = SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im)
    if is_phase(raw_re) && is_const_zero(raw_im)
        return phase_coeff(phase_argument(raw_re))
    end
    return cnum(re, im)
end

function mul_by_im(c::Coeff)
    tail = c.tail
    tail isa Native && return mul_native(native_scalar(c), EXACT_IM)
    tail isa Poly && return tiered(poly_scale, tail, EXACT_IM)
    return raw_coeff(im * tail.expr, false)
end

cnum(re::Num, im::Num) =
    to_cnum(SymbolicUtils.unwrap(re)) + mul_by_im(to_cnum(SymbolicUtils.unwrap(im)))

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
@inline is_big(parts::Vector{Poly}) = any(is_big, parts)
@inline is_big(_) = false

@inline as_tier(::Type{E}, p::Poly) where {E <: ExactScalar} = as_tier(E, p.terms)

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

Base.isone(c::Coeff) = isequal(c, CNUM_ONE)
qadjoint(c::Coeff) = inner_adjoint(c)
function inner_adjoint(c::Coeff)
    t = c.tail
    t isa RawSymbolicCoeff || return conj(c)
    return to_cnum(inner_adjoint(t.expr))
end
Base.adjoint(c::Coeff) = conj(c)
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
        return CNUM_ONE / c
    end
    return raw_leaf(normalize_phase(inv(tail.expr)), tail.known_real)
end

Base.:-(a::Coeff, b::Coeff) = a + (-b)
function Base.:/(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return native_div(native_scalar(a), native_scalar(b))
    ta, tb = a.tail, b.tail
    tb isa Native && ta isa Poly && return tiered(divided_terms, ta, native_scalar(b))
    if tb isa Poly && length(tb.terms) == 1
        inverse = inv(b)
        inverse_tail = inverse.tail
        inverse_tail isa Poly || return a * inverse
        ta isa Native && return tiered(poly_scale, inverse_tail, native_scalar(a))
        ta isa Poly && return tiered(poly_mul, ta, inverse_tail)
        any(is_radical_atom, only(inverse_tail.terms).syms) && return a * inverse
    end
    return raw_leaf(normalize_phase(raw_binary(/, a, b)), cnum_is_real(a) && cnum_is_real(b))
end

const CoeffOperand = Union{Number, SymbolicUtils.BasicSymbolic}
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::Coeff, b::CoeffOperand) = $op(a, to_cnum(b))
        Base.$op(a::CoeffOperand, b::Coeff) = $op(to_cnum(a), b)
    end
end

raw_conj(x::Number) = conj(x)
raw_conj(x::RawExpression) = SymbolicUtils.symtype(x) <: Real ? x : fold_tree(conj_leaf, x)
function conj_leaf(x::RawExpression)
    SymbolicUtils.symtype(x) <: Real && return x
    is_phase(x) && return expim_expanded(-phase_argument(x))
    is_call_of(x, conj) && return only(SymbolicUtils.arguments(x))
    if is_call_of(x, complex)
        args = SymbolicUtils.arguments(x)
        return raw_complex(raw_conj(args[1]), -raw_conj(args[2]))
    end
    return conj(x)
end

conj_poly(p::Poly) = tiered(conj_terms, p)

function conj_terms(::Type{E}, source::Vector{Monomial{E}}) where {E}
    terms = Vector{Monomial{E}}(undef, length(source))
    @inbounds for k in eachindex(source)
        m = source[k]
        nsyms = Vector{RawExpression}(undef, length(m.syms))
        nexps = copy(m.exps)
        for i in eachindex(m.syms)
            s = m.syms[i]
            if is_phase(s)
                nsyms[i] = s
                nexps[i] = -nexps[i]
            else
                nsyms[i] = conj_leaf(s)
            end
        end
        sort_factors!(nsyms, nexps)
        terms[k] = Monomial{E}(conj(monomial_scalar(m)), nsyms, nexps)
    end
    return canonical_terms!(terms)
end

@inline function Base.conj(c::Coeff)
    t = c.tail
    t isa Native && return native(conj(native_scalar(c)))
    t isa Poly && return conj_poly(t)
    return raw_leaf(normalize_phase(raw_conj(t.expr)), t.known_real)
end

@inline function iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, NUM_ZERO)
end

@inline function Base.iszero(c::Coeff)
    is_native(c) && return iszero(native_scalar(c))
    c.tail isa Poly && return false
    return is_const_zero(c.tail.expr)
end

@noinline wide_native(op::F, a::ExactComplex, b::ExactComplex) where {F} =
    scalar_coeff(op(BigExactComplex(a), BigExactComplex(b)))

@inline function native_op(op::F, checked::G, a::NativeScalar, b::NativeScalar) where {F, G}
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = checked(a, b)
        return overflow ? wide_native(op, a, b) : native(z)
    end
    return native(op(to_float_scalar(a), to_float_scalar(b)))
end
@inline mul_native(a::NativeScalar, b::NativeScalar) = native_op(*, mul_checked, a, b)
@inline add_native(a::NativeScalar, b::NativeScalar) = native_op(+, add_checked, a, b)

@inline function Base.:*(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return mul_native(native_scalar(a), native_scalar(b))
    return mul_coeff_slow(a, b)
end

@noinline function mul_coeff_slow(a::Coeff, b::Coeff)
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
    return raw_coeff(raw, cnum_is_real(a) && cnum_is_real(b))
end

@inline is_polynomial_coeff(c::Coeff) = !(c.tail isa RawSymbolicCoeff)

@inline negated_terms(terms::Vector{Monomial{E}}) where {E} =
    Monomial{E}[with_scalar(m, normalize_scalar(-monomial_scalar(m))) for m in terms]

@inline function Base.:-(a::Coeff)
    t = a.tail
    t isa Native && return native(-native_scalar(a))
    t isa Poly && return poly_coeff(Poly(negated_terms(t.terms)))
    return raw_coeff(-t.expr, t.known_real)
end

@inline function Base.:+(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return add_native(native_scalar(a), native_scalar(b))
    is_native(a) && iszero(native_scalar(a)) && return b
    is_native(b) && iszero(native_scalar(b)) && return a
    return add_coeff_slow(a, b)
end

@noinline function add_coeff_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return tiered(poly_add, ta, tb)
    elseif ta isa Poly && tb isa Native
        return tiered(constant_added_terms, ta, native_scalar(b))
    elseif tb isa Poly && ta isa Native
        return tiered(constant_added_terms, tb, native_scalar(a))
    end
    sum_is_zero(a, b) && return CNUM_ZERO
    raw = raw_binary(+, a, b)
    return raw_coeff(raw, cnum_is_real(a) && cnum_is_real(b))
end

@inline constant_added_terms(::Type{E}, terms::Vector{Monomial{E}}, z::CoeffScalar) where {E} =
    poly_add(terms, Monomial{E}[constant_term(E, z)])

@inline function pow_cnum_nonnegative(base::Coeff, n::Int)
    result = CNUM_ONE
    while n > 0
        isodd(n) && (result = result * base)
        n >>= 1
        n > 0 && (base = base * base)
    end
    return result
end

function pow_cnum_integer(base::Coeff, n::Int)
    n >= 0 && return pow_cnum_nonnegative(base, n)
    n == typemin(Int) &&
        return CNUM_ONE / (pow_cnum_nonnegative(base, typemax(Int)) * base)
    return CNUM_ONE / pow_cnum_nonnegative(base, -n)
end

Base.:^(base::Coeff, n::Integer) = pow_cnum_integer(base, Int(n))

poly_is_real(p::Poly) = terms_are_real(p.terms)
function terms_are_real(terms::Vector{<:Monomial})
    for monomial in terms
        scalar_is_real(monomial) || return false
        for symbol in monomial.syms
            (is_phase(symbol) || is_imaginary_unit(symbol)) && return false
            SymbolicUtils.symtype(symbol) <: Real || return false
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

LinearAlgebra.dot(a::Coeff, b::Coeff) = conj(a) * b
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
        iszero(A[top, j]) && continue
        term = A[top, j] * LinearAlgebra.det(A[rest, filter(!=(j), cols)])
        acc = isodd(k) ? acc + term : acc - term
    end
    return acc
end

@inline function raw_depends_on(x, variable)
    isequal(x, variable) && return true
    x isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(x) || return false
    return any(argument -> raw_depends_on(argument, variable), SymbolicUtils.arguments(x))
end

function coefficient_depends_on(c::Coeff, variable)
    tail = c.tail
    tail isa Native && return false
    tail isa Poly && return any(m -> any(s -> raw_depends_on(s, variable), m.syms), tail.terms)
    return raw_depends_on(tail.expr, variable)
end
