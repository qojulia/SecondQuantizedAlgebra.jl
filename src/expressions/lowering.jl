@inline num_from_scalar(x::Union{Int, Float64}) = Num(x)
@inline num_from_scalar(x::Rational{Int}) =
    denominator(x) == 1 ? Num(Int(numerator(x))) : Num(x)
function num_from_scalar(x::Rational{BigInt})
    fits_small(x) && return num_from_scalar(Rational{Int}(x))
    return denominator(x) == 1 ? Num(numerator(x)) : Num(x)
end

is_complex_atom(s::RawExpression) =
    !(SymbolicUtils.symtype(s) <: Real) && !is_phase(s) && !is_imaginary_unit(s)
has_complex_atom(p::Poly) = any(m -> any(is_complex_atom, m.syms), p.terms)

function term_to_num(m::Monomial)
    raw = factors_to_raw(m, num_power)
    scalar = monomial_scalar(m)
    if any(is_complex_atom, m.syms)
        re, im = raw_realimag(raw)
        return Complex(num_from_scalar(real(scalar)), num_from_scalar(imag(scalar))) * Complex(Num(re), Num(im))
    end
    prod = Num(raw)
    real_scalar, imag_scalar = real(scalar), imag(scalar)
    re = iszero(real_scalar) ? NUM_ZERO : num_from_scalar(real_scalar) * prod
    imag_ = iszero(imag_scalar) ? NUM_ZERO : num_from_scalar(imag_scalar) * prod
    return Complex(re, imag_)
end

poly_to_num(p::Poly) = terms_to_num(p.terms)

function poly_part(p::Poly, imaginary::Bool)
    result = monomial_part(p.terms[1], imaginary)
    @inbounds for i in 2:length(p.terms)
        result += monomial_part(p.terms[i], imaginary)
    end
    return result
end

function monomial_part(m::Monomial, imaginary::Bool)
    i = phase_factor_index(m.syms)
    if i == 0
        a = term_to_num(m)
        return imaginary ? imag(a) : real(a)
    end
    re, im = if length(m.syms) == 1
        s = monomial_scalar(m)
        (num_from_scalar(real(s)), num_from_scalar(imag(s)))
    else
        reim(term_to_num(replace_exps(m, i, 0 // 1, 0, 0 // 1)))
    end
    angle, flipped = canonical_angle(Num(only(SymbolicUtils.arguments(m.syms[i]))), numerator(m.exps[i]))
    imaginary && return combine(scale_trig(im, cos, angle), scale_trig(re, sin, angle), flipped)
    return combine(scale_trig(re, cos, angle), scale_trig(im, sin, angle), !flipped)
end

function scale_trig(a::Num, trig::F, angle::Num) where {F}
    iszero_num(a) && return NUM_ZERO
    isequal(a, NUM_ONE) && return trig(angle)
    isequal(a, NUM_NEG_ONE) && return -trig(angle)
    return a * trig(angle)
end

function combine(x::Num, y::Num, subtract::Bool)
    iszero_num(y) && return x
    iszero_num(x) && return subtract ? -y : y
    return subtract ? x - y : x + y
end
function terms_to_num(terms::Vector{<:Monomial})
    isempty(terms) && return Complex(NUM_ZERO, NUM_ZERO)
    acc = term_to_num(terms[1])
    @inbounds for i in 2:length(terms)
        acc = acc + term_to_num(terms[i])
    end
    return acc
end

@inline function with_raw_number(g::G, z::ExactComplex) where {G}
    if iszero(z.im)
        isone(z.den) && return g(z.re)
        return g(z.re // z.den)
    end
    isone(z.den) && return g(Complex(z.re, z.im))
    return g(Complex(z.re // z.den, z.im // z.den))
end
@inline function with_raw_number(g::G, z::BigExactComplex) where {G}
    fits_small(z) && return with_raw_number(g, as_tier(ExactComplex, z))
    iszero(z.im) && return g(real(z))
    return g(Complex(real(z), imag(z)))
end
@inline with_raw_number(g::G, z::ComplexF64) where {G} = iszero(imag(z)) ? g(real(z)) : g(z)

function raw_power(symbol::RawExpression, exponent::Rational{Int})
    denominator(exponent) == 1 && return symbol^numerator(exponent)
    return symbol^exponent
end

function num_power(symbol::RawExpression, exponent::Rational{Int})
    if exponent < 0 && is_phase(symbol)
        symbol = SymbolicUtils.unwrap(expim_symbolic(-only(SymbolicUtils.arguments(symbol))))
        exponent = -exponent
    end
    q = div(numerator(exponent), denominator(exponent))
    whole = q >= 0 ? symbol^q : 1 / symbol^(-q)
    f = exponent - q
    f == 1 // 2 && return whole * sqrt(symbol)
    f == -1 // 2 && return whole / sqrt(symbol)
    f == 0 && return whole
    return whole * symbol^f
end

function factors_to_raw(m::Monomial, power::F = raw_power) where {F}
    factors = RawExpression[]
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        is_radical_atom(symbol) && continue
        push!(factors, power(symbol, m.exps[i]))
    end
    for (f, n) in radical_factors(m)
        push!(factors, radical_expression(n, f))
    end
    isempty(factors) && return RAW_ONE
    return reduce(*, factors)
end

function term_to_raw(m::Monomial)
    factors = factors_to_raw(m)
    scalar = monomial_scalar(m)
    isone(scalar) && return factors
    return with_raw_number(v -> v * factors, scalar)
end

raw_sum(terms::Vector{Monomial{E}}, range::AbstractUnitRange{Int}) where {E} =
    reduce(+, (term_to_raw(terms[i]) for i in range))

poly_to_raw(p::Poly) = terms_to_raw(p.terms)
function terms_to_raw(terms::Vector{Monomial{E}}) where {E}
    isempty(first(terms).syms) || return raw_sum(terms, eachindex(terms))
    constant = monomial_scalar(first(terms))
    length(terms) == 1 &&
        return with_raw_number(SymbolicUtils.Const{SymbolicUtils.SymReal}, constant)
    variable = raw_sum(terms, 2:lastindex(terms))
    return with_raw_number(v -> variable + v, constant)
end

@inline function raw_tail(c::Coeff)
    tail = c.tail
    tail isa Poly && return poly_to_raw(tail)
    tail isa RawSymbolicCoeff && return tail.expr
    return with_raw_number(SymbolicUtils.Const{SymbolicUtils.SymReal}, native_scalar(c))
end

@inline function raw_binary(f::F, a::Coeff, b::Coeff) where {F}
    if is_native(a)
        y = raw_tail(b)
        return with_raw_number(v -> f(v, y), native_scalar(a))
    elseif is_native(b)
        x = raw_tail(a)
        return with_raw_number(v -> f(x, v), native_scalar(b))
    end
    return f(raw_tail(a), raw_tail(b))
end

function raw_realimag(x::RawExpression)
    SymbolicUtils.isconst(x) && return (real(x), imag(x))
    is_imaginary_unit(x) && return (RAW_ZERO, RAW_ONE)
    SymbolicUtils.symtype(x) <: Real && return (x, RAW_ZERO)
    is_phase(x) && begin
        argument = only(SymbolicUtils.arguments(x))
        return (cos(argument), sin(argument))
    end
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        op === complex && return (args[1], args[2])
        if op === (+)
            re, im = RAW_ZERO, RAW_ZERO
            for argument in args
                ar, ai = raw_realimag(argument)
                re += ar
                im += ai
            end
            return (re, im)
        elseif op === (*)
            re, im = RAW_ONE, RAW_ZERO
            for argument in args
                ar, ai = raw_realimag(argument)
                re, im = re * ar - im * ai, re * ai + im * ar
            end
            return (re, im)
        elseif op === (/) && length(args) == 2 &&
                SymbolicUtils.symtype(args[2]) <: Real
            re, im = raw_realimag(args[1])
            return (re / args[2], im / args[2])
        elseif op === conj
            re, im = raw_realimag(only(args))
            return (re, -im)
        end
    end
    return (real(x), imag(x))
end

function raw_parts(c::Coeff)
    expr = normalize_phase(c.tail.expr)
    return cnum_is_real(c) ? (expr, RAW_ZERO) : raw_realimag(expr)
end

function coefficient_part(c::Coeff, imaginary::Bool)
    t = c.tail
    if t isa Native
        z = native_scalar(c)
        return num_from_scalar(imaginary ? imag(z) : real(z))
    end
    t isa Poly && return poly_part(t, imaginary)
    return Num(normalize_phase(raw_parts(c)[imaginary ? 2 : 1]))
end
Base.real(c::Coeff) = coefficient_part(c, false)
Base.imag(c::Coeff) = coefficient_part(c, true)
@inline realimag(c::Coeff) = reim(to_num(c))

for f in (:abs, :abs2)
    @eval function Base.$f(c::Coeff)
        is_native(c) && return num_from_scalar($f(native_scalar(c)))
        is_pure_phase(c) || throw(MethodError($f, (c,)))
        return NUM_ONE
    end
end

"""
    to_num(c::Coeff)

Coefficient as a `Complex{Num}` of Symbolics.

Use it to pass a coefficient read from an expression to code built on Symbolics, such as
`build_function`. It is public but not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: exponential_form, to_num

julia> @variables θ;

julia> to_num(exponential_form(cos(θ))) isa Complex
true
```

See also [`Coeff`](@ref).
"""
function to_num(c::Coeff)
    t = c.tail
    t isa Native && return Complex(num_from_scalar(real(native_scalar(c))), num_from_scalar(imag(native_scalar(c))))
    t isa Poly && return poly_to_num(t)
    re, im = raw_parts(c)
    return Complex(Num(normalize_phase(re)), Num(normalize_phase(im)))
end


Base.show(io::IO, c::Coeff) = show(io, to_num(c))

@inline function rewrite_complex_slots(x)
    x isa RawExpression || return x
    is_call_of(x, complex) || return x
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return x
    return args[1] + im * args[2]
end

const COMPLEX_SLOT_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(rewrite_complex_slots),
)

@inline lower_complex_slots(x) = COMPLEX_SLOT_REWRITER(x)

(D::Symbolics.Differential)(c::Coeff) =
    is_native(c) ? D(to_num(c)) : D(raw_rule(c))

function Symbolics.derivative(c::Coeff, var; simplify = false, kwargs...)
    is_native(c) && return CNUM_ZERO
    input = raw_rule(c)
    differentiated = Symbolics.expand_derivatives(
        Symbolics.Differential(var)(input), simplify; kwargs...,
    )
    return to_cnum(differentiated)
end
