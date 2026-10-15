@inline num_from_scalar(x::Union{Int, Float64}) = Num(x)
@inline num_from_scalar(x::Rational{Int}) =
    denominator(x) == 1 ? Num(Int(numerator(x))) : Num(x)
function num_from_scalar(x::Rational{BigInt})
    fits_small(x) && return num_from_scalar(Rational{Int}(x))
    return denominator(x) == 1 ? Num(numerator(x)) : Num(x)
end

function radical_factors(m::Monomial)
    factors = Tuple{Rational{Int}, RadicalPrime}[]
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        is_radical_atom(s) || continue
        e = m.exps[i]
        p = s.val::RadicalPrime
        mergeable = p isa Int && p < RADICAL_TRIAL_BOUND
        k = mergeable ? findfirst(f -> f[1] == e && f[2] isa Int, factors) : nothing
        if k === nothing
            push!(factors, (e, p))
        else
            product, overflow = Base.mul_with_overflow(factors[k][2]::Int, p::Int)
            overflow ? push!(factors, (e, p)) : (factors[k] = (e, product))
        end
    end
    return insertion_sort!(factors, radical_factor_less)
end
@inline radical_factor_less(a::Tuple{Rational{Int}, RadicalPrime}, b::Tuple{Rational{Int}, RadicalPrime}) =
    a[1] == b[1] ? a[2] < b[2] : a[1] < b[1]

function radical_expression(n::RadicalPrime, f::Rational{Int})
    base = SymbolicUtils.Const{SymbolicUtils.SymReal}(n)
    f == 1 // 2 && return sqrt(base)
    f == 1 // 3 && return cbrt(base)
    return SymbolicUtils.term(^, base, f; type = Real)
end

function term_to_num(m::Monomial)
    prod = NUM_ONE
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        e = m.exps[i]
        is_radical_atom(s) && continue
        if e < 0 && is_phase(s)
            s = expim_symbolic(-only(SymbolicUtils.arguments(s)))
            e = -e
        end
        base = Num(s)
        q = div(numerator(e), denominator(e))
        if q >= 0
            for _ in 1:q
                prod = prod * base
            end
        else
            for _ in 1:(-q)
                prod = prod / base
            end
        end
        f = e - q
        if f == 1 // 2
            prod = prod * sqrt(base)
        elseif f == -1 // 2
            prod = prod / sqrt(base)
        elseif f != 0
            prod = prod * Num(SymbolicUtils.unwrap(base)^f)
        end
    end
    for (f, n) in radical_factors(m)
        prod = prod * Num(radical_expression(n, f))
    end
    scalar = monomial_scalar(m)
    real_scalar, imag_scalar = real(scalar), imag(scalar)
    re = iszero(real_scalar) ? NUM_ZERO : num_from_scalar(real_scalar) * prod
    imag_ = iszero(imag_scalar) ? NUM_ZERO : num_from_scalar(imag_scalar) * prod
    return Complex(re, imag_)
end

poly_to_num(p::Poly) = terms_to_num(p.terms)
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
@inline with_raw_number(g::G, z::ComplexF64) where {G} = g(z)

function raw_power(symbol::RawExpression, exponent::Rational{Int})
    denominator(exponent) == 1 && return symbol^numerator(exponent)
    return symbol^exponent
end

function factors_to_raw(m::Monomial)
    factors = RawExpression[]
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        is_radical_atom(symbol) && continue
        push!(factors, raw_power(symbol, m.exps[i]))
    end
    for (f, n) in radical_factors(m)
        push!(factors, radical_expression(n, f))
    end
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

function Base.real(c::Coeff)
    is_native(c) && return num_from_scalar(real(native_scalar(c)))
    tail = c.tail
    if tail isa Poly
        is_pure_phase(tail) && return pure_phase_real(tail)
        return real(poly_to_num(tail))
    end
    expr = normalize_phase(raw_tail(c))
    re, _ = cnum_is_real(c) ? (expr, RAW_ZERO) : raw_realimag(expr)
    return Num(normalize_phase(re))
end

function Base.imag(c::Coeff)
    is_native(c) && return num_from_scalar(imag(native_scalar(c)))
    tail = c.tail
    if tail isa Poly
        is_pure_phase(tail) && return pure_phase_imag(tail)
        return imag(poly_to_num(tail))
    end
    expr = normalize_phase(raw_tail(c))
    _, im = cnum_is_real(c) ? (RAW_ZERO, RAW_ZERO) : raw_realimag(expr)
    return Num(normalize_phase(im))
end

for f in (:abs, :abs2)
    @eval function Base.$f(c::Coeff)
        is_native(c) && return num_from_scalar($f(native_scalar(c)))
        is_pure_phase(c) || throw(MethodError($f, (c,)))
        return NUM_ONE
    end
end

@inline function realimag(c::Coeff)
    is_native(c) && return (num_from_scalar(real(native_scalar(c))), num_from_scalar(imag(native_scalar(c))))
    return reim(to_num(c))
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
    expr = normalize_phase(t.expr)
    re, im = cnum_is_real(c) ? (expr, RAW_ZERO) : raw_realimag(expr)
    return Complex(Num(normalize_phase(re)), Num(normalize_phase(im)))
end


Base.show(io::IO, c::Coeff) = show(io, to_num(c))

@inline function rewrite_complex_slots(x)
    x isa SymbolicUtils.BasicSymbolic || return x
    SymbolicUtils.iscall(x) || return x
    SymbolicUtils.operation(x) === complex || return x
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return x
    return args[1] + im * args[2]
end

const COMPLEX_SLOT_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(rewrite_complex_slots),
)

@inline lower_complex_slots(x) = COMPLEX_SLOT_REWRITER(x)
@inline derivative_expression(c::Coeff) = lower_complex_slots(raw_tail(c))

(D::Symbolics.Differential)(c::Coeff) =
    is_native(c) ? D(to_num(c)) : D(derivative_expression(c))

function Symbolics.derivative(c::Coeff, var; simplify = false, kwargs...)
    is_native(c) && return CNUM_ZERO
    input = derivative_expression(c)
    differentiated = Symbolics.expand_derivatives(
        Symbolics.Differential(var)(input), simplify; kwargs...,
    )
    return to_cnum(differentiated)
end
