"""
    expim(x)

Exact unit phase with real argument `x`.

`expim(x)` stands for ``e^{ix}`` and is kept as an exact symbol rather than a
generic exponential. Multiplying phases adds their arguments and opposite phases
cancel, so the factors produced by drives and frame changes do not pile up.

The argument has to be real; for a complex exponent use `exp`. Fractional powers
are refused, since they would have to pick a branch.

# Examples

```jldoctest
julia> @variables ω t θ;

julia> expim(ω * t) * expim(θ)
exp(im*(θ + t*ω))

julia> conj(expim(ω * t)) * expim(ω * t)
1

julia> a = Destroy(FockSpace(:cavity), :a);

julia> expim(ω * t) * expim(-ω * t) * a
a
```

See also [`trigonometric_form`](@ref), [`exponential_form`](@ref).
"""
expim(x::Real) = exp(im * x)
expim(x::Num) = phase_coeff(x)
expim(x::SymbolicUtils.BasicSymbolic) = phase_coeff(x)

@noinline function nonreal_phase_argument(x)
    throw(ArgumentError("`expim` requires a provably real argument; got `$x`"))
end

expim(x::Number) = nonreal_phase_argument(x)

@inline is_call_of(x::SymbolicUtils.BasicSymbolic, op::F) where {F} =
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === op

@inline is_phase(b) = b isa RawExpression && is_call_of(b, expim)

@inline is_imaginary_unit(x) = isequal(x, Symbolics.IM)

@inline function phase_factor_index(syms)
    @inbounds for i in eachindex(syms)
        is_phase(syms[i]) && return i
    end
    return 0
end

const SCALAR_SHAPE = UnitRange{Int}[]
raw_complex(re, im) = SymbolicUtils.term(
    complex, SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im);
    type = Complex{Real}, shape = SCALAR_SHAPE,
)
expim_expanded(x) = SymbolicUtils.term(
    expim, SymbolicUtils.unwrap(x); type = Complex{Real}, shape = SCALAR_SHAPE,
)
expim_symbolic(x) = expim_expanded(expand(x))

SymbolicUtils.promote_symtype(::typeof(expim), ::SymbolicUtils.TypeT) = Complex{Real}
SymbolicUtils.promote_shape(::typeof(expim), ::SymbolicUtils.ShapeT) =
    SymbolicUtils.ShapeVecT()

Symbolics.@register_derivative expim(x) 1 im * expim_symbolic(x)

leading_sign(v) = 0
leading_sign(v::Num) = leading_sign(SymbolicUtils.unwrap(v))
leading_sign(v::Real) = isnan(v) ? 0 : Int(sign(v))
function leading_sign(u::RawExpression)
    if SymbolicUtils.iscall(u)
        op = SymbolicUtils.operation(u)
        args = SymbolicUtils.arguments(u)
        (op === (+) && !isempty(args)) && return leading_sign(args[1])
        if op === (*)
            for f in args
                s = const_sign(f)
                iszero(s) || return s
            end
            return 0
        end
    end
    return const_sign(u)
end

function scale_expanded(a::Num, n::Int)
    isone(n) && return a
    raw = SymbolicUtils.unwrap(a)
    is_call_of(raw, +) || return n * a
    terms = SymbolicUtils.arguments(raw)
    result = n * Num(first(terms))
    for i in 2:length(terms)
        result += n * Num(terms[i])
    end
    return result
end

function phase_coeff_expanded(a::Num)
    u = SymbolicUtils.unwrap(a)
    if is_real_const(u)
        v = real_const_float(u)
        return iszero(v) ? CNUM_ONE : native(cis(v))
    end
    SymbolicUtils.symtype(u) <: Real || return nonreal_phase_argument(a)
    angle, flipped = canonical_angle(a, 1)
    c = atom_coeff(expim_expanded(angle))
    return flipped ? conj(c) : c
end

function canonical_angle(argument::Num, n::Int)
    angle = scale_expanded(argument, n)
    leading_sign(angle) < 0 || return (angle, false)
    return (scale_expanded(angle, -1), true)
end

phase_coeff(x) = phase_coeff_expanded(Num(expand(x)))

@inline function is_phase_power(x::RawExpression)
    is_phase(x) && return true
    is_call_of(x, ^) || return false
    args = SymbolicUtils.arguments(x)
    return length(args) == 2 && is_phase(args[1]) && is_const_int(args[2])
end
@inline function phase_power(x::RawExpression)
    is_phase(x) && return (only(SymbolicUtils.arguments(x)), 1)
    args = SymbolicUtils.arguments(x)
    return (only(SymbolicUtils.arguments(args[1])), const_int(args[2]))
end

function split_phase_factors(factors)
    ordinary = RawExpression[]
    angle = RAW_ZERO
    phase_factors = 0
    changed = false
    for factor in factors
        if is_phase_power(factor)
            argument, exponent = phase_power(factor)
            angle += exponent * argument
            phase_factors += 1
            changed |= exponent != 1
        else
            push!(ordinary, factor)
        end
    end
    return (ordinary, angle, phase_factors, changed)
end

function rewrite_phase(x::RawExpression)
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if length(args) == 1 && is_phase(args[1])
        argument = only(SymbolicUtils.arguments(args[1]))
        op === conj && return expim_expanded(-argument)
        op === real && return cos(argument)
        op === imag && return sin(argument)
        (op === abs || op === abs2) && return RAW_ONE
    elseif op === (/) && length(args) == 2
        denominator = args[2]
        factors = is_call_of(denominator, *) ? SymbolicUtils.arguments(denominator) : (denominator,)
        any(is_phase_power, factors) || return x
        ordinary, angle, _, _ = split_phase_factors(factors)
        amplitude = isempty(ordinary) ? args[1] : args[1] / foldl(*, ordinary; init = RAW_ONE)
        return amplitude * expim_expanded(SymbolicUtils.expand(-angle))
    elseif op === (*)
        ordinary, angle, phase_factors, changed = split_phase_factors(args)
        (phase_factors >= 2 || changed) || return x
        expanded_angle = SymbolicUtils.expand(angle)
        no_phase = is_const_zero(expanded_angle)
        isempty(ordinary) && return no_phase ? RAW_ONE : expim_expanded(expanded_angle)
        amplitude = foldl(*, ordinary; init = RAW_ONE)
        no_phase && return amplitude
        return raw_complex(amplitude, 0 // 1) * expim_expanded(expanded_angle)
    elseif op === (^) && length(args) == 2 && is_phase(args[1])
        is_const_int(args[2]) || return x
        exponent = const_int(args[2])
        argument = only(SymbolicUtils.arguments(args[1]))
        iszero(exponent) && return RAW_ONE
        return expim_expanded(SymbolicUtils.expand(exponent * argument))
    end
    return x
end

const PHASE_NORMALIZER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Fixpoint(SymbolicUtils.Rewriters.Postwalk(rewrite_phase)),
)

normalize_phase(x::RawExpression) = PHASE_NORMALIZER(x)
function simplify_raw(x; kwargs...)
    normalized = normalize_phase(x)
    simplified = SymbolicUtils.simplify(normalized; kwargs...)
    return normalize_phase(simplified)
end

@noinline function canonical_phase_monomial!(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        syms::Vector{RawExpression},
        exps::Vector{Rational{Int}},
    ) where {E <: ExactScalar}
    flip = all(i -> !is_phase(syms[i]) || exps[i] < 0, eachindex(syms))
    angle = NUM_ZERO
    kept = 0
    @inbounds for i in eachindex(syms)
        symbol = syms[i]
        exponent = exps[i]
        if is_phase(symbol)
            denominator(exponent) == 1 ||
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            argument = Num(only(SymbolicUtils.arguments(symbol)))
            term = scale_expanded(argument, flip ? -numerator(exponent) : numerator(exponent))
            angle = iszero_num(angle) ? term : angle + term
        else
            kept += 1
            syms[kept] = symbol
            exps[kept] = exponent
        end
    end
    resize!(syms, kept)
    resize!(exps, kept)
    phase = phase_coeff_expanded(angle)
    if phase.tail isa Native
        unit = flip ? conj(native_scalar(phase)) : native_scalar(phase)
        return Monomial{E}(scalar_mul(scalar, as_tier(E, unit)), syms, exps)
    end
    phase_term = only(as_tier(E, phase.tail::Poly))
    append!(syms, phase_term.syms)
    append!(exps, flip ? -phase_term.exps : phase_term.exps)
    sort_factors!(syms, exps)
    return Monomial{E}(
        scalar_mul(scalar, monomial_scalar(phase_term)), syms, exps,
    )
end

@inline unit_modulus(z::ComplexF64) = abs2(z) == 1.0
@inline unit_modulus(z::ExactComplex) =
    widemul(z.re, z.re) + widemul(z.im, z.im) == widemul(z.den, z.den)
@inline unit_modulus(z::BigExactComplex) = isone(abs2(z))

function is_pure_phase(c::Coeff)
    tail = c.tail
    return tail isa Poly && is_pure_phase(tail)
end
function is_pure_phase(tail::Poly)
    length(tail.terms) == 1 || return false
    monomial = only(tail.terms)
    (length(monomial.syms) == 1 && is_phase(only(monomial.syms))) || return false
    denominator(only(monomial.exps)) == 1 || return false
    return unit_modulus(monomial_scalar(monomial))
end
