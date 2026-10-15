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

function negate_expanded(a::Num)
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result -= Num(term)
        end
        return result
    end
    return -a
end

function scale_expanded(a::Num, n::Int)
    isone(n) && return a
    n == -1 && return negate_expanded(a)
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result += n * Num(term)
        end
        return result
    end
    return n * a
end

function phase_coeff_expanded(a::Num)
    u = SymbolicUtils.unwrap(a)
    if is_real_const(u)
        v = real_const_float(u)
        return iszero(v) ? CNUM_ONE : native(cis(v))
    end
    SymbolicUtils.symtype(u) <: Real || return nonreal_phase_argument(a)
    neg = leading_sign(a) < 0
    c = atom_coeff(expim_expanded(neg ? negate_expanded(a) : a))
    return neg ? conj_cnum(c) : c
end

phase_coeff(x) = phase_coeff_expanded(Num(expand(x)))

@inline function is_phase_power(x::RawExpression)
    is_phase(x) && return true
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === (^) || return false
    args = SymbolicUtils.arguments(x)
    return length(args) == 2 && is_phase(args[1]) && is_const_int(args[2])
end
@inline function phase_power(x::RawExpression)
    is_phase(x) && return (only(SymbolicUtils.arguments(x)), 1)
    args = SymbolicUtils.arguments(x)
    return (only(SymbolicUtils.arguments(args[1])), const_int(args[2]))
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
        angle = RAW_ZERO
        found_phase = false
        denominator = args[2]
        if is_phase_power(denominator)
            argument, exponent = phase_power(denominator)
            return args[1] * expim_expanded(SymbolicUtils.expand(-exponent * argument))
        elseif SymbolicUtils.iscall(denominator) && SymbolicUtils.operation(denominator) === (*)
            factors = SymbolicUtils.arguments(denominator)
            ordinary_count = 0
            for factor in factors
                if is_phase_power(factor)
                    argument, exponent = phase_power(factor)
                    angle += exponent * argument
                    found_phase = true
                else
                    ordinary_count += 1
                end
            end
            found_phase || return x
            if ordinary_count == 0
                return args[1] * expim_expanded(SymbolicUtils.expand(-angle))
            end
            ordinary_denominator = RAW_ONE
            for factor in factors
                is_phase_power(factor) ||
                    (ordinary_denominator = ordinary_denominator * factor)
            end
            return (args[1] / ordinary_denominator) *
                expim_expanded(SymbolicUtils.expand(-angle))
        end
    elseif op === (*)
        ordinary = RawExpression[]
        angle = RAW_ZERO
        phase_factors = 0
        changed = false
        for factor in args
            if !is_phase_power(factor)
                push!(ordinary, factor)
                continue
            end
            argument, exponent = phase_power(factor)
            angle += exponent * argument
            phase_factors += 1
            changed |= exponent != 1
        end
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

@noinline function canonical_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        syms::Vector{RawExpression},
        exps::Vector{Rational{Int}},
    ) where {E <: ExactScalar}
    ordinary_syms = RawExpression[]
    ordinary_exps = Rational{Int}[]
    sizehint!(ordinary_syms, length(syms) - 1)
    sizehint!(ordinary_exps, length(exps) - 1)
    angle = NUM_ZERO
    @inbounds for i in eachindex(syms)
        symbol = syms[i]
        exponent = exps[i]
        if is_phase(symbol)
            denominator(exponent) == 1 ||
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            argument = Num(only(SymbolicUtils.arguments(symbol)))
            angle += scale_expanded(argument, numerator(exponent))
        else
            push!(ordinary_syms, symbol)
            push!(ordinary_exps, exponent)
        end
    end
    ordinary_syms, ordinary_exps = merge_factor_list(ordinary_syms, ordinary_exps)
    phase = phase_coeff_expanded(angle)
    phase_tail = phase.tail
    if phase_tail isa Native
        return Monomial{E}(scalar_mul(scalar, as_tier(E, native_scalar(phase))), ordinary_syms, ordinary_exps)
    end
    phase_tail isa Poly || throw(ArgumentError("expected a polynomial coefficient"))
    phase_term = only(as_tier(E, phase_tail))
    append!(ordinary_syms, phase_term.syms)
    append!(ordinary_exps, phase_term.exps)
    sort_factors!(ordinary_syms, ordinary_exps)
    return Monomial{E}(
        scalar_mul(scalar, monomial_scalar(phase_term)), ordinary_syms, ordinary_exps,
    )
end

@noinline function scaled_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        symbol::Num,
        exponent::Rational{Int},
    ) where {E <: ExactScalar}
    denominator(exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(symbol))))
    phase = phase_coeff_expanded(scale_expanded(argument, numerator(exponent)))
    phase_tail = phase.tail
    phase_tail isa Native && return constant_term(E, scalar_mul(scalar, as_tier(E, native_scalar(phase))))
    phase_tail isa Poly || throw(ArgumentError("expected a polynomial coefficient"))
    phase_term = only(as_tier(E, phase_tail))
    return Monomial{E}(
        scalar_mul(scalar, monomial_scalar(phase_term)),
        phase_term.syms,
        phase_term.exps,
    )
end

@noinline function merged_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        left::Num,
        left_exponent::Rational{Int},
        right::Num,
        right_exponent::Rational{Int},
    ) where {E <: ExactScalar}
    denominator(left_exponent) == denominator(right_exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    left_argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(left))))
    right_argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(right))))
    left_power = numerator(left_exponent)
    right_power = numerator(right_exponent)
    common_negative = left_power < 0 && right_power < 0
    if common_negative
        left_power = -left_power
        right_power = -right_power
    end
    angle = scale_expanded(left_argument, left_power) +
        scale_expanded(right_argument, right_power)
    phase = phase_coeff_expanded(angle)
    phase_tail = phase.tail
    phase_tail isa Native && return constant_term(E, scalar_mul(scalar, as_tier(E, native_scalar(phase))))
    phase_tail isa Poly || throw(ArgumentError("expected a polynomial coefficient"))
    phase_term = only(as_tier(E, phase_tail))
    exponents = common_negative ? -phase_term.exps : phase_term.exps
    return Monomial{E}(
        scalar_mul(scalar, monomial_scalar(phase_term)),
        phase_term.syms,
        exponents,
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

function pure_phase_parts(tail::Poly)
    monomial = only(tail.terms)
    symbol = only(monomial.syms)
    exponent = only(monomial.exps)
    scalar = monomial_scalar(monomial)
    argument = Num(only(SymbolicUtils.arguments(symbol)))
    angle = scale_expanded(argument, numerator(exponent))
    flip_sine = false
    if leading_sign(angle) < 0
        angle = negate_expanded(angle)
        flip_sine = true
    end
    return (scalar, angle, flip_sine)
end

@inline function scaled_trig(a::Real, trig, angle::Num)
    iszero(a) && return NUM_ZERO
    isone(a) && return trig(angle)
    a == -1 && return -trig(angle)
    return num_from_scalar(a) * trig(angle)
end

function pure_phase_real(tail::Poly)
    scalar, angle, flip_sine = pure_phase_parts(tail)
    sine = flip_sine ? -imag(scalar) : imag(scalar)
    return scaled_trig(real(scalar), cos, angle) - scaled_trig(sine, sin, angle)
end

function pure_phase_imag(tail::Poly)
    scalar, angle, flip_sine = pure_phase_parts(tail)
    sine = flip_sine ? -real(scalar) : real(scalar)
    return scaled_trig(sine, sin, angle) + scaled_trig(imag(scalar), cos, angle)
end
