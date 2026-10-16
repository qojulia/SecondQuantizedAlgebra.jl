@inline factor_cnum(s::SymbolicUtils.BasicSymbolic, e::Rational{Int}) =
    tiered(monomial_terms, EXACT_ONE, RawExpression[s], Rational{Int}[e])

@inline function euler_cos(argument)
    phase = phase_coeff(argument)
    return CNUM_HALF * (phase + conj(phase))
end

@inline function euler_sin(argument)
    phase = phase_coeff(argument)
    difference = phase - conj(phase)
    return CNUM_NEG_IM * CNUM_HALF * difference
end

function exponential_leaf_phase(u::RawExpression, argument)
    real_part, imaginary_part = realimag(to_cnum(argument))
    (iszero_num(real_part) && SymbolicUtils.symtype(SymbolicUtils.unwrap(imaginary_part)) <: Real) ||
        return to_cnum(u)
    return phase_coeff(imaginary_part)
end

function exponential_tree(x)
    u = SymbolicUtils.unwrap(x)
    u isa RawExpression || return to_cnum(u)
    return fold_tree(exponential_leaf, u)
end

function exponential_leaf(u::RawExpression)
    SymbolicUtils.iscall(u) || return to_cnum(u)
    op = SymbolicUtils.operation(u)
    args = SymbolicUtils.arguments(u)
    length(args) == 1 || return to_cnum(u)
    argument = only(args)
    op === exp && return exponential_leaf_phase(u, argument)
    op === cis && return phase_coeff(argument)
    op === cos && return euler_cos(argument)
    op === sin && return euler_sin(argument)
    return to_cnum(u)
end

function exponential_factor(symbol::RawExpression, exponent::Rational{Int})
    (denominator(exponent) == 1 && (is_call_of(symbol, cos) || is_call_of(symbol, sin))) ||
        return factor_cnum(symbol, exponent)
    return pow_cnum_integer(exponential_leaf(symbol), numerator(exponent))
end

function map_factors(f::F, terms::Vector{<:Monomial}) where {F}
    result = CNUM_ZERO
    for m in terms
        term = to_cnum(monomial_scalar(m))
        @inbounds for i in eachindex(m.syms)
            term = term * f(m.syms[i], m.exps[i])
        end
        result = result + term
    end
    return result
end

function exponential_form(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    tail isa Poly && return map_factors(exponential_factor, tail.terms)
    return exponential_tree(tail.expr)
end

function trigonometric_factor(symbol::RawExpression, exponent::Rational{Int})
    (is_phase(symbol) && denominator(exponent) == 1) || return factor_cnum(symbol, exponent)
    angle, flipped = canonical_angle(Num(phase_argument(symbol)), numerator(exponent))
    sine = to_cnum(sin(angle))
    return to_cnum(cos(angle)) + (flipped ? CNUM_NEG_IM : CNUM_IM) * sine
end

@inline function rewrite_phase_to_trig(x)
    is_phase(x) || return x
    argument = phase_argument(x)
    return cos(argument) + im * sin(argument)
end

const TRIGONOMETRIC_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(rewrite_phase_to_trig),
)

function trigonometric_form(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    if tail isa RawSymbolicCoeff
        expanded = SymbolicUtils.expand(TRIGONOMETRIC_REWRITER(tail.expr))
        return raw_leaf(simplify_raw(expanded))
    end
    result = map_factors(trigonometric_factor, tail.terms)
    return result.tail isa RawSymbolicCoeff ? raw_leaf(normalize_phase(result.tail.expr)) : result
end

"""
    exponential_form(x)

Rewrite cosines and sines in `x` as unit phases.

`x` can be a coefficient or a whole operator expression. The function is not
exported.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: exponential_form

julia> @variables ω t;

julia> a = Destroy(FockSpace(:cavity), :a);

julia> exponential_form(cos(ω * t) * a)
((1//2)*exp(im*t*ω) + (1//2)*exp(-im*t*ω)) * a
```

See also [`trigonometric_form`](@ref), [`expim`](@ref).
"""
exponential_form(x::Num) = exponential_tree(x)
exponential_form(x::SymbolicUtils.BasicSymbolic) = exponential_tree(x)
exponential_form(x::Number) = x
exponential_form(x::Complex{Num}) = exponential_form(to_cnum(x))

"""
    PhaseTerm

One Fourier component of a coefficient, as returned by [`phase_terms`](@ref).

It stands for `amplitude * expim(phase)`, with the fields `amplitude`, a [`Coeff`](@ref),
and `phase`, a real symbolic expression.

See also [`phase_terms`](@ref).
"""
struct PhaseTerm
    amplitude::Coeff
    phase::Num
end

@inline constant_phase_term(c::Coeff) = PhaseTerm(c, NUM_ZERO)

@noinline no_phase_decomposition_error() = throw(
    ArgumentError("a phase occurs inside an operation with no finite phase decomposition"),
)
@noinline phase_denominator_error() = throw(
    ArgumentError("a phase-bearing denominator does not have a finite phase decomposition"),
)

function monomial_phase_term(m::Monomial)
    phase_index = phase_factor_index(m.syms)
    for i in eachindex(m.syms)
        i == phase_index && continue
        contains_phase(m.syms[i]) && no_phase_decomposition_error()
    end
    phase_index == 0 && return constant_phase_term(from_poly([m]))

    exponent = m.exps[phase_index]
    denominator(exponent) == 1 || fractional_phase_error()
    argument = phase_argument(m.syms[phase_index])
    phase = scale_expanded(Num(argument), numerator(exponent))
    amplitude = from_poly([replace_exps(m, phase_index, 0 // 1, 0, 0 // 1)])
    return PhaseTerm(amplitude, phase)
end

contains_phase(x) = has_subterm(is_phase, x)

function multiply_phase_terms(left::Vector{PhaseTerm}, right::Vector{PhaseTerm})
    result = PhaseTerm[]
    sizehint!(result, length(left) * length(right))
    for a in left, b in right
        push!(result, PhaseTerm(a.amplitude * b.amplitude, expand(a.phase + b.phase)))
    end
    return merge_phase_terms!(result)
end

function merge_phase_terms!(terms::Vector{PhaseTerm})
    slot = Dict{Num, Int}()
    kept = 0
    for term in terms
        k = get(slot, term.phase, 0)
        if iszero(k)
            kept += 1
            terms[kept] = term
            slot[term.phase] = kept
        else
            terms[k] = PhaseTerm(terms[k].amplitude + term.amplitude, term.phase)
        end
    end
    resize!(terms, kept)
    return filter!(term -> !iszero(term.amplitude), terms)
end

function phase_power_terms(base, exponent::Int)
    exponent < 0 && contains_phase(base) && phase_denominator_error()
    exponent < 0 && return PhaseTerm[constant_phase_term(pow_cnum_integer(to_cnum(base), exponent))]

    result = PhaseTerm[constant_phase_term(CNUM_ONE)]
    factor = raw_phase_terms(base)
    for _ in 1:exponent
        result = multiply_phase_terms(result, factor)
    end
    return result
end

raw_phase_terms(x::Number) = PhaseTerm[constant_phase_term(to_cnum(x))]
function raw_phase_terms(x::RawExpression)
    is_numeric_const(x) && return PhaseTerm[constant_phase_term(const_coeff(x))]
    if is_phase_power(x)
        argument, exponent = phase_power(x)
        return PhaseTerm[PhaseTerm(CNUM_ONE, scale_expanded(Num(argument), exponent))]
    end

    SymbolicUtils.iscall(x) || return PhaseTerm[constant_phase_term(to_cnum(x))]
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        result = PhaseTerm[]
        for argument in args
            append!(result, raw_phase_terms(argument))
        end
        return result
    elseif op === (*)
        result = PhaseTerm[constant_phase_term(CNUM_ONE)]
        for factor in args
            result = multiply_phase_terms(result, raw_phase_terms(factor))
        end
        return result
    elseif op === (/) && length(args) == 2
        contains_phase(args[2]) && phase_denominator_error()
        denominator_coeff = to_cnum(args[2])
        return PhaseTerm[PhaseTerm(term.amplitude / denominator_coeff, term.phase) for term in raw_phase_terms(args[1])]
    elseif op === (^) && length(args) == 2
        is_const_int(args[2]) && return phase_power_terms(args[1], const_int(args[2]))
        contains_phase(args[1]) && fractional_phase_error()
        return PhaseTerm[constant_phase_term(to_cnum(x))]
    elseif op === conj && length(args) == 1
        return PhaseTerm[PhaseTerm(conj(term.amplitude), -term.phase) for term in raw_phase_terms(only(args))]
    end
    contains_phase(x) && no_phase_decomposition_error()
    return PhaseTerm[constant_phase_term(to_cnum(x))]
end

"""
    phase_terms(c)

Split a coefficient into its Fourier components.

The coefficient is written as a sum of `amplitude * expim(phase)`, and each term
of that sum is returned with the fields `amplitude` and `phase`. Cosines and sines
have to be turned into phases first with [`exponential_form`](@ref). A coefficient
that is not a finite sum of phases throws an `ArgumentError`. The function is not
exported.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: phase_terms, exponential_form

julia> @variables ω t g;

julia> c = exponential_form(g * cos(ω * t));

julia> terms = phase_terms(c); length(terms)
2

julia> terms[2].amplitude, terms[2].phase
((1//2)*g, t*ω)
```

See also [`exponential_form`](@ref), [`expim`](@ref).
"""
function phase_terms(c::Coeff)
    normalized = exponential_form(c)
    tail = normalized.tail
    if tail isa Native
        return iszero(normalized) ? PhaseTerm[] : PhaseTerm[constant_phase_term(normalized)]
    elseif tail isa Poly
        return merge_phase_terms!(PhaseTerm[monomial_phase_term(monomial) for monomial in tail.terms])
    end
    return merge_phase_terms!(raw_phase_terms(tail.expr))
end

"""
    trigonometric_form(x)

Rewrite unit phases in `x` as cosines and sines.

Each phase [`expim`](@ref)`(θ)` becomes ``\\cos θ + i \\sin θ``, so opposite phases
combine into real functions. `x` can be a coefficient or a whole operator
expression. The function is not exported.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: trigonometric_form

julia> @variables ω t;

julia> a = Destroy(FockSpace(:cavity), :a);

julia> trigonometric_form((expim(ω * t) + expim(-ω * t)) * a)
2cos(t*ω) * a
```

See also [`exponential_form`](@ref), [`expim`](@ref).
"""
trigonometric_form(x::Union{Num, SymbolicUtils.BasicSymbolic, Complex{Num}}) = trigonometric_form(to_cnum(x))
trigonometric_form(x::Number) = x
