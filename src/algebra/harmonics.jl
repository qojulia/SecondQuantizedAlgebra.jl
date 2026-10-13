@noinline function harmonic_error(message::String)
    throw(ArgumentError(message))
end

function harmonic_frequency_coeff(ω::Real, t::Num)
    SymbolicUtils.issym(SymbolicUtils.unwrap(t)) ||
        harmonic_error("the time argument must be a symbolic variable; got `$t`")
    frequency = to_cnum(ω)
    iszero_cnum(frequency) && harmonic_error("the base frequency must be nonzero")
    coefficient_depends_on(frequency, SymbolicUtils.unwrap(t)) &&
        harmonic_error("the base frequency `$ω` depends on `$t`")
    return frequency
end

function exact_integer(c::CNum)
    c.tail isa Native || return false, 0
    z = native_scalar(c.slot)
    if z isa ComplexF64
        iszero(imag(z)) && isinteger(real(z)) && abs(real(z)) < 2.0^62 ||
            return false, 0
        return true, Int(real(z))
    end
    (iszero(z.im) && isone(z.den)) || return false, 0
    return true, Int(z.re)
end

# The order is read from the ratio of the phase rate to ω and then confirmed by exact
# subtraction, so a ratio the simplifier cannot reduce is refused rather than misassigned.
function harmonic_order(phase::Num, ω::Num, frequency::CNum, t::Num)
    variable = SymbolicUtils.unwrap(t)
    rate = expand(Symbolics.derivative(phase, t))
    raw_depends_on(SymbolicUtils.unwrap(rate), variable) &&
        harmonic_error("the phase `$phase` is not linear in `$t`")
    rate_coeff = to_cnum(rate)
    iszero_cnum(rate_coeff) && return 0
    ratio = div_cnum(rate_coeff, frequency)
    ok, n = exact_integer(ratio)
    ok || ((ok, n) = exact_integer(simplify_prefactor(ratio)))
    (ok && iszero_cnum(rate_coeff - n * frequency)) || harmonic_error(
        "the phase `$phase` is not an integer multiple of `$ω * $t`",
    )
    return n
end

function coefficient_harmonics(c::CNum, ω::Num, frequency::CNum, t::Num)
    out = Dict{Int, CNum}()
    variable = SymbolicUtils.unwrap(t)
    for term in phase_terms(c)
        coefficient_depends_on(term.amplitude, variable) && harmonic_error(
            "the amplitude `$(to_num(term.amplitude))` depends on `$t` outside a phase " *
                "`expim(n*ω*$t)`, so it is not a finite sum of harmonics",
        )
        n = harmonic_order(term.phase, ω, frequency, t)
        offset = expand(term.phase - n * ω * t)
        component = iszero(offset) ? term.amplitude : term.amplitude * expim(offset)
        total = add_cnum(get(out, n, CNUM_ZERO), component)
        iszero_cnum(total) ? delete!(out, n) : (out[n] = total)
    end
    return out
end

"""
    harmonics(x, ω, t) -> Dict{Int}

Decompose `x` into its Fourier components in the base frequency `ω`, so that

```julia
x == sum(xₙ * expim(n * ω * t) for (n, xₙ) in harmonics(x, ω, t))
```

after [`exponential_form`](@ref). `x` is a coefficient, a quantum expression, a
[`SuperOperator`](@ref) or a [`Lindbladian`](@ref); the components have the type of `x` (a
`Coeff` for a scalar, a `QAdd` for an operator). Zero components are omitted.

Trigonometric factors are converted to `expim` first, as in [`phase_terms`](@ref). Every
phase must be `n * ω * t + φ` with an exact integer `n` and a `t`-independent offset `φ`,
which stays in the component as `expim(φ)`. Anything else throws an `ArgumentError`
instead of being dropped or assigned to a guessed harmonic: a frequency that is not an exact
integer multiple of `ω` (such as `ω'` or `ω/2`), a phase nonlinear in `t`, or a
`t`-dependence outside a phase. Choose `ω` as the common base frequency of all phases.

For a `Lindbladian` the projection acts on the GKSL form, `H`, the rates and the jumps, and
returns a `Lindbladian`; a jump `J = Σ Jₘ expim(m ω t)` with partner `K = Σ Kₗ expim(l ω t)`
contributes `𝒟[Jₘ, Kₗ]` to harmonic `m - l`. Only the secular (`n = 0`) part is again
Hermiticity preserving in general.

See also [`harmonic`](@ref), [`secular`](@ref).
"""
harmonics(x::Union{Number, Coeff}, ω::Real, t::Num) =
    coefficient_harmonics(to_cnum(x), Num(ω), harmonic_frequency_coeff(ω, t), t)
harmonics(op::QSym, ω::Real, t::Num) =
    harmonics(single_qadd(CNUM_ONE, Op[op]), ω, t)

function harmonics(q::QAdd, ω::Real, t::Num)
    frequency = harmonic_frequency_coeff(ω, t)
    parts = Dict{Int, QTermDict}()
    for (term, c) in q.arguments, (n, cn) in coefficient_harmonics(c, Num(ω), frequency, t)
        addto_key!(get!(QTermDict, parts, n), copy_key(term), cn)
    end
    return Dict{Int, QAdd}(n => QAdd(d, drop_unused_indices(d, q.indices)) for (n, d) in parts)
end

"""
    harmonic(x, ω, t, n::Integer)

The coefficient of `expim(n * ω * t)` in `x`, with the same type as `x` (a `Coeff` for a
scalar). Every phase of `x` must be an exact integer multiple of `ω * t` plus a constant; see
[`harmonics`](@ref) for the precise contract. A component that is not of this form throws
an `ArgumentError` even when `n` is a different order, because it could not be excluded
from harmonic `n` without assuming a relation between frequencies.

# Examples

```julia
@variables ω0::Real ωd::Real F::Real t::Real
@qnumbers a::Destroy(FockSpace(:a))
H = ω0 * a' * a - F * (a + a') * (expim(ωd * t) + expim(-ωd * t)) / 2
Ht = transform(H, Rotation(a, ωd * t, t))
secular(Ht, ωd, t)      # (ω0 - ωd) a'a - F/2 (a + a')
harmonic(Ht, ωd, t, 2)  # -F/2 a'
```

See also [`secular`](@ref), [`phase_terms`](@ref).
"""
harmonic(x::Union{Number, Coeff}, ω::Real, t::Num, n::Integer) =
    get(harmonics(x, ω, t), Int(n), CNUM_ZERO)
harmonic(op::QSym, ω::Real, t::Num, n::Integer) =
    harmonic(single_qadd(CNUM_ONE, Op[op]), ω, t, n)
function harmonic(q::QAdd, ω::Real, t::Num, n::Integer)
    frequency = harmonic_frequency_coeff(ω, t)
    order = Int(n)
    return map_coefficients(c -> get(coefficient_harmonics(c, Num(ω), frequency, t), order, CNUM_ZERO), q)
end
