# === Fock transformations ===

"""
    Displace(a, α[, t])

Displace a Fock mode by the scalar amplitude `α`, so `a ↦ a + α`. The timed form
stores the complete c-number gauge of the moving displacement.
"""
function Displace(a::Op, α::Coefficient)
    d = _fock_or_throw(a, "`Displace`")
    return _fock_displacement(d, _to_cnum(α))
end

function _fock_displacement(d::Op, c::CNum)
    return _static_transform(
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (c, Op[]))),
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (_neg_cnum(c), Op[]))),
    )
end

function Displace(a::Op, α::Coefficient, t::Num)
    d = _fock_or_throw(a, "`Displace`")
    tt = _time_or_throw(t)
    c = _to_cnum(α)
    derivative = _dt(c, tt)
    derivative_adjoint = _conj_cnum(derivative)
    gauge = _fock_displacement_gauge(
        d, c, derivative, derivative_adjoint,
        _mul_cnum(_CNUM_NEG_IM, derivative),
        _mul_cnum(_CNUM_IM, derivative_adjoint),
    )
    return _timed_transform(_fock_displacement(d, c), gauge, tt)
end

function _fock_displacement_gauge(
        d::Op, c::CNum, derivative::CNum, derivative_adjoint::CNum,
        raising_coefficient::CNum, lowering_coefficient::CNum,
    )
    scalar = _mul_cnum(
        _CNUM_NEG_IM,
        _mul_cnum(
            _CNUM_HALF,
            _add_cnum(
                _mul_cnum(_conj_cnum(c), derivative),
                _neg_cnum(_mul_cnum(c, derivative_adjoint)),
            ),
        ),
    )
    return _rule_qadd(
        (raising_coefficient, Op[adjoint(d)]),
        (lowering_coefficient, Op[d]),
        (scalar, Op[]),
    )
end

struct _FockLinearReference
    frequency::CNum
    drive::CNum
    drive_adjoint::CNum
    normalized_drive::CNum
end

function _linear_fock_hamiltonian(H::QAdd, d::Op)
    isempty(H.indices) || _unitary_error(
        "automatic `Displace` does not support summed reference Hamiltonians",
    )
    frequency = _CNUM_ZERO
    drive = _CNUM_ZERO
    drive_adjoint = _CNUM_ZERO
    scalar = _CNUM_ZERO
    raising = adjoint(d)
    for (term, coefficient) in H
        isempty(term.ne) || _unitary_error(
            "automatic `Displace` does not support constrained reference terms",
        )
        if isempty(term.ops)
            scalar = _add_cnum(scalar, coefficient)
        elseif term.ops == Op[raising, d]
            frequency = _add_cnum(frequency, coefficient)
        elseif term.ops == Op[raising]
            drive = _add_cnum(drive, coefficient)
        elseif term.ops == Op[d]
            drive_adjoint = _add_cnum(drive_adjoint, coefficient)
        else
            unsupported = _rule_qadd((coefficient, copy(term.ops)))
            _unitary_error(
                "automatic `Displace` expects only scalar, `$raising * $d`, `$raising`, " *
                    "and `$d` terms; got `$(sprint(show, unsupported))`",
            )
        end
    end

    isequal(frequency, _conj_cnum(frequency)) ||
        _unitary_error(
            "automatic `Displace` requires a real oscillator frequency; " *
                "got `$(to_num(frequency))`",
        )
    isequal(scalar, _conj_cnum(scalar)) || _unitary_error(
        "automatic `Displace` requires a Hermitian scalar reference term",
    )
    normalized_drive = exponential_form(drive)
    normalized_adjoint = exponential_form(drive_adjoint)
    isequal(normalized_adjoint, _conj_cnum(normalized_drive)) || _unitary_error(
        "automatic `Displace` requires adjoint coefficients for `$raising` and `$d`",
    )
    return _FockLinearReference(
        frequency, drive, drive_adjoint, normalized_drive,
    )
end

function _validate_fock_time(data::_FockLinearReference, t::Num)
    variable = SymbolicUtils.unwrap(t)
    _coefficient_depends_on(data.frequency, variable) && _unitary_error(
        "automatic `Displace` requires a time-independent oscillator frequency; " *
            "got `$(to_num(data.frequency))`",
    )
    return nothing
end

function _harmonic_frequency(monomial::Monomial, t::Num)
    phase_index = 0
    variable = SymbolicUtils.unwrap(t)
    @inbounds for i in eachindex(monomial.syms)
        factor = monomial.syms[i]
        if _is_phase(factor)
            phase_index == 0 || _unitary_error(
                "automatic `Displace` found more than one phase in one drive component",
            )
            phase_index = i
        elseif _raw_depends_on(factor, variable)
            _unitary_error(
                "automatic `Displace` supports finite harmonic drives, not the " *
                    "time-dependent envelope `$(Num(factor))`",
            )
        end
    end
    phase_index == 0 && return _CNUM_ZERO

    exponent = monomial.exps[phase_index]
    denominator(exponent) == 1 || _unitary_error(
        "automatic `Displace` requires integer phase powers; got exponent `$exponent`",
    )
    phase = monomial.syms[phase_index]
    argument = Num(only(SymbolicUtils.arguments(phase))) * numerator(exponent)
    frequency = Symbolics.derivative(argument, t)
    _raw_depends_on(SymbolicUtils.unwrap(frequency), variable) &&
        _unitary_error(
            "automatic `Displace` requires a phase linear in `$t`; " *
                "got `$(Num(only(SymbolicUtils.arguments(phase))))`",
        )
    return _to_cnum(frequency)
end

function _bounded_harmonic_displacement(drive::CNum, frequency::CNum, t::Num)
    _iszero_cnum(drive) && return _CNUM_ZERO
    tail = drive.tail
    if tail isa Native
        _iszero_cnum(frequency) && _unitary_error(
            "automatic `Displace` found a resonant constant drive with zero frequency",
        )
        return _neg_cnum(drive) / frequency
    elseif !(tail isa Poly)
        _unitary_error(
            "automatic `Displace` supports only finite harmonic drive coefficients; " *
                "got `$(to_num(drive))`",
        )
    end

    amplitude = _CNUM_ZERO
    for monomial in tail.terms
        harmonic = _harmonic_frequency(monomial, t)
        divisor = _add_cnum(frequency, harmonic)
        component = _from_poly(Monomial[monomial])
        _iszero_cnum(divisor) && _unitary_error(
            "automatic `Displace` found a resonant drive component `$(to_num(component))`",
        )
        amplitude = _add_cnum(amplitude, _neg_cnum(component) / divisor)
    end
    return amplitude
end

function _static_fock_displacement(data::_FockLinearReference)
    _iszero_cnum(data.drive) && return _CNUM_ZERO
    _iszero_cnum(data.frequency) && _unitary_error(
        "automatic `Displace` found a resonant constant drive with zero frequency",
    )
    return _neg_cnum(data.drive) / data.frequency
end

"""
    Displace(a, Hlin)

Construct the static equilibrium displacement generated by the one-mode linear reference
Hamiltonian `Hlin = ω*a'a + η*a' + conj(η)*a + c`. Symbolic coefficients are treated as
time-independent parameters. A structurally zero response divisor is rejected.
"""
function Displace(a::Op, Hlin::QAdd)
    d = _fock_or_throw(a, "`Displace`")
    data = _linear_fock_hamiltonian(Hlin, d)
    return _fock_displacement(d, _static_fock_displacement(data))
end

"""
    Displace(a, Hlin, t)

Construct the bounded moving displacement generated by the one-mode linear reference
Hamiltonian `Hlin = ω*a'a + η(t)*a' + conj(η(t))*a + c(t)`. The drive must be a finite
sum of exact harmonic phases. The free homogeneous solution is set to zero; use
`Displace(a, α, t)` to supply a transient or arbitrary displacement field explicitly.

A structurally resonant harmonic is rejected. Symbolic divisors are retained as exact
quotients and therefore describe the response away from their resonance surfaces.
"""
function Displace(a::Op, Hlin::QAdd, t::Num)
    d = _fock_or_throw(a, "`Displace`")
    tt = _time_or_throw(t)
    data = _linear_fock_hamiltonian(Hlin, d)
    _validate_fock_time(data, tt)
    amplitude = _bounded_harmonic_displacement(
        data.normalized_drive, data.frequency, tt,
    )

    conjugate_amplitude = _conj_cnum(amplitude)
    derivative = _mul_cnum(
        _CNUM_NEG_IM,
        _add_cnum(_mul_cnum(data.frequency, amplitude), data.normalized_drive),
    )
    derivative_adjoint = _mul_cnum(
        _CNUM_IM,
        _add_cnum(
            _mul_cnum(data.frequency, conjugate_amplitude),
            _conj_cnum(data.normalized_drive),
        ),
    )
    gauge = _fock_displacement_gauge(
        d, amplitude, derivative, derivative_adjoint,
        _neg_cnum(_add_cnum(data.drive, _mul_cnum(data.frequency, amplitude))),
        _neg_cnum(
            _add_cnum(
                data.drive_adjoint,
                _mul_cnum(data.frequency, conjugate_amplitude),
            ),
        ),
    )
    return _timed_transform(_fock_displacement(d, amplitude), gauge, tt)
end

"""
    Rotation(a, θ)
    Rotation(a, θ, t)

Rotate a Fock mode by `exp(-im*θ*a'a)`, so `a ↦ exp(-im*θ)*a`. The positional timed
form requires a symbolic moving angle and a symbolic time variable.
"""
function Rotation(a::Op, θ::Real)
    d = _fock_or_throw(a, "`Rotation`")
    return _static_transform(
        _with_adjoint(d, _scaled(_conj_phase(θ), d)),
        _with_adjoint(d, _scaled(_phase(θ), d)),
    )
end

function Rotation(a::Op, θ::Num, t::Num)
    d = _fock_or_throw(a, "`Rotation`")
    tt = _time_or_throw(t)
    return _timed_transform(Rotation(a, θ), _gauge(adjoint(d) * d, θ, tt), tt)
end

"""
    Squeeze(a, r, ϕ = 0)
    Squeeze(a, r, ϕ, t)

Single-mode squeezing with `a ↦ cosh(r)*a + exp(im*ϕ)*sinh(r)*a'`. The timed form
supports both a moving magnitude and a moving phase.
"""
function Squeeze(a::Op, r::Real, ϕ::Real = 0)
    d = _fock_or_throw(a, "`Squeeze`")
    ch = _to_cnum(cosh(r))
    sh = _mul_cnum(_phase(ϕ), _to_cnum(sinh(r)))
    return _static_transform(
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (sh, Op[adjoint(d)]))),
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (_neg_cnum(sh), Op[adjoint(d)]))),
        ParamRelation[_hyp_rel(r)],
    )
end

function _squeeze_gauge(d::Op, r::Real, ϕ::Real, t::Num)
    sh = _to_cnum(sinh(r))
    ch = _to_cnum(cosh(r))
    imaginary_dr = _mul_cnum(_CNUM_IM, _dt(r, t))
    dϕ = _dt(ϕ, t)
    diagonal = _mul_cnum(dϕ, _mul_cnum(sh, sh))
    mixing = _mul_cnum(dϕ, _mul_cnum(sh, ch))
    return _rule_qadd(
        (diagonal, Op[adjoint(d), d]),
        (_mul_cnum(_CNUM_HALF, diagonal), Op[]),
        (
            _mul_cnum(
                _CNUM_HALF,
                _mul_cnum(_phase(ϕ), _add_cnum(mixing, _neg_cnum(imaginary_dr))),
            ),
            Op[adjoint(d), adjoint(d)],
        ),
        (
            _mul_cnum(
                _CNUM_HALF,
                _mul_cnum(_conj_phase(ϕ), _add_cnum(mixing, imaginary_dr)),
            ),
            Op[d, d],
        ),
    )
end

function Squeeze(a::Op, r::Real, ϕ::Real, t::Num)
    d = _fock_or_throw(a, "`Squeeze`")
    tt = _time_or_throw(t)
    return _timed_transform(Squeeze(a, r, ϕ), _squeeze_gauge(d, r, ϕ, tt), tt)
end

# === Two-mode and phase-space transformations ===

function _two_modes(a::Op, b::Op, what::AbstractString)
    x = _fock_or_throw(a, what)
    y = _fock_or_throw(b, what)
    site_key(x) == site_key(y) && _unitary_error("$what needs two distinct modes")
    return x, y
end

function _phase_pair(x::Op, p::Op, what::AbstractString)
    (is_position(x) && is_momentum(p)) || _unitary_error(
        "$what expects a `(Position, Momentum)` pair in that order",
    )
    site_key(x) == site_key(p) || _unitary_error(
        "$what expects both quadratures on the same site",
    )
    return nothing
end

function _beamsplitter(a::Op, b::Op, θ::Real)
    x, y = _two_modes(a, b, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (s, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (negative_s, Op[x]))),
        ),
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (negative_s, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (s, Op[x]))),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function _quadrature_rotation(x::Op, p::Op, θ::Real)
    _phase_pair(x, p, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        _pair_rules(
            x, p, _rule_qadd((c, Op[x]), (s, Op[p])),
            _rule_qadd((c, Op[p]), (negative_s, Op[x])),
        ),
        _pair_rules(
            x, p, _rule_qadd((c, Op[x]), (negative_s, Op[p])),
            _rule_qadd((c, Op[p]), (s, Op[x])),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

"""Mix two Fock modes (beamsplitter) or rotate a canonical quadrature pair."""
Rotation(a::Op, b::Op, θ::Real) =
    _is_phase_space(a) ? _quadrature_rotation(a, b, θ) : _beamsplitter(a, b, θ)

function Rotation(a::Op, b::Op, θ::Real, t::Num)
    tt = _time_or_throw(t)
    if _is_phase_space(a)
        U = _quadrature_rotation(a, b, θ)
        generator = (a * a + b * b) * (1 // 2)
    else
        x, y = _two_modes(a, b, "`Rotation`")
        U = _beamsplitter(x, y, θ)
        generator = im * (adjoint(x) * y - adjoint(y) * x)
    end
    return _timed_transform(U, _gauge(generator, θ, tt), tt)
end

function _two_mode_squeeze(a::Op, b::Op, r::Real)
    x, y = _two_modes(a, b, "`Squeeze`")
    u = _to_cnum(cosh(r))
    v = _to_cnum(sinh(r))
    negative_v = _neg_cnum(v)
    return _static_transform(
        merge(
            _with_adjoint(x, _rule_qadd((u, Op[x]), (v, Op[adjoint(y)]))),
            _with_adjoint(y, _rule_qadd((u, Op[y]), (v, Op[adjoint(x)]))),
        ),
        merge(
            _with_adjoint(x, _rule_qadd((u, Op[x]), (negative_v, Op[adjoint(y)]))),
            _with_adjoint(y, _rule_qadd((u, Op[y]), (negative_v, Op[adjoint(x)]))),
        ),
        ParamRelation[_hyp_rel(r)],
    )
end

function _quadrature_squeeze(x::Op, p::Op, r::Real)
    _phase_pair(x, p, "`Squeeze`")
    up = _to_cnum(exp(r))
    down = _to_cnum(inv(exp(r)))
    return _static_transform(
        _pair_rules(x, p, _scaled(up, x), _scaled(down, p)),
        _pair_rules(x, p, _scaled(down, x), _scaled(up, p)),
    )
end

"""Squeeze two Fock modes or a canonical quadrature pair."""
Squeeze(a::Op, b::Op, r::Real) =
    _is_phase_space(a) ? _quadrature_squeeze(a, b, r) : _two_mode_squeeze(a, b, r)

function Squeeze(a::Op, b::Op, r::Real, t::Num)
    tt = _time_or_throw(t)
    if _is_phase_space(a)
        U = _quadrature_squeeze(a, b, r)
        generator = (a * b + b * a) * (1 // 2)
    else
        x, y = _two_modes(a, b, "`Squeeze`")
        U = _two_mode_squeeze(x, y, r)
        generator = im * (adjoint(x) * adjoint(y) - y * x)
    end
    return _timed_transform(U, _gauge(generator, r, tt), tt)
end

function _quadrature_displacement(x::Op, p::Op, cx::CNum, cp::CNum)
    return _static_transform(
        _pair_rules(
            x, p, _rule_qadd((_CNUM_ONE, Op[x]), (cx, Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (cp, Op[])),
        ),
        _pair_rules(
            x, p, _rule_qadd((_CNUM_ONE, Op[x]), (_neg_cnum(cx), Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (_neg_cnum(cp), Op[])),
        ),
    )
end

"""Displace a canonical quadrature pair by real scalar shifts."""
function Displace(x::Op, p::Op, dx::Real, dp::Real)
    _phase_pair(x, p, "`Displace`")
    return _quadrature_displacement(x, p, _to_cnum(dx), _to_cnum(dp))
end

function _quadrature_displacement_gauge(
        x::Op, p::Op, cx::CNum, cp::CNum,
        derivative_x::CNum, derivative_p::CNum,
        x_coefficient::CNum, p_coefficient::CNum,
    )
    scalar = _mul_cnum(
        _neg_cnum(_CNUM_HALF),
        _add_cnum(
            _mul_cnum(cp, derivative_x),
            _neg_cnum(_mul_cnum(cx, derivative_p)),
        ),
    )
    return _rule_qadd(
        (x_coefficient, Op[x]), (p_coefficient, Op[p]), (scalar, Op[]),
    )
end

function Displace(x::Op, p::Op, dx::Real, dp::Real, t::Num)
    _phase_pair(x, p, "`Displace`")
    tt = _time_or_throw(t)
    cx = _to_cnum(dx)
    cp = _to_cnum(dp)
    derivative_x = _dt(cx, tt)
    derivative_p = _dt(cp, tt)
    gauge = _quadrature_displacement_gauge(
        x, p, cx, cp, derivative_x, derivative_p,
        derivative_p, _neg_cnum(derivative_x),
    )
    return _timed_transform(_quadrature_displacement(x, p, cx, cp), gauge, tt)
end

struct _QuadratureLinearReference
    quadratic_x::CNum
    quadratic_cross::CNum
    quadratic_p::CNum
    drive_x::CNum
    drive_p::CNum
    normalized_drive_x::CNum
    normalized_drive_p::CNum
end

@noinline function _quadrature_reference_error(coefficient::CNum, name::AbstractString)
    _unitary_error(
        "automatic `Displace` requires a real $name coefficient; " *
            "got `$(to_num(coefficient))`",
    )
end

function _linear_quadrature_hamiltonian(H::QAdd, x::Op, p::Op)
    isempty(H.indices) || _unitary_error(
        "automatic `Displace` does not support summed reference Hamiltonians",
    )
    half_quadratic_x = _CNUM_ZERO
    quadratic_cross = _CNUM_ZERO
    half_quadratic_p = _CNUM_ZERO
    drive_x = _CNUM_ZERO
    drive_p = _CNUM_ZERO
    for (term, coefficient) in H
        isempty(term.ne) || _unitary_error(
            "automatic `Displace` does not support constrained reference terms",
        )
        if isempty(term.ops)
            continue
        elseif term.ops == Op[x, x]
            half_quadratic_x = _add_cnum(half_quadratic_x, coefficient)
        elseif term.ops == Op[x, p]
            quadratic_cross = _add_cnum(quadratic_cross, coefficient)
        elseif term.ops == Op[p, p]
            half_quadratic_p = _add_cnum(half_quadratic_p, coefficient)
        elseif term.ops == Op[x]
            drive_x = _add_cnum(drive_x, coefficient)
        elseif term.ops == Op[p]
            drive_p = _add_cnum(drive_p, coefficient)
        else
            unsupported = _rule_qadd((coefficient, copy(term.ops)))
            _unitary_error(
                "automatic `Displace` expects only scalar, `$x * $x`, `$x * $p`, " *
                    "`$p * $p`, `$x`, and `$p` terms; got " *
                    "`$(sprint(show, unsupported))`",
            )
        end
    end

    isequal(H, H') || _unitary_error(
        "automatic `Displace` requires a Hermitian quadrature reference Hamiltonian",
    )
    quadratic_x = _add_cnum(half_quadratic_x, half_quadratic_x)
    quadratic_p = _add_cnum(half_quadratic_p, half_quadratic_p)
    for (coefficient, name) in (
            (quadratic_x, "x²"),
            (quadratic_cross, "symmetrized x-p"),
            (quadratic_p, "p²"),
        )
        isequal(coefficient, _conj_cnum(coefficient)) ||
            _quadrature_reference_error(coefficient, name)
    end
    normalized_drive_x = exponential_form(drive_x)
    normalized_drive_p = exponential_form(drive_p)
    isequal(normalized_drive_x, _conj_cnum(normalized_drive_x)) ||
        _quadrature_reference_error(normalized_drive_x, "x-drive")
    isequal(normalized_drive_p, _conj_cnum(normalized_drive_p)) ||
        _quadrature_reference_error(normalized_drive_p, "p-drive")
    return _QuadratureLinearReference(
        quadratic_x, quadratic_cross, quadratic_p,
        drive_x, drive_p, normalized_drive_x, normalized_drive_p,
    )
end

function _validate_quadrature_time(data::_QuadratureLinearReference, t::Num)
    variable = SymbolicUtils.unwrap(t)
    for (coefficient, name) in (
            (data.quadratic_x, "x²"),
            (data.quadratic_cross, "symmetrized x-p"),
            (data.quadratic_p, "p²"),
        )
        _coefficient_depends_on(coefficient, variable) && _unitary_error(
            "automatic `Displace` requires a time-independent $name coefficient; " *
                "got `$(to_num(coefficient))`",
        )
    end
    return nothing
end

@noinline function _quadrature_resonance_error(
        drive_x::CNum, drive_p::CNum, harmonic::CNum,
    )
    _unitary_error(
        "automatic `Displace` found a resonant quadrature drive component " *
            "`($(to_num(drive_x)), $(to_num(drive_p)))` at frequency " *
            "`$(to_num(harmonic))`",
    )
end

function _quadrature_response_component(
        data::_QuadratureLinearReference, drive_x::CNum, drive_p::CNum,
        harmonic::CNum,
    )
    imaginary_harmonic = _mul_cnum(_CNUM_IM, harmonic)
    upper_right = _add_cnum(data.quadratic_cross, imaginary_harmonic)
    lower_left = _add_cnum(
        data.quadratic_cross, _neg_cnum(imaginary_harmonic),
    )
    determinant = _add_cnum(
        _mul_cnum(data.quadratic_x, data.quadratic_p),
        _neg_cnum(_mul_cnum(upper_right, lower_left)),
    )
    _iszero_cnum(determinant) &&
        _quadrature_resonance_error(drive_x, drive_p, harmonic)

    displacement_x = _neg_cnum(
        _add_cnum(
            _mul_cnum(data.quadratic_p, drive_x),
            _neg_cnum(_mul_cnum(upper_right, drive_p)),
        ),
    ) / determinant
    displacement_p = _add_cnum(
        _mul_cnum(lower_left, drive_x),
        _neg_cnum(_mul_cnum(data.quadratic_x, drive_p)),
    ) / determinant
    return (x = displacement_x, p = displacement_p)
end

function _static_quadrature_displacement(data::_QuadratureLinearReference)
    (_iszero_cnum(data.drive_x) && _iszero_cnum(data.drive_p)) &&
        return (x = _CNUM_ZERO, p = _CNUM_ZERO)
    return _quadrature_response_component(
        data, data.drive_x, data.drive_p, _CNUM_ZERO,
    )
end

function _bounded_quadrature_drive(
        drive::CNum, data::_QuadratureLinearReference, t::Num, ::Val{axis},
    ) where {axis}
    _iszero_cnum(drive) && return (x = _CNUM_ZERO, p = _CNUM_ZERO)
    tail = drive.tail
    displacement_x = _CNUM_ZERO
    displacement_p = _CNUM_ZERO
    if tail isa Native
        response = axis === :x ?
            _quadrature_response_component(data, drive, _CNUM_ZERO, _CNUM_ZERO) :
            _quadrature_response_component(data, _CNUM_ZERO, drive, _CNUM_ZERO)
        return response
    elseif !(tail isa Poly)
        _unitary_error(
            "automatic `Displace` supports only finite harmonic drive coefficients; " *
                "got `$(to_num(drive))`",
        )
    end

    for monomial in tail.terms
        harmonic = _harmonic_frequency(monomial, t)
        component = _from_poly(Monomial[monomial])
        response = axis === :x ?
            _quadrature_response_component(data, component, _CNUM_ZERO, harmonic) :
            _quadrature_response_component(data, _CNUM_ZERO, component, harmonic)
        displacement_x = _add_cnum(displacement_x, response.x)
        displacement_p = _add_cnum(displacement_p, response.p)
    end
    return (x = displacement_x, p = displacement_p)
end

function _bounded_quadrature_displacement(
        data::_QuadratureLinearReference, t::Num,
    )
    from_x = _bounded_quadrature_drive(data.normalized_drive_x, data, t, Val(:x))
    from_p = _bounded_quadrature_drive(data.normalized_drive_p, data, t, Val(:p))
    return (
        x = _add_cnum(from_x.x, from_p.x),
        p = _add_cnum(from_x.p, from_p.p),
    )
end

"""
    Displace(x, p, Hlin)

Construct the static equilibrium displacement generated by a Hermitian one-mode quadratic
quadrature reference Hamiltonian. Symbolic coefficients are treated as time-independent
parameters. A structurally singular driven response is rejected.
"""
function Displace(x::Op, p::Op, Hlin::QAdd)
    _phase_pair(x, p, "`Displace`")
    data = _linear_quadrature_hamiltonian(Hlin, x, p)
    displacement = _static_quadrature_displacement(data)
    return _quadrature_displacement(x, p, displacement.x, displacement.p)
end

"""
    Displace(x, p, Hlin, t)

Construct the bounded moving displacement generated by a Hermitian one-mode quadratic
quadrature reference Hamiltonian with finite harmonic linear drives. The homogeneous
solution is set to zero. Structural resonances are rejected; symbolic determinants are
retained as quotients valid away from their zero surfaces.
"""
function Displace(x::Op, p::Op, Hlin::QAdd, t::Num)
    _phase_pair(x, p, "`Displace`")
    tt = _time_or_throw(t)
    data = _linear_quadrature_hamiltonian(Hlin, x, p)
    _validate_quadrature_time(data, tt)
    displacement = _bounded_quadrature_displacement(data, tt)

    derivative_x = _add_cnum(
        _add_cnum(
            _mul_cnum(data.quadratic_cross, displacement.x),
            _mul_cnum(data.quadratic_p, displacement.p),
        ),
        data.normalized_drive_p,
    )
    derivative_p = _neg_cnum(
        _add_cnum(
            _add_cnum(
                _mul_cnum(data.quadratic_x, displacement.x),
                _mul_cnum(data.quadratic_cross, displacement.p),
            ),
            data.normalized_drive_x,
        ),
    )
    transformed_x_coefficient = _add_cnum(
        _add_cnum(
            _mul_cnum(data.quadratic_x, displacement.x),
            _mul_cnum(data.quadratic_cross, displacement.p),
        ),
        data.drive_x,
    )
    transformed_p_coefficient = _add_cnum(
        _add_cnum(
            _mul_cnum(data.quadratic_cross, displacement.x),
            _mul_cnum(data.quadratic_p, displacement.p),
        ),
        data.drive_p,
    )
    gauge = _quadrature_displacement_gauge(
        x, p, displacement.x, displacement.p, derivative_x, derivative_p,
        _neg_cnum(transformed_x_coefficient),
        _neg_cnum(transformed_p_coefficient),
    )
    static = _quadrature_displacement(x, p, displacement.x, displacement.p)
    return _timed_transform(static, gauge, tt)
end

# === Pauli and spin transformations ===

_axis_or_throw(axis::Integer) = 1 <= axis <= 3 ? Int32(axis) :
    _unitary_error("axis must be 1, 2, or 3, got $axis")
_axis_op(o::Op, axis::Integer) =
    Op(o.kind, o.name_id, o.space_index, o.index, Int32(axis), 0, 0, 0)

function _triple_or_throw(S::Op)
    (is_pauli(S) || is_spin(S)) || _unitary_error(
        "`Rotation` expects a Pauli or Spin operator, got $(S.kind)",
    )
    return nothing
end

"""Rotate a Pauli or spin triple by `θ` around axis 1, 2, or 3."""
function Rotation(S::Op, axis::Integer, θ::Real)
    _triple_or_throw(S)
    main_axis = Int(_axis_or_throw(axis))
    u = _axis_op(S, mod1(main_axis + 1, 3))
    v = _axis_op(S, mod1(main_axis + 2, 3))
    fixed = _axis_op(S, main_axis)
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    negative_s = _neg_cnum(s)
    return _static_transform(
        Dict{Op, QAdd}(
            u => _rule_qadd((c, Op[u]), (negative_s, Op[v])),
            v => _rule_qadd((c, Op[v]), (s, Op[u])),
            fixed => _scaled(_CNUM_ONE, fixed),
        ),
        Dict{Op, QAdd}(
            u => _rule_qadd((c, Op[u]), (s, Op[v])),
            v => _rule_qadd((c, Op[v]), (negative_s, Op[u])),
            fixed => _scaled(_CNUM_ONE, fixed),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function Rotation(S::Op, axis::Integer, θ::Real, t::Num)
    _triple_or_throw(S)
    tt = _time_or_throw(t)
    fixed = _axis_op(S, Int(_axis_or_throw(axis)))
    generator = S.kind === OP_PAULI ? _scaled(_CNUM_HALF, fixed) :
        _scaled(_CNUM_ONE, fixed)
    return _timed_transform(Rotation(S, axis, θ), _gauge(generator, θ, tt), tt)
end

# === Ordinary N-level basis transformations ===

_transition_op(o::Op, i::Integer, j::Integer) =
    Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)

function _nlevel_or_throw(σ::Op, W::AbstractMatrix)
    is_transition(σ) || _unitary_error(
        "`Rotation` expects an ordinary `Transition` operator, got $(σ.kind)",
    )
    n = Int(σ.nlev)
    size(W) == (n, n) || _unitary_error(
        "`W` must be $n×$n for a $n-level space, got $(size(W))",
    )
    return n
end

function _coefficient_matrix(W::AbstractMatrix)
    n, m = size(W)
    converted = Matrix{Coeff}(undef, n, m)
    for j in 1:m, i in 1:n
        converted[i, j] = _to_cnum(W[i, j])
    end
    return converted
end

function _exact_unitary_or_throw(W::Matrix{Coeff}, Wdagger::Matrix{Coeff})
    n = size(W, 1)
    scratch = ParamRelation[]
    for j in 1:n, k in 1:n
        residual = j == k ? _CNUM_NEG1 : _CNUM_ZERO
        for i in 1:n
            residual = _add_cnum(
                residual, _mul_cnum(Wdagger[j, i], W[i, k]),
            )
        end
        reduced = _reduce_all(residual, ParamRelation[], false, scratch)
        _iszero_cnum(reduced) || _unitary_error(
            "`W` must be provably unitary exactly; entry ($j, $k) of `W'W - I` is " *
                "`$(to_num(reduced))`",
        )
    end
    return W
end

function _matrix_unit_rules(
        σ::Op, W::Matrix{Coeff}, Wdagger::Matrix{Coeff} = _dagger_coefficients(W),
    )
    n = size(W, 1)
    operator_terms = QTerm[
        QTerm(Op[_transition_op(σ, k, l)], _EMPTY_NE) for k in 1:n, l in 1:n
    ]
    nonzero_rows = [
        [(column, W[row, column]) for column in 1:n if !_iszero_cnum(W[row, column])]
            for row in 1:n
    ]
    rules = Dict{Op, QAdd}()
    for i in 1:n, j in 1:n
        terms = QTermDict()
        for (k, _) in nonzero_rows[i], (l, wjl) in nonzero_rows[j]
            coefficient = _mul_cnum(Wdagger[k, i], wjl)
            _iszero_cnum(coefficient) ||
                _addto_key!(terms, operator_terms[k, l], coefficient)
        end
        rules[_transition_op(σ, i, j)] = QAdd(terms, _EMPTY_INDICES)
    end
    return rules
end

function _dagger_coefficients(W::Matrix{Coeff})
    n = size(W, 1)
    out = Matrix{Coeff}(undef, n, n)
    for j in 1:n, i in 1:n
        out[i, j] = _conj_cnum(W[j, i])
    end
    return out
end

function _nlevel_rotation(σ::Op, W::AbstractMatrix)
    _nlevel_or_throw(σ, W)
    coefficients = _coefficient_matrix(W)
    dagger = _dagger_coefficients(coefficients)
    _exact_unitary_or_throw(coefficients, dagger)
    U = _static_transform(
        _matrix_unit_rules(σ, coefficients, dagger),
        _matrix_unit_rules(σ, dagger, coefficients),
    )
    return U, coefficients
end

"""
    Rotation(σ, W)
    Rotation(σ, W, t)

Rotate an ordinary N-level basis by a matrix `W` satisfying `W'W = I`. The timed form
derives the Hamiltonian gauge `im*Ẇ'W` entrywise with respect to `t`.
"""
function Rotation(σ::Op, W::AbstractMatrix)
    U, _ = _nlevel_rotation(σ, W)
    return U
end

function _nlevel_gauge(σ::Op, W::Matrix{Coeff}, t::Num)
    n = size(W, 1)
    conjugated_derivatives = Matrix{Coeff}(undef, n, n)
    for j in 1:n, k in 1:n
        conjugated_derivatives[k, j] = _conj_cnum(_dt(W[k, j], t))
    end
    gauge = QTermDict()
    for j in 1:n, l in 1:n
        coefficient = _CNUM_ZERO
        for k in 1:n
            coefficient = _add_cnum(
                coefficient, _mul_cnum(conjugated_derivatives[k, j], W[k, l]),
            )
        end
        coefficient = _mul_cnum(_CNUM_IM, coefficient)
        _iszero_cnum(coefficient) ||
            _addto!(gauge, Op[_transition_op(σ, j, l)], coefficient)
    end
    return QAdd(gauge, _EMPTY_INDICES)
end

function Rotation(σ::Op, W::AbstractMatrix, t::Num)
    tt = _time_or_throw(t)
    U, coefficients = _nlevel_rotation(σ, W)
    return _timed_transform(U, _nlevel_gauge(σ, coefficients, tt), tt)
end
