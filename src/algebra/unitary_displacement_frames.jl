struct FockLinearReference
    frequency::CNum
    drive::CNum
    drive_adjoint::CNum
    normalized_drive::CNum
end

function linear_fock_hamiltonian(H::QAdd, d::Op)
    isempty(H.indices) || unitary_error(
        "`DisplacementFrame` does not support summed reference Hamiltonians",
    )
    frequency = CNUM_ZERO
    drive = CNUM_ZERO
    drive_adjoint = CNUM_ZERO
    scalar = CNUM_ZERO
    raising = adjoint(d)
    for (term, coefficient) in H
        isempty(term.ne) || unitary_error(
            "`DisplacementFrame` does not support constrained reference terms",
        )
        if isempty(term.ops)
            scalar = add_cnum(scalar, coefficient)
        elseif term.ops == Op[raising, d]
            frequency = add_cnum(frequency, coefficient)
        elseif term.ops == Op[raising]
            drive = add_cnum(drive, coefficient)
        elseif term.ops == Op[d]
            drive_adjoint = add_cnum(drive_adjoint, coefficient)
        else
            unsupported = rule_qadd((coefficient, copy(term.ops)))
            unitary_error(
                "`DisplacementFrame` expects only scalar, `$raising * $d`, `$raising`, " *
                    "and `$d` terms; got `$(sprint(show, unsupported))`",
            )
        end
    end

    isequal(frequency, conj_cnum(frequency)) || unitary_error(
        "`DisplacementFrame` requires a real oscillator frequency; " *
            "got `$(to_num(frequency))`",
    )
    isequal(scalar, conj_cnum(scalar)) || unitary_error(
        "`DisplacementFrame` requires a Hermitian scalar reference term",
    )
    normalized_drive = exponential_form(drive)
    normalized_adjoint = exponential_form(drive_adjoint)
    isequal(normalized_adjoint, conj_cnum(normalized_drive)) || unitary_error(
        "`DisplacementFrame` requires adjoint coefficients for `$raising` and `$d`",
    )
    return FockLinearReference(frequency, drive, drive_adjoint, normalized_drive)
end

function validate_fock_frame_time(data::FockLinearReference, t::Num)
    variable = SymbolicUtils.unwrap(t)
    coefficient_depends_on(data.frequency, variable) && unitary_error(
        "`DisplacementFrame` requires a time-independent oscillator frequency; " *
            "got `$(to_num(data.frequency))`",
    )
    return data
end

function harmonic_frequency(monomial::Monomial, t::Num)
    phase_index = 0
    variable = SymbolicUtils.unwrap(t)
    @inbounds for i in eachindex(monomial.syms)
        factor = monomial.syms[i]
        if is_phase(factor)
            phase_index == 0 || unitary_error(
                "`DisplacementFrame` found more than one phase in one drive component",
            )
            phase_index = i
        elseif raw_depends_on(factor, variable)
            unitary_error(
                "`DisplacementFrame` supports finite harmonic drives, not the " *
                    "time-dependent envelope `$(Num(factor))`",
            )
        end
    end
    phase_index == 0 && return CNUM_ZERO

    exponent = monomial.exps[phase_index]
    denominator(exponent) == 1 || unitary_error(
        "`DisplacementFrame` requires integer phase powers; got exponent `$exponent`",
    )
    phase_factor = monomial.syms[phase_index]
    argument = Num(only(SymbolicUtils.arguments(phase_factor))) * numerator(exponent)
    frequency = Symbolics.derivative(argument, t)
    raw_depends_on(SymbolicUtils.unwrap(frequency), variable) && unitary_error(
        "`DisplacementFrame` requires a phase linear in `$t`; " *
            "got `$(Num(only(SymbolicUtils.arguments(phase_factor))))`",
    )
    return to_cnum(frequency)
end

function bounded_harmonic_displacement(drive::CNum, frequency::CNum, t::Num)
    iszero_cnum(drive) && return CNUM_ZERO
    tail = drive.tail
    if tail isa Native
        iszero_cnum(frequency) && unitary_error(
            "`DisplacementFrame` found a resonant constant drive with zero frequency",
        )
        return neg_cnum(drive) / frequency
    elseif !(tail isa Poly)
        unitary_error(
            "`DisplacementFrame` supports only finite harmonic drive coefficients; " *
                "got `$(to_num(drive))`",
        )
    end

    amplitude = CNUM_ZERO
    for monomial in tail.terms
        harmonic = harmonic_frequency(monomial, t)
        divisor = add_cnum(frequency, harmonic)
        component = from_poly([monomial])
        iszero_cnum(divisor) && unitary_error(
            "`DisplacementFrame` found a resonant drive component `$(to_num(component))`",
        )
        amplitude = add_cnum(amplitude, neg_cnum(component) / divisor)
    end
    return amplitude
end

function static_fock_displacement(data::FockLinearReference)
    iszero_cnum(data.drive) && return CNUM_ZERO
    iszero_cnum(data.frequency) && unitary_error(
        "`DisplacementFrame` found a resonant constant drive with zero frequency",
    )
    return neg_cnum(data.drive) / data.frequency
end

"""
    DisplacementFrame(a, H)
    DisplacementFrame(x, p, H)
    DisplacementFrame(a, H, t)
    DisplacementFrame(x, p, H, t)
    DisplacementFrame(a, L)
    DisplacementFrame(a, L, t)

Displacement to the classical steady state of a driven mode.

`H` is a Hamiltonian of one mode that is at most linear in the drive, such as
``\\Delta a^\\dagger a + \\eta a^\\dagger + \\eta^* a``. The transformation displaces the
mode by its steady-state amplitude, here ``-\\eta/\\Delta``, which removes the drive. The
quadrature form accepts any Hermitian quadratic `H`. With a time variable `t` the drive may
be a finite sum of phases [`expim`](@ref)`(ω t)`, and the frame follows the steady driven
response, without the transient free oscillation. A drive on resonance throws an
`ArgumentError`, and symbolic denominators are kept as written.

With a damped linear [`Lindbladian`](@ref) `L` in place of `H`, the amplitude is the fixed
point of the drift `transpose(L)(a) = λ a + c`, so loss and heating enter it as well.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables Δ η::Number;

julia> H = Δ * a' * a + η * a' + conj(η) * a;

julia> conjugate(H, DisplacementFrame(a, H))
(-conj(η)*η) / Δ + Δ * a' * a
```

See also [`Displace`](@ref), [`transform`](@ref).
"""
function DisplacementFrame(a::Op, Hlin::QAdd)
    d = fock_or_throw(a, "`DisplacementFrame`")
    data = linear_fock_hamiltonian(Hlin, d)
    return fock_displacement(d, static_fock_displacement(data))
end

function DisplacementFrame(a::Op, Hlin::QAdd, t::Num)
    d = fock_or_throw(a, "`DisplacementFrame`")
    tt = time_or_throw(t)
    data = linear_fock_hamiltonian(Hlin, d)
    validate_fock_frame_time(data, tt)
    amplitude = bounded_harmonic_displacement(data.normalized_drive, data.frequency, tt)

    conjugate_amplitude = conj_cnum(amplitude)
    derivative = mul_cnum(
        CNUM_NEG_IM,
        add_cnum(mul_cnum(data.frequency, amplitude), data.normalized_drive),
    )
    derivative_adjoint = mul_cnum(
        CNUM_IM,
        add_cnum(
            mul_cnum(data.frequency, conjugate_amplitude),
            conj_cnum(data.normalized_drive),
        ),
    )
    gauge = fock_displacement_gauge(
        d, amplitude, derivative, derivative_adjoint,
        neg_cnum(add_cnum(data.drive, mul_cnum(data.frequency, amplitude))),
        neg_cnum(
            add_cnum(
                data.drive_adjoint,
                mul_cnum(data.frequency, conjugate_amplitude),
            ),
        ),
    )
    return timed_transform(fock_displacement(d, amplitude), gauge, tt)
end

struct QuadratureLinearReference
    quadratic_x::CNum
    quadratic_cross::CNum
    quadratic_p::CNum
    drive_x::CNum
    drive_p::CNum
    normalized_drive_x::CNum
    normalized_drive_p::CNum
end

@enum QuadratureDriveAxis::UInt8 begin
    QUADRATURE_DRIVE_X
    QUADRATURE_DRIVE_P
end

@noinline function quadrature_reference_error(coefficient::CNum, name::AbstractString)
    unitary_error(
        "`DisplacementFrame` requires a real $name coefficient; " *
            "got `$(to_num(coefficient))`",
    )
end

function linear_quadrature_hamiltonian(H::QAdd, x::Op, p::Op)
    isempty(H.indices) || unitary_error(
        "`DisplacementFrame` does not support summed reference Hamiltonians",
    )
    half_quadratic_x = CNUM_ZERO
    quadratic_cross = CNUM_ZERO
    half_quadratic_p = CNUM_ZERO
    drive_x = CNUM_ZERO
    drive_p = CNUM_ZERO
    for (term, coefficient) in H
        isempty(term.ne) || unitary_error(
            "`DisplacementFrame` does not support constrained reference terms",
        )
        if isempty(term.ops)
            continue
        elseif term.ops == Op[x, x]
            half_quadratic_x = add_cnum(half_quadratic_x, coefficient)
        elseif term.ops == Op[x, p]
            quadratic_cross = add_cnum(quadratic_cross, coefficient)
        elseif term.ops == Op[p, p]
            half_quadratic_p = add_cnum(half_quadratic_p, coefficient)
        elseif term.ops == Op[x]
            drive_x = add_cnum(drive_x, coefficient)
        elseif term.ops == Op[p]
            drive_p = add_cnum(drive_p, coefficient)
        else
            unsupported = rule_qadd((coefficient, copy(term.ops)))
            unitary_error(
                "`DisplacementFrame` expects only scalar, `$x * $x`, `$x * $p`, " *
                    "`$p * $p`, `$x`, and `$p` terms; got " *
                    "`$(sprint(show, unsupported))`",
            )
        end
    end

    isequal(H, H') || unitary_error(
        "`DisplacementFrame` requires a Hermitian quadrature reference Hamiltonian",
    )
    quadratic_x = add_cnum(half_quadratic_x, half_quadratic_x)
    quadratic_p = add_cnum(half_quadratic_p, half_quadratic_p)
    for (coefficient, name) in (
            (quadratic_x, "x²"),
            (quadratic_cross, "symmetrized x-p"),
            (quadratic_p, "p²"),
        )
        isequal(coefficient, conj_cnum(coefficient)) ||
            quadrature_reference_error(coefficient, name)
    end
    normalized_drive_x = exponential_form(drive_x)
    normalized_drive_p = exponential_form(drive_p)
    isequal(normalized_drive_x, conj_cnum(normalized_drive_x)) ||
        quadrature_reference_error(normalized_drive_x, "x-drive")
    isequal(normalized_drive_p, conj_cnum(normalized_drive_p)) ||
        quadrature_reference_error(normalized_drive_p, "p-drive")
    return QuadratureLinearReference(
        quadratic_x, quadratic_cross, quadratic_p,
        drive_x, drive_p, normalized_drive_x, normalized_drive_p,
    )
end

function validate_quadrature_frame_time(data::QuadratureLinearReference, t::Num)
    variable = SymbolicUtils.unwrap(t)
    for (coefficient, name) in (
            (data.quadratic_x, "x²"),
            (data.quadratic_cross, "symmetrized x-p"),
            (data.quadratic_p, "p²"),
        )
        coefficient_depends_on(coefficient, variable) && unitary_error(
            "`DisplacementFrame` requires a time-independent $name coefficient; " *
                "got `$(to_num(coefficient))`",
        )
    end
    return data
end

@noinline function quadrature_resonance_error(
        drive_x::CNum, drive_p::CNum, harmonic::CNum,
    )
    unitary_error(
        "`DisplacementFrame` found a resonant quadrature drive component " *
            "`($(to_num(drive_x)), $(to_num(drive_p)))` at frequency " *
            "`$(to_num(harmonic))`",
    )
end

function quadrature_response_component(
        data::QuadratureLinearReference, drive_x::CNum, drive_p::CNum,
        harmonic::CNum,
    )
    imaginary_harmonic = mul_cnum(CNUM_IM, harmonic)
    upper_right = add_cnum(data.quadratic_cross, imaginary_harmonic)
    lower_left = add_cnum(data.quadratic_cross, neg_cnum(imaginary_harmonic))
    determinant = add_cnum(
        mul_cnum(data.quadratic_x, data.quadratic_p),
        neg_cnum(mul_cnum(upper_right, lower_left)),
    )
    iszero_cnum(determinant) &&
        quadrature_resonance_error(drive_x, drive_p, harmonic)

    displacement_x = neg_cnum(
        add_cnum(
            mul_cnum(data.quadratic_p, drive_x),
            neg_cnum(mul_cnum(upper_right, drive_p)),
        ),
    ) / determinant
    displacement_p = add_cnum(
        mul_cnum(lower_left, drive_x),
        neg_cnum(mul_cnum(data.quadratic_x, drive_p)),
    ) / determinant
    return (x = displacement_x, p = displacement_p)
end

function static_quadrature_displacement(data::QuadratureLinearReference)
    (iszero_cnum(data.drive_x) && iszero_cnum(data.drive_p)) &&
        return (x = CNUM_ZERO, p = CNUM_ZERO)
    return quadrature_response_component(data, data.drive_x, data.drive_p, CNUM_ZERO)
end

function bounded_quadrature_drive(
        drive::CNum, data::QuadratureLinearReference, t::Num, axis::QuadratureDriveAxis,
    )
    iszero_cnum(drive) && return (x = CNUM_ZERO, p = CNUM_ZERO)
    tail = drive.tail
    displacement_x = CNUM_ZERO
    displacement_p = CNUM_ZERO
    if tail isa Native
        return axis === QUADRATURE_DRIVE_X ?
            quadrature_response_component(data, drive, CNUM_ZERO, CNUM_ZERO) :
            quadrature_response_component(data, CNUM_ZERO, drive, CNUM_ZERO)
    elseif !(tail isa Poly)
        unitary_error(
            "`DisplacementFrame` supports only finite harmonic drive coefficients; " *
                "got `$(to_num(drive))`",
        )
    end

    for monomial in tail.terms
        harmonic = harmonic_frequency(monomial, t)
        component = from_poly([monomial])
        response = axis === QUADRATURE_DRIVE_X ?
            quadrature_response_component(data, component, CNUM_ZERO, harmonic) :
            quadrature_response_component(data, CNUM_ZERO, component, harmonic)
        displacement_x = add_cnum(displacement_x, response.x)
        displacement_p = add_cnum(displacement_p, response.p)
    end
    return (x = displacement_x, p = displacement_p)
end

function bounded_quadrature_displacement(data::QuadratureLinearReference, t::Num)
    from_x = bounded_quadrature_drive(
        data.normalized_drive_x, data, t, QUADRATURE_DRIVE_X,
    )
    from_p = bounded_quadrature_drive(
        data.normalized_drive_p, data, t, QUADRATURE_DRIVE_P,
    )
    return (
        x = add_cnum(from_x.x, from_p.x),
        p = add_cnum(from_x.p, from_p.p),
    )
end

function DisplacementFrame(x::Op, p::Op, Hlin::QAdd)
    phase_pair(x, p, "`DisplacementFrame`")
    data = linear_quadrature_hamiltonian(Hlin, x, p)
    displacement = static_quadrature_displacement(data)
    return quadrature_displacement(x, p, displacement.x, displacement.p)
end

function DisplacementFrame(x::Op, p::Op, Hlin::QAdd, t::Num)
    phase_pair(x, p, "`DisplacementFrame`")
    tt = time_or_throw(t)
    data = linear_quadrature_hamiltonian(Hlin, x, p)
    validate_quadrature_frame_time(data, tt)
    displacement = bounded_quadrature_displacement(data, tt)

    derivative_x = add_cnum(
        add_cnum(
            mul_cnum(data.quadratic_cross, displacement.x),
            mul_cnum(data.quadratic_p, displacement.p),
        ),
        data.normalized_drive_p,
    )
    derivative_p = neg_cnum(
        add_cnum(
            add_cnum(
                mul_cnum(data.quadratic_x, displacement.x),
                mul_cnum(data.quadratic_cross, displacement.p),
            ),
            data.normalized_drive_x,
        ),
    )
    transformed_x_coefficient = add_cnum(
        add_cnum(
            mul_cnum(data.quadratic_x, displacement.x),
            mul_cnum(data.quadratic_cross, displacement.p),
        ),
        data.drive_x,
    )
    transformed_p_coefficient = add_cnum(
        add_cnum(
            mul_cnum(data.quadratic_cross, displacement.x),
            mul_cnum(data.quadratic_p, displacement.p),
        ),
        data.drive_p,
    )
    gauge = quadrature_displacement_gauge(
        x, p, displacement.x, displacement.p, derivative_x, derivative_p,
        neg_cnum(transformed_x_coefficient),
        neg_cnum(transformed_p_coefficient),
    )
    static = quadrature_displacement(x, p, displacement.x, displacement.p)
    return timed_transform(static, gauge, tt)
end
