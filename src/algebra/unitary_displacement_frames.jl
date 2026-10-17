struct LinearDrift
    basis::Vector{Op}
    primary::Int
    rate::Matrix{CNum}
    drive::Vector{CNum}
    scale::CNum
end

function fock_frame_basis(modes::FockModes)
    lowering_modes = fock_modes(modes, "`DisplacementFrame`")
    return vcat(lowering_modes, adjoint.(lowering_modes)), length(lowering_modes)
end

function hamiltonian_drift(basis::Vector{Op}, primary::Int, H::QAdd, scale::CNum)
    isempty(H.indices) || unitary_error(
        "`DisplacementFrame` does not support summed reference Hamiltonians",
    )
    displaced = Set(basis)
    for (term, _) in H
        isempty(term.ne) || unitary_error(
            "`DisplacementFrame` does not support constrained reference terms",
        )
        for operator in term.ops
            operator in displaced || unitary_error(
                "`DisplacementFrame` reference acts on `$operator` outside the displaced modes",
            )
        end
    end
    adjoint_H = adjoint(H)
    (isequal(H, adjoint_H) || isequal(exponential_form(H), exponential_form(adjoint_H))) || unitary_error(
        "`DisplacementFrame` requires a Hermitian reference Hamiltonian",
    )
    rate, drive = affine_commutator_data(H, basis)
    isone(scale) && return LinearDrift(basis, primary, rate, drive, scale)
    return LinearDrift(basis, primary, scale .* rate, scale .* drive, scale)
end

function frequency_slot(components::Vector{Tuple{CNum, Vector{CNum}}}, frequency::CNum, n::Int)
    for (k, (other, _)) in enumerate(components)
        isequal(other, frequency) && return k
    end
    push!(components, (frequency, fill(CNUM_ZERO, n)))
    return length(components)
end

function harmonic_drives(drive::Vector{CNum}, t::Num)
    variable = SymbolicUtils.unwrap(t)
    components = Tuple{CNum, Vector{CNum}}[]
    for (row, entry) in enumerate(drive), term in phase_terms(entry)
        coefficient_depends_on(term.amplitude, variable) && unitary_error(
            "`DisplacementFrame` supports finite harmonic drives, not the " *
                "time-dependent envelope `$(to_num(term.amplitude))`",
        )
        frequency = expand(Symbolics.derivative(term.phase, t))
        raw_depends_on(SymbolicUtils.unwrap(frequency), variable) && unitary_error(
            "`DisplacementFrame` requires a phase linear in `$t`; got `$(term.phase)`",
        )
        component = iszero(term.phase) ? term.amplitude : term.amplitude * phase_coeff(term.phase)
        k = frequency_slot(components, to_cnum(frequency), length(drive))
        rhs = last(components[k])
        rhs[row] = rhs[row] + component
    end
    return components
end

function cramer_column(matrix::Matrix{CNum}, column::Int, rhs::Vector{CNum})
    replaced = copy(matrix)
    replaced[:, column] = rhs
    return LinearAlgebra.det(replaced)
end

resonance_check(denominator::CNum, frequency::CNum) = iszero(denominator) && unitary_error(
    "`DisplacementFrame` found a resonant drive at frequency `$(to_num(frequency))`",
)

function steady_shift(drift::LinearDrift, components::Vector{Tuple{CNum, Vector{CNum}}})
    n = length(drift.basis)
    m = drift.primary
    shift = fill(CNUM_ZERO, n)
    for block in connected_components(drift.rate)
        first(block) <= m || continue
        for (frequency, drive) in components
            all(row -> iszero(drive[row]), block) && continue
            if length(block) == 1
                row = only(block)
                denominator = drift.scale * frequency - drift.rate[row, row]
                resonance_check(denominator, frequency)
                shift[row] = shift[row] + drive[row] / denominator
                continue
            end
            rhs = drive[block]
            response = CNum[
                (i == j ? drift.scale * frequency : CNUM_ZERO) - drift.rate[i, j]
                    for i in block, j in block
            ]
            denominator = LinearAlgebra.det(response)
            resonance_check(denominator, frequency)
            for (k, row) in enumerate(block)
                row <= m || continue
                shift[row] = shift[row] + cramer_column(response, k, rhs) / denominator
            end
        end
    end
    for row in (m + 1):n
        shift[row] = conj(shift[row - m])
    end
    return shift
end

function drift_velocity(drift::LinearDrift, shift::Vector{CNum}, drive::Vector{CNum})
    factor = CNUM_IM / drift.scale
    velocity = similar(drive)
    for row in eachindex(drive)
        total = drive[row]
        for column in eachindex(shift)
            entry = drift.rate[row, column]
            iszero(entry) || (total = entry * shift[column] + total)
        end
        velocity[row] = factor * total
    end
    return velocity
end

displacement_frame(drift::LinearDrift) =
    shift_transform(drift.basis, steady_shift(drift, [(CNUM_ZERO, drift.drive)]))

function displacement_frame(drift::LinearDrift, t::Num)
    tt = time_or_throw(t)
    variable = SymbolicUtils.unwrap(tt)
    for entry in drift.rate
        coefficient_depends_on(entry, variable) && unitary_error(
            "`DisplacementFrame` requires a time-independent linear response; " *
                "got `$(to_num(entry))`",
        )
    end
    normalized = map(exponential_form, drift.drive)
    shift = steady_shift(drift, harmonic_drives(normalized, tt))
    gauge = displacement_gauge(
        drift.basis, shift,
        drift_velocity(drift, shift, drift.drive), drift_velocity(drift, shift, normalized),
    )
    return timed_transform(shift_transform(drift.basis, shift), gauge, tt)
end

"""
    DisplacementFrame(a, H)
    DisplacementFrame(modes, H)
    DisplacementFrame(x, p, H)
    DisplacementFrame(a, H, t)
    DisplacementFrame(modes, H, t)
    DisplacementFrame(x, p, H, t)
    DisplacementFrame(modes, L)
    DisplacementFrame(modes, L, t)

Displacement to the classical steady state of driven modes.

`H` is a Hermitian Hamiltonian of the displaced modes that is at most quadratic, such as
``\\Delta a^\\dagger a + \\eta a^\\dagger + \\eta^* a``. The transformation displaces each
mode by its steady-state amplitude, here ``-\\eta/\\Delta``, which removes the drive.
`modes` is a mode or a collection of distinct modes, and the quadrature form takes one
`(Position, Momentum)` pair. With a time variable `t` the drive may be a finite sum of
phases [`expim`](@ref)`(ω t)`, and the frame follows the steady driven response, without the
transient free oscillation. A drive on resonance throws an `ArgumentError`, and symbolic
denominators are kept as written.

With a damped linear [`Lindbladian`](@ref) `L` in place of `H`, the amplitudes are the fixed
point of the drift `transpose(L)(a) = λ a + c`, so loss and heating enter them as well.

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
DisplacementFrame(modes::FockModes, H::QAdd) =
    displacement_frame(hamiltonian_drift(fock_frame_basis(modes)..., H, CNUM_ONE))
DisplacementFrame(modes::FockModes, H::QAdd, t::Num) =
    displacement_frame(hamiltonian_drift(fock_frame_basis(modes)..., H, CNUM_ONE), t)

quadrature_frame_basis(x::Op, p::Op) = (phase_pair(x, p, "`DisplacementFrame`"); Op[x, p])

DisplacementFrame(x::Op, p::Op, H::QAdd) =
    displacement_frame(hamiltonian_drift(quadrature_frame_basis(x, p), 2, H, CNUM_IM))
DisplacementFrame(x::Op, p::Op, H::QAdd, t::Num) =
    displacement_frame(hamiltonian_drift(quadrature_frame_basis(x, p), 2, H, CNUM_IM), t)
