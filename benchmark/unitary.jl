function gaussian_family_workflow(a, b, t, Δa, Δb, J, η, κ, θ, r, ϕ, ω, α)
    H =
        Δa * a' * a +
        Δb * b' * b +
        J * (a' * b + b' * a) +
        η * (a + a') +
        κ * (a' * b' + b * a)

    static_transforms = (
        Displace(a, α),
        Rotation(a, θ),
        Squeeze(a, r, ϕ),
        Rotation(a, b, θ),
        Squeeze(a, b, r),
    )
    timed_transforms = (
        Displace(a, (η + im * ω * t) * t, t),
        Rotation(a, ω * t, t),
        Squeeze(a, r * t, ϕ + ω * t, t),
        Rotation(a, b, θ * t, t),
        Squeeze(a, b, r * t, t),
    )

    static_results = map(U -> transform(H, U), static_transforms)
    timed_results = map(U -> transform(H, U), timed_transforms)
    gauges = map(gauge_term, timed_transforms)
    return static_results, timed_results, gauges
end

function frame_composition_workflow(a, b, t, Δa, Δb, J, η, θ, ω)
    H = Δa * a' * a + Δb * b' * b + J * (a' * b + b' * a) + η * (b + b')

    static = Rotation(a, b, θ) * Displace(a, η) * Rotation(b, θ / 2)
    timed = Rotation(a, b, J * t, t) * Displace(b, η * t, t) * Rotation(a, ω * t, t)

    static_result = transform(H, static)
    timed_result = transform(H, timed)
    round_trip = conjugate(conjugate(a + b, static), inv(static))
    return static_result, timed_result, simplify(round_trip), generators(static), gauge_term(timed)
end

function phase_space_spin_workflow(
        x, p, Sx, Sy, Sz, σx, σy, σz, t, ω, Ω, δ, g, dx, dp, r,
    )
    H =
        (ω / 2) * (x * x + p * p) +
        Ω * Sz +
        (δ / 2) * σz +
        g * x * (Sx + σx) +
        (g / 3) * Sy * σy

    transforms = (
        Displace(x, p, dx * t, dp * t^2, t),
        Rotation(x, p, ω * t, t),
        Squeeze(x, p, r * t, t),
        Rotation(Sx, 1, Ω * t, t),
        Rotation(Sx, 2, Ω * t, t),
        Rotation(Sx, 3, Ω * t, t),
        Rotation(σx, 1, δ * t, t),
        Rotation(σx, 2, δ * t, t),
        Rotation(σx, 3, δ * t, t),
    )
    results = map(U -> transform(H, U), transforms)
    gauges = map(gauge_term, transforms)
    return results, gauges
end

function moving_two_level_workflow(σ, transitions, t, Ω, Δ1, Δ2, g)
    angle = Ω * t
    W = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    σ11, σ22, σ12, σ21 = transitions
    H = Δ1 * σ11 + Δ2 * σ22 + g * (σ12 + σ21)

    static = Rotation(σ, [0 1; 1 0])
    moving = Rotation(σ, W, t)
    return transform(H, static), simplify(transform(H, moving)), gauge_term(moving)
end


function phase_harmonic_workflow(a, t, ω, ν, r)
    argument = ω * t
    series = 0 * a
    merged = SecondQuantizedAlgebra.expim(argument)
    for harmonic in -16:16
        x = (ω + harmonic * ν) * t
        phase = SecondQuantizedAlgebra.expim(x)
        weight = harmonic + 17
        series += weight * (phase * a' + conj(phase) * a)
        harmonic > 0 && (merged *= phase)
    end

    identity = cos(argument)^2 + sin(argument)^2 + cosh(r)^2 - sinh(r)^2
    reduced = simplify(identity * (a + a'))
    differentiated = derivative(merged, t)
    evaluated = substitute(merged, Dict(ω => 2.0, ν => 0.25, t => 0.5))
    pair =
        SecondQuantizedAlgebra.expim(argument) * a' +
        conj(SecondQuantizedAlgebra.expim(argument)) * a
    trigonometric = SecondQuantizedAlgebra.trigonometric_form(pair)
    reconstructed = SecondQuantizedAlgebra.exponential_form(trigonometric)
    rendered = sprint(show, series)
    return series, merged, reduced, differentiated, evaluated, trigonometric, reconstructed,
        rendered
end

function benchmark_unitary!(SUITE)
    # CI compares complete workflows rather than individual phase products or rule
    # lookups. Millisecond-scale leaves collect enough samples to be stable while
    # retaining construction, application, gauges, composition, and rendering in
    # the measured paths.
    @variables θ::Real ϕ::Real r::Real
    @variables ω::Real ν::Real Ω::Real δ::Real t::Real
    @variables Δa::Real Δb::Real Δ1::Real Δ2::Real
    @variables J::Real η::Real κ::Real g::Real dx::Real dp::Real
    @variables α::Number

    fock = FockSpace(:a) ⊗ FockSpace(:b)
    a = Destroy(fock, :a, 1)
    b = Destroy(fock, :b, 2)

    composite = PhaseSpace(:oscillator) ⊗ SpinSpace(:spin) ⊗ PauliSpace(:qubit)
    x = Position(composite, :x, 1)
    p = Momentum(composite, :p, 1)
    Sx = Spin(composite, :S, 1, 2)
    Sy = Spin(composite, :S, 2, 2)
    Sz = Spin(composite, :S, 3, 2)
    σx = Pauli(composite, :σ, 1, 3)
    σy = Pauli(composite, :σ, 2, 3)
    σz = Pauli(composite, :σ, 3, 3)

    atom = NLevelSpace(:atom, 2)
    σ = Transition(atom, :σ, 1, 2)
    transitions = (
        Transition(atom, :σ, 1, 1),
        Transition(atom, :σ, 2, 2),
        Transition(atom, :σ, 1, 2),
        Transition(atom, :σ, 2, 1),
    )

    group = SUITE["Unitary and exact phase workflows"]
    group["Fock Gaussian constructor family"] = @benchmarkable gaussian_family_workflow(
        $a, $b, $t, $Δa, $Δb, $J, $η, $κ, $θ, $r, $ϕ, $ω, $α,
    ) seconds = 3 evals = 1
    group["Static and timed frame composition"] = @benchmarkable frame_composition_workflow(
        $a, $b, $t, $Δa, $Δb, $J, $η, $θ, $ω,
    ) seconds = 3 evals = 1
    group["Phase-space, spin, and Pauli family"] =
        @benchmarkable phase_space_spin_workflow(
        $x, $p, $Sx, $Sy, $Sz, $σx, $σy, $σz, $t, $ω, $Ω, $δ, $g, $dx, $dp, $r,
    ) seconds = 3 evals = 1
    group["Static and moving two-level basis"] = @benchmarkable moving_two_level_workflow(
        $σ, $transitions, $t, $Ω, $Δ1, $Δ2, $g,
    ) seconds = 3 evals = 1
    group["Thirty-three-sideband exact phase pipeline"] = @benchmarkable phase_harmonic_workflow(
        $a, $t, $ω, $ν, $r,
    ) seconds = 3 evals = 1
    return SUITE
end
