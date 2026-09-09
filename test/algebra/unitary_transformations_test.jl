using SecondQuantizedAlgebra
using LinearAlgebra: exp
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis
using Test
using Symbolics: Num, @variables
import SecondQuantizedAlgebra: expim, exponential_form

@testset "Unitary transformation workflows" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)
    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)
    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)
    spin = SpinSpace(:spin)
    Sx = Spin(spin, :S, 1)
    Sy = Spin(spin, :S, 2)
    Sz = Spin(spin, :S, 3)
    pauli = PauliSpace(:pauli)
    σx = Pauli(pauli, :σ, 1)
    σy = Pauli(pauli, :σ, 2)
    σz = Pauli(pauli, :σ, 3)
    atom = NLevelSpace(:atom, 2)
    σ11 = Transition(atom, :σ, 1, 1)
    σ12 = Transition(atom, :σ, 1, 2)
    σ21 = adjoint(σ12)
    σ22 = Transition(atom, :σ, 2, 2)
    @variables α::Number z::Number θ ϕ r ω η g t s dx dp Ω F ω₀
    number = left' * left + right' * right
    weighted_number = Complex(Num(z), Num(z)) * number

    @testset "static closed forms" begin
        U = Displace(a, α)
        @test iszero(simplify(conjugate(a, U) - (a + α)))
        @test iszero(simplify(conjugate(a', inv(U)) - (a' - conj(α))))
        @test iszero(simplify(conjugate(a' * a, U) - (a' + conj(α)) * (a + α)))

        U = Rotation(a, θ)
        @test iszero(simplify(conjugate(a, U) - expim(-θ) * a))
        @test iszero(simplify(conjugate(a', inv(U)) - expim(-θ) * a'))
        @test iszero(simplify(gauge_term(U)))

        spectator = Destroy(FockSpace(:spectator), :spectator)
        @test iszero(simplify(conjugate(spectator, U) - spectator))

        U = Squeeze(a, r, ϕ)
        image = cosh(r) * a + expim(ϕ) * sinh(r) * a'
        @test iszero(simplify(conjugate(a, U) - image))
        @test iszero(
            simplify(
                conjugate(a, inv(U)) -
                    (cosh(r) * a - expim(ϕ) * sinh(r) * a')
            )
        )

        beamsplitter = Rotation(left, right, θ)
        @test iszero(
            simplify(
                conjugate(left, beamsplitter) -
                    (cos(θ) * left + sin(θ) * right)
            )
        )
        @test iszero(
            simplify(
                conjugate(right, beamsplitter) -
                    (cos(θ) * right - sin(θ) * left)
            )
        )
        @test iszero(
            simplify(
                conjugate(left * right, inv(beamsplitter)) -
                    (cos(θ) * left - sin(θ) * right) *
                    (cos(θ) * right + sin(θ) * left)
            )
        )
        @test iszero(simplify(conjugate(number, beamsplitter) - number))
        @test iszero(simplify(conjugate(weighted_number, beamsplitter) - weighted_number))

        two_mode_squeeze = Squeeze(left, right, r)
        @test iszero(
            simplify(
                conjugate(left, two_mode_squeeze) -
                    (cosh(r) * left + sinh(r) * right'),
            )
        )
        @test iszero(
            simplify(
                conjugate(right, two_mode_squeeze) -
                    (cosh(r) * right + sinh(r) * left'),
            )
        )
    end

    @testset "phase-space and spin families" begin
        displacement = Displace(x, p, dx, dp)
        @test iszero(simplify(conjugate(x, displacement) - (x + dx)))
        @test iszero(simplify(conjugate(p, displacement) - (p + dp)))

        quadrature_rotation = Rotation(x, p, θ)
        @test iszero(
            simplify(
                conjugate(x, quadrature_rotation) -
                    (cos(θ) * x + sin(θ) * p)
            )
        )
        @test iszero(
            simplify(
                conjugate(p, quadrature_rotation) -
                    (cos(θ) * p - sin(θ) * x)
            )
        )

        quadrature_squeeze = Squeeze(x, p, r)
        @test iszero(simplify(conjugate(x, quadrature_squeeze) - exp(r) * x))
        @test iszero(simplify(conjugate(p, quadrature_squeeze) - exp(-r) * p))

        spin_rotation = Rotation(Sx, 3, θ)
        @test iszero(
            simplify(
                conjugate(Sx, spin_rotation) -
                    (cos(θ) * Sx - sin(θ) * Sy)
            )
        )
        @test iszero(
            simplify(
                conjugate(Sy, spin_rotation) -
                    (cos(θ) * Sy + sin(θ) * Sx)
            )
        )
        @test iszero(simplify(conjugate(Sz, spin_rotation) - Sz))

        pauli_rotation = Rotation(σx, 2, θ)
        @test iszero(
            simplify(
                conjugate(σz, pauli_rotation) -
                    (cos(θ) * σz - sin(θ) * σx)
            )
        )
        @test iszero(
            simplify(
                conjugate(σx, pauli_rotation) -
                    (cos(θ) * σx + sin(θ) * σz)
            )
        )
        @test iszero(simplify(conjugate(σy, pauli_rotation) - σy))

        for axis in 1:3
            spin_axis_rotation = Rotation(Sx, axis, θ)
            pauli_axis_rotation = Rotation(σx, axis, θ)
            for observable in (Sx, Sy, Sz)
                @test iszero(
                    simplify(
                        conjugate(
                            conjugate(observable, spin_axis_rotation),
                            inv(spin_axis_rotation),
                        ) - observable,
                    ),
                )
            end
            for observable in (σx, σy, σz)
                @test iszero(
                    simplify(
                        conjugate(
                            conjugate(observable, pauli_axis_rotation),
                            inv(pauli_axis_rotation),
                        ) - observable,
                    ),
                )
            end
        end
    end

    @testset "N-level matrix rotations" begin
        W = [0 1; 1 0]
        U = Rotation(σ12, W)
        @test iszero(simplify(conjugate(σ11, U) - σ22))
        @test iszero(simplify(conjugate(σ12, U) - σ21))
        @test iszero(simplify(conjugate(σ11 + 2σ22, U) - (σ22 + 2σ11)))

        phases = [expim(-ω * t) 0; 0 expim(ω * t)]
        timed = Rotation(σ12, phases, t)
        @test iszero(simplify(conjugate(σ12, timed) - expim(2ω * t) * σ12))
        @test iszero(simplify(gauge_term(timed) - (-ω * σ11 + ω * σ22)))
        @test iszero(
            simplify(
                transform(Ω * (σ11 + σ22), timed) -
                    ((Ω - ω) * σ11 + (Ω + ω) * σ22),
            )
        )

        # As with raw Bogoliubov maps, mathematical unitarity is a caller precondition.
        assumed = @inferred Rotation(σ12, [1 0; 0 2])
        @test iszero(simplify(conjugate(σ12, assumed) - 2σ12))
    end

    @testset "numeric oracles validate symbolic transformations" begin
        # Convert the public result and compare it with an independently
        # assembled finite-dimensional matrix transformation. This catches a
        # pair of self-consistent but incorrect symbolic rules.
        fock_dimension = 24
        lowering = zeros(ComplexF64, fock_dimension, fock_dimension)
        for n in 2:fock_dimension
            lowering[n - 1, n] = sqrt(n - 1)
        end
        amplitude = 0.04 + 0.03im
        displacement = exp(amplitude * lowering' - conj(amplitude) * lowering)
        fock_basis = FockBasis(fock_dimension)
        transformed_fock = Matrix(
            to_numeric(conjugate(a, Displace(a, amplitude)), fock_basis).data,
        )
        expected_fock = displacement' * lowering * displacement
        @test transformed_fock[1:6, 1:6] ≈ expected_fock[1:6, 1:6]

        spin_basis = SpinBasis(1 // 2)
        θ₀ = 0.37
        spin_matrix(op) = Matrix(to_numeric(op, spin_basis).data)
        sy = spin_matrix(σy)
        spin_unitary = exp(-0.5im * θ₀ * sy)
        expected_spin = spin_unitary' * spin_matrix(σz) * spin_unitary
        actual_spin = spin_matrix(conjugate(σz, Rotation(σx, 2, θ₀)))
        @test actual_spin ≈ expected_spin

        atom_basis = NLevelBasis(2)
        W = [3 // 5 -4 // 5; 4 // 5 3 // 5]
        level_matrix(op) = Matrix(to_numeric(op, atom_basis).data)
        expected_level = ComplexF64.(W)' * level_matrix(σ11) * ComplexF64.(W)
        actual_level = level_matrix(conjugate(σ11, Rotation(σ12, W)))
        @test actual_level ≈ expected_level

        frequency = 1.7
        instant = 0.31
        step = 1.0e-6
        Wfun(τ) = ComplexF64[
            exp(-im * frequency * τ) 0
            0 exp(im * frequency * τ)
        ]
        derivative = (Wfun(instant + step) - Wfun(instant - step)) / (2step)
        expected_gauge = im * derivative' * Wfun(instant)
        phases = [expim(-ω * t) 0; 0 expim(ω * t)]
        timed_level = Rotation(σ12, phases, t)
        actual_gauge = Matrix(
            to_numeric(
                gauge_term(timed_level), atom_basis;
                parameter = Dict(ω => frequency),
            ).data,
        )
        @test actual_gauge ≈ expected_gauge atol = 1.0e-9
    end

    @testset "timed gauges" begin
        moving_rotation = Rotation(a, ω * t, t)
        @test iszero(simplify(conjugate(a, moving_rotation) - expim(-ω * t) * a))
        @test iszero(simplify(gauge_term(moving_rotation) + ω * a' * a))
        @test iszero(
            simplify(
                transform(Ω * a' * a, moving_rotation) - (Ω - ω) * a' * a,
            )
        )

        moving_squeeze = Squeeze(a, r * t, 0, t)
        @test iszero(
            simplify(
                gauge_term(moving_squeeze) - (im * r / 2) * (a * a - a' * a'),
            )
        )

        moving_beamsplitter = Rotation(left, right, θ * t, t)
        @test iszero(simplify(conjugate(number, moving_beamsplitter) - number))
        @test iszero(
            simplify(conjugate(weighted_number, moving_beamsplitter) - weighted_number),
        )
        @test iszero(
            simplify(conjugate(number^2, moving_beamsplitter) - number^2),
        )
        @test iszero(
            simplify(
                gauge_term(moving_beamsplitter) + im * θ * (left' * right - right' * left),
            )
        )

        moving_two_mode_squeeze = Squeeze(left, right, r * t, t)
        @test iszero(
            simplify(
                gauge_term(moving_two_mode_squeeze) + im * r * (left' * right' - right * left),
            )
        )

        moving_quadratures = Displace(x, p, η * t, ω * t^2, t)
        @test iszero(
            simplify(
                gauge_term(moving_quadratures) -
                    (2ω * t * x - η * p + η * ω * t^2 / 2),
            )
        )
        @test iszero(simplify(gauge_term(Rotation(Sx, 3, θ * t, t)) + θ * Sz))
        @test iszero(
            simplify(
                gauge_term(Rotation(σx, 3, θ * t, t)) + θ * σz / 2,
            )
        )
    end

    @testset "composition and algebraic laws" begin
        U = Rotation(a, θ) * Squeeze(a, r, ϕ)
        observable = a + 2a'
        @test iszero(
            simplify(
                conjugate(a * a', U) - conjugate(a, U) * conjugate(a', U),
            )
        )
        @test iszero(simplify(conjugate(conjugate(observable, U), inv(U)) - observable))
        @test iszero(
            simplify(
                conjugate(observable, U * inv(U)) - observable,
            )
        )

        independent = Destroy(FockSpace(:independent), :b)
        disjoint = Rotation(a, θ) * Rotation(independent, ϕ)
        @test length(generators(disjoint)) == 4
        @test iszero(
            simplify(
                conjugate(a + independent, disjoint) -
                    (expim(-θ) * a + expim(-ϕ) * independent),
            )
        )

        timed = Rotation(a, ω * t, t) * Displace(a, η * t, t)
        H = Ω * a' * a + g * (a + a')
        @test iszero(
            simplify(
                transform(H, timed) - transform(
                    transform(H, Rotation(a, ω * t, t)),
                    Displace(a, η * t, t)
                ),
            )
        )
        @test iszero(
            simplify(
                transform(Ω * a' * a, U) -
                    conjugate(Ω * a' * a, U)
            )
        )
        @test iszero(simplify(transform(H, timed) - adjoint(transform(H, timed))))

        listed = generators(U)
        empty!(listed)
        @test !isempty(generators(U))

        diagonal = Rotation(a, θ) * Rotation(a, ϕ)
        @test iszero(simplify(conjugate(a, diagonal) - expim(-(θ + ϕ)) * a))
        @test iszero(
            simplify(
                conjugate(a, adjoint(diagonal)) - conjugate(a, inv(diagonal)),
            ),
        )
        timed_diagonal = Rotation(a, θ * t, t) * Rotation(a, ϕ * t, t)
        @test iszero(
            simplify(
                transform(a' * a, timed_diagonal) -
                    transform(a' * a, Rotation(a, (θ + ϕ) * t, t)),
            ),
        )

        numeric_diagonal = Rotation(a, 1) * Rotation(a, -1)
        @test iszero(simplify(conjugate(a, numeric_diagonal) - a))
        static_dynamic = Rotation(a, 1) * Rotation(a, ω * t, t)
        @test iszero(
            simplify(
                conjugate(a, static_dynamic) - expim(-1) * expim(-ω * t) * a,
            ),
        )
        static_symbolic_dynamic = Rotation(a, θ) * Rotation(a, ω * t, t)
        @test iszero(
            simplify(
                conjugate(a, static_symbolic_dynamic) -
                    expim(-θ) * expim(-ω * t) * a,
            ),
        )
        dynamic_static = Rotation(a, ω * t, t) * Rotation(a, 1)
        @test iszero(
            simplify(
                conjugate(a, dynamic_static) - expim(-ω * t) * expim(-1) * a,
            ),
        )

        moving_quadrature_rotation = Rotation(x, p, θ * t, t)
        @test iszero(
            simplify(
                gauge_term(moving_quadrature_rotation) + θ * (x * x + p * p) / 2,
            ),
        )
        moving_quadrature_squeeze = Squeeze(x, p, r * t, t)
        @test iszero(
            simplify(
                gauge_term(moving_quadrature_squeeze) - (im * r / 2 - r * x * p),
            ),
        )
    end

    @testset "invalid transformations fail at the public boundary" begin
        @test_throws ArgumentError Displace(σx, 1)
        @test_throws ArgumentError Rotation(x, θ)
        @test_throws ArgumentError Squeeze(Sx, r)
        @test_throws ArgumentError Rotation(x, x, θ)
        @test_throws ArgumentError Rotation(left, left, θ)
        @test_throws ArgumentError Rotation(Sx, 0, θ)
        @test_throws ArgumentError Rotation(Sx, 4, θ)
        @test_throws ArgumentError Rotation(σ12, [1 0 0; 0 1 0])
        @test_throws MethodError Rotation(a, 1 + im)
        @test_throws MethodError Squeeze(a, 1 + im)
        @test_throws MethodError Displace(x, p, 1 + im, 0)
        @test_throws MethodError Rotation(a, ω * t, 1.0)
        @test_throws ArgumentError Rotation(a, ω * t, t) * Rotation(a, ω * s, s)
        @test_throws ArgumentError Rotation(a, ω * t) * Rotation(a, ω * t, t)
    end

    @testset "public transformation boundaries are inferable" begin
        @inferred Displace(a, α)
        @inferred Rotation(a, θ)
        @inferred Squeeze(a, r, ϕ)
        @inferred Rotation(left, right, θ)
        @inferred conjugate(a' * a, Rotation(a, θ))
        @inferred transform(a' * a, Rotation(a, ω * t, t))
        @inferred inv(Rotation(a, ω * t, t))
        @test !occursin(
            "cos(", string(
                exponential_form(
                    transform(
                        a + a', Rotation(a, ω * t, t),
                    )
                )
            )
        )
    end
end
