using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim

@testset "Affine unitary transformation contracts" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    atom = NLevelSpace(:atom, 2)
    σ11 = Transition(atom, :σ, 1, 1)
    σ12 = Transition(atom, :σ, 1, 2)

    @variables θ r ω Ω η g t dx dp ωd K
    @variables α::Number
    @variables envelope(t)

    @testset "explicit named constructors share the compatibility maps" begin
        beam = BeamSplitter(left, right, θ)
        rotation = Rotation(left, right, θ)
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, beam) - conjugate(op, rotation)))
            @test iszero(
                simplify(conjugate(op, inv(beam)) - conjugate(op, inv(rotation))),
            )
        end

        two_mode = TwoModeSqueeze(left, right, r)
        squeeze = Squeeze(left, right, r)
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, two_mode) - conjugate(op, squeeze)))
            @test iszero(
                simplify(conjugate(op, inv(two_mode)) - conjugate(op, inv(squeeze))),
            )
        end

        W = [0 1; 1 0]
        basis_rotation = BasisRotation(σ12, W)
        compatibility = Rotation(σ12, W)
        @test iszero(
            simplify(conjugate(σ11, basis_rotation) - conjugate(σ11, compatibility)),
        )
        @test iszero(
            simplify(conjugate(σ12, basis_rotation) - conjugate(σ12, compatibility)),
        )
    end

    @testset "affine IR validates basis layouts" begin
        SQA = SecondQuantizedAlgebra
        @test_throws ArgumentError SQA.AffineAction(
            SQA.GenericAffine(), Op[a], [1 0; 0 1], [0],
        )
        @test_throws ArgumentError SQA.AffineAction(
            SQA.GenericAffine(), Op[a], reshape([1], 1, 1), [0, 0],
        )
        @test_throws ArgumentError SQA.AffineAction(
            SQA.GenericAffine(), Op[a, a], [1 0; 0 1], [0, 0],
        )
        @test_throws ArgumentError SQA.AffineAction(
            Op[a], reshape([1], 1, 1), [0],
        )
    end

    @testset "legacy compiled-rule transforms remain composable" begin
        SQA = SecondQuantizedAlgebra
        source = Rotation(a, θ)
        legacy = SQA.validated_transform(
            copy(source.rules), copy(source.inverse_rules), gauge_term(source),
            SQA.StaticTime(),
        )
        @test legacy.action === nothing
        @test iszero(simplify(conjugate(a, legacy) - conjugate(a, source)))
        @test iszero(simplify(conjugate(conjugate(a, legacy), inv(legacy)) - a))

        doubled = legacy * legacy
        @test iszero(
            simplify(
                conjugate(a, doubled) - conjugate(conjugate(a, legacy), legacy),
            ),
        )

        affine = Displace(a, α)
        legacy_first = legacy * affine
        affine_first = affine * legacy
        @test iszero(
            simplify(
                conjugate(a, legacy_first) - conjugate(conjugate(a, legacy), affine),
            ),
        )
        @test iszero(
            simplify(
                conjugate(a, affine_first) - conjugate(conjugate(a, affine), legacy),
            ),
        )
    end

    @testset "affine composition preserves family semantics" begin
        phase_composed = Rotation(x, p, θ) * Displace(x, p, dx, dp)
        for op in (x, p)
            sequential = conjugate(conjugate(op, Rotation(x, p, θ)), Displace(x, p, dx, dp))
            @test iszero(simplify(conjugate(op, phase_composed) - sequential))
        end
        phase_inverse = inv(phase_composed)
        @test iszero(simplify(conjugate(conjugate(x + p, phase_composed), phase_inverse) - x - p))

        level_swap = BasisRotation(σ12, [0 1; 1 0])
        mixed = Rotation(a, θ) * level_swap
        mixed_inverse = inv(mixed)
        for op in (a, a', σ11, σ12)
            @test iszero(simplify(conjugate(conjugate(op, mixed), mixed_inverse) - op))
        end
        @test iszero(
            simplify(
                conjugate(a + σ11, mixed) -
                    conjugate(conjugate(a + σ11, Rotation(a, θ)), level_swap),
            ),
        )
    end

    @testset "static Fock displacement frame" begin
        reference = ω * a' * a + η * (a + a') + g
        U = DisplacementFrame(a, reference)

        @test U isa UnitaryTransform
        @test iszero(simplify(conjugate(a, U) - (a - η / ω)))
        @test iszero(
            simplify(transform(reference, U) - (ω * a' * a - η^2 / ω + g)),
        )
        @test iszero(simplify(conjugate(a, inv(U)) - (a + η / ω)))

        undriven = DisplacementFrame(a, ω * a' * a)
        @test iszero(simplify(conjugate(a, undriven) - a))
    end

    @testset "bounded harmonic Fock displacement frame" begin
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')
        U = DisplacementFrame(a, reference, t)

        expected =
            im * Ω / (2 * (ωd - ω)) * expim(-ωd * t) -
            im * Ω / (2 * (ωd + ω)) * expim(ωd * t)
        @test iszero(simplify(conjugate(a, U) - a - expected))

        transformed = transform(reference, U)
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops == Op[a', a]
        end

        constant = DisplacementFrame(a, ω * a' * a + η * (a + a'), t)
        @test iszero(simplify(conjugate(a, constant) - (a - η / ω)))

        numeric_constant = DisplacementFrame(a, ω * a' * a + a + a', t)
        @test iszero(simplify(conjugate(a, numeric_constant) - (a - 1 / ω)))

        multitone_drive = η * cos(ωd * t) + g * sin(2ωd * t)
        multitone_reference = ω * a' * a + multitone_drive * (a + a')
        multitone = DisplacementFrame(a, multitone_reference, t)
        for (term, _) in transform(multitone_reference, multitone)
            isempty(term.ops) && continue
            @test term.ops == Op[a', a]
        end

        Kerr = (K / 2) * a'^2 * a^2
        @test iszero(
            simplify(
                transform(reference + Kerr, U) -
                    transform(reference, U) - conjugate(Kerr, U),
            ),
        )
    end

    @testset "Fock frame refusal is exact" begin
        other = Destroy(FockSpace(:other), :b)
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')

        @test_throws ArgumentError DisplacementFrame(a, reference + K * a'^2 * a^2, t)
        @test_throws ArgumentError DisplacementFrame(a, reference + g * other, t)
        @test_throws ArgumentError DisplacementFrame(a, ω * a' * a + η * a', t)
        @test_throws ArgumentError DisplacementFrame(a, im * a' * a + η * (a + a'))
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + η * (a + a') + im,
        )
        @test_throws ArgumentError DisplacementFrame(
            a, (ω + g * cos(ωd * t)) * a' * a + η * (a + a'), t,
        )
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + envelope * (a + a'), t,
        )

        indexed_fock = FockSpace(:indexed_reference)
        i = Index(indexed_fock, :i, 3, indexed_fock)
        indexed_a = IndexedOperator(Destroy(indexed_fock, :c), i)
        @test_throws ArgumentError DisplacementFrame(
            a, Σ(indexed_a' * indexed_a, i),
        )

        nonlinear_phase = expim(t^2)
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + nonlinear_phase * a' + conj(nonlinear_phase) * a, t,
        )

        resonant = expim(-ω * t)
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + resonant * a' + conj(resonant) * a, t,
        )
        @test_throws ArgumentError DisplacementFrame(a, η * (a + a'), t)
        @test_throws ArgumentError DisplacementFrame(a, η * (a + a'))
    end

    @testset "static quadrature displacement frame" begin
        reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 + η * x + dx * p
        U = DisplacementFrame(x, p, reference)
        determinant = ω * Ω - g^2

        @test U isa UnitaryTransform
        @test iszero(
            simplify(conjugate(x, U) - x - (g * dx - Ω * η) / determinant),
        )
        @test iszero(
            simplify(conjugate(p, U) - p - (g * η - ω * dx) / determinant),
        )

        transformed = simplify(transform(reference, U))
        allowed = Set((Op[x, x], Op[x, p], Op[p, p]))
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops in allowed
        end
    end

    @testset "bounded harmonic quadrature displacement frame" begin
        reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        U = DisplacementFrame(x, p, reference, t)
        @test U isa UnitaryTransform

        transformed = transform(reference, U)
        allowed = Set((Op[x, x], Op[p, p]))
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops in allowed
        end

        numeric_constant = DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + x + p, t,
        )
        @test numeric_constant isa UnitaryTransform

        multitone_reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 +
            (η * cos(ωd * t) + dx * sin(2ωd * t)) * x +
            dp * cos(3ωd * t) * p
        multitone = DisplacementFrame(x, p, multitone_reference, t)
        transformed_multitone = transform(multitone_reference, multitone)
        allowed_multitone = Set((Op[x, x], Op[x, p], Op[p, p]))
        for (term, _) in transformed_multitone
            isempty(term.ops) && continue
            @test term.ops in allowed_multitone
        end
    end

    @testset "quadrature frame refusal is exact" begin
        reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        other_phase = PhaseSpace(:selected) ⊗ PhaseSpace(:other)
        other_x = Position(other_phase, :other_x, 2)
        nonlinear_phase = expim(t^2)

        @test_throws ArgumentError DisplacementFrame(p, x, reference, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + other_x, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + K * x^3, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + im * x, t)
        @test_throws ArgumentError DisplacementFrame(
            x, p,
            ((ω + g * cos(ωd * t)) / 2) * x^2 + (Ω / 2) * p^2 + η * x,
            t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + envelope * x, t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p,
            (ω / 2) * (x^2 + p^2) +
                (nonlinear_phase + conj(nonlinear_phase)) * x,
            t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + η * cos(ω * t) * x, t,
        )
        @test_throws ArgumentError DisplacementFrame(x, p, (ω / 2) * x^2 + η * p)

        nonunique = DisplacementFrame(x, p, (ω / 2) * x^2)
        @test iszero(simplify(conjugate(x + p, nonunique) - (x + p)))
    end
end
