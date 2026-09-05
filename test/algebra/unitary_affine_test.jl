using SecondQuantizedAlgebra
using Test
using Symbolics: @variables

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

    @variables θ r dx dp
    @variables α::Number

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
end
