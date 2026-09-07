using SecondQuantizedAlgebra
using Test
using Symbolics: @variables

@testset "Affine block architecture" begin
    SQA = SecondQuantizedAlgebra

    fock = FockSpace(:fock)
    a = Destroy(fock, :a)
    b = Destroy(fock, :b)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    atom = NLevelSpace(:atom, 2)
    σ = Transition(atom, :σ, 1, 2)

    @variables θ ϕ η dx dp

    @testset "disjoint algebras remain separate blocks" begin
        phase_rotation = Rotation(x, p, θ)
        level_rotation = Rotation(σ, [0 1; 1 0])
        mixed = phase_rotation * level_rotation

        @test mixed.action isa SQA.AffineAction
        @test length(mixed.action.blocks) == 2
        @test Set(block.structure for block in mixed.action.blocks) == Set(
            (
                SQA.SymplecticPhaseSpace(), SQA.UnitaryLinearAction(),
            )
        )

        inverse = inv(mixed)
        @test length(inverse.action.blocks) == 2
        for op in (x, p, σ)
            @test iszero(simplify(conjugate(conjugate(op, mixed), inverse) - op))
        end
    end

    @testset "overlapping transforms merge within one algebra" begin
        composed = Rotation(x, p, θ) * Displace(x, p, dx, dp)
        @test length(composed.action.blocks) == 1
        block = only(composed.action.blocks)
        @test block.structure === SQA.SymplecticPhaseSpace()
        for op in (x, p)
            sequential = conjugate(conjugate(op, Rotation(x, p, θ)), Displace(x, p, dx, dp))
            @test iszero(simplify(conjugate(op, composed) - sequential))
        end
    end

    @testset "subset composition matches sequential application" begin
        first = Rotation(a, b, θ)
        second = Displace(a, η)
        third = Rotation(b, ϕ)
        composed = first * second * third

        @test length(composed.action.blocks) == 1
        @test only(composed.action.blocks).structure === SQA.BosonicNambu()

        inverse = inv(composed)
        for op in (a, b, adjoint(a), adjoint(b))
            sequential = conjugate(conjugate(conjugate(op, first), second), third)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "mixed raw affine bases are not flattened" begin
        @test_throws ArgumentError SQA.AffineAction(
            Op[a, x], [1 0; 0 1], [0, 0],
        )
    end

    @testset "scalar substitution recompiles affine metadata" begin
        U = Rotation(a, θ)
        resolved = @inferred substitute(U, Dict(θ => 0))
        @test resolved isa UnitaryTransform
        @test iszero(simplify(conjugate(a, resolved) - a))
        @test iszero(simplify(conjugate(a, inv(resolved)) - a))
    end
end
