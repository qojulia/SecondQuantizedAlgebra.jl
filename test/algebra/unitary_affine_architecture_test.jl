using SecondQuantizedAlgebra
using Test
using Symbolics: @variables

@testset "Affine transformation composition" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)
    b = Destroy(fock, :b)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    atom = NLevelSpace(:atom, 2)
    σ = Transition(atom, :σ, 1, 2)

    @variables θ ϕ η dx dp

    @testset "disjoint algebra transformations compose independently" begin
        phase_rotation = Rotation(x, p, θ)
        level_rotation = Rotation(σ, [0 1; 1 0])
        mixed = phase_rotation * level_rotation
        inverse = inv(mixed)

        for op in (x, p, σ)
            sequential = conjugate(conjugate(op, phase_rotation), level_rotation)
            @test iszero(simplify(conjugate(op, mixed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, mixed), inverse) - op))
        end
    end

    @testset "overlapping transformations agree with sequential application" begin
        rotation = Rotation(x, p, θ)
        displacement = Displace(x, p, dx, dp)
        composed = rotation * displacement
        inverse = inv(composed)

        for op in (x, p)
            sequential = conjugate(conjugate(op, rotation), displacement)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "subset composition matches sequential application" begin
        first = Rotation(a, b, θ)
        second = Displace(a, η)
        third = Rotation(b, ϕ)
        composed = first * second * third
        inverse = inv(composed)

        for op in (a, b, adjoint(a), adjoint(b))
            sequential = conjugate(conjugate(conjugate(op, first), second), third)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "scalar substitution recompiles transformation semantics" begin
        U = Rotation(a, θ)
        resolved = @inferred substitute(U, Dict(θ => 0))
        @test resolved isa UnitaryTransform
        @test iszero(simplify(conjugate(a, resolved) - a))
        @test iszero(simplify(conjugate(a, inv(resolved)) - a))
    end
end
