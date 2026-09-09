using SecondQuantizedAlgebra
using Test
using Symbolics: @variables

@testset "Affine transformation composition" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)
    b = Destroy(fock, :b)
    c = Destroy(fock, :c)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    atom = NLevelSpace(:atom, 2)
    σ = Transition(atom, :σ, 1, 2)

    @variables θ ϕ η α β γ dx dp t

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

    @testset "partially overlapping mode rotations extend the affine basis" begin
        first = Rotation(a, b, θ)
        second = Rotation(b, c, ϕ)
        composed = first * second
        inverse = inv(composed)

        for op in (a, b, c, adjoint(a), adjoint(b), adjoint(c))
            sequential = conjugate(conjugate(op, first), second)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "one affine block can merge several existing blocks" begin
        displacements = Displace(a, α) * Displace(b, β) * Displace(c, γ)
        identity = [1 0 0; 0 1 0; 0 0 1]
        zeros3 = zeros(Int, 3, 3)
        joint = Bogoliubov((a, b, c), identity, zeros3)
        composed = displacements * joint
        inverse = inv(composed)

        for op in (a, b, c, adjoint(a), adjoint(b), adjoint(c))
            sequential = conjugate(conjugate(op, displacements), joint)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "one block can overlap one member of a multi-block action" begin
        displacements = Displace(a, α) * Displace(b, β)
        rotation = Rotation(a, θ)
        composed = displacements * rotation
        inverse = inv(composed)

        for op in (a, b, adjoint(a), adjoint(b))
            sequential = conjugate(conjugate(op, displacements), rotation)
            @test iszero(simplify(conjugate(op, composed) - sequential))
            @test iszero(simplify(conjugate(conjugate(op, composed), inverse) - op))
        end
    end

    @testset "scalar substitution recompiles transformation semantics" begin
        U = Rotation(a, b, θ)
        resolved = @inferred substitute(U, Dict(θ => 0))
        @test resolved isa UnitaryTransform
        for op in (a, b)
            @test iszero(simplify(conjugate(op, resolved) - op))
            @test iszero(simplify(conjugate(op, inv(resolved)) - op))
        end
    end

    @testset "moving transforms keep their differentiation variable" begin
        U = Rotation(a, θ * t, t)
        @test iszero(simplify(transform(a, U) - conjugate(a, U) - gauge_term(U)))
        @test_throws ArgumentError substitute(U, Dict(t => 0))
    end
end
