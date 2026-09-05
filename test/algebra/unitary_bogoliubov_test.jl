using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim

@testset "Exact bosonic Bogoliubov transformations" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)

    @variables r ϕ θ

    @testset "single-mode squeezing is a Bogoliubov map" begin
        S = [
            cosh(r) expim(ϕ) * sinh(r)
            expim(-ϕ) * sinh(r) cosh(r)
        ]
        raw = Bogoliubov(a, S)
        named = Squeeze(a, r, ϕ)

        for op in (a, a')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
            @test iszero(
                simplify(conjugate(op, inv(raw)) - conjugate(op, inv(named))),
            )
        end
        @test iszero(simplify(conjugate(a, raw) - cosh(r) * a - expim(ϕ) * sinh(r) * a'))
    end

    @testset "passive two-mode mixing uses the same Nambu core" begin
        U = [cos(θ) sin(θ); -sin(θ) cos(θ)]
        V = [0 0; 0 0]
        raw = Bogoliubov((left, right), U, V)
        named = BeamSplitter(left, right, θ)

        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
    end

    @testset "two-mode squeezing uses the same Nambu core" begin
        U = [cosh(r) 0; 0 cosh(r)]
        V = [0 sinh(r); sinh(r) 0]
        raw = Bogoliubov(Op[left, right], U, V)
        named = TwoModeSqueeze(left, right, r)

        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
        @test iszero(
            simplify(
                conjugate(left, raw) - cosh(r) * left - sinh(r) * right',
            ),
        )
    end

    @testset "exact canonicality is enforced" begin
        @test_throws ArgumentError Bogoliubov(Op[], Matrix{Int}(undef, 0, 0))
        @test_throws ArgumentError Bogoliubov(a, [1 1; 1 1])
        @test_throws ArgumentError Bogoliubov(a, [1 0; 0 2])
        @test_throws ArgumentError Bogoliubov(a, [1 1; 0 1])
        @test_throws ArgumentError Bogoliubov(a, [1 0 0; 0 1 0])
        @test_throws ArgumentError Bogoliubov((left, left), [1 0; 0 1], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1], [1 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1; 0 0], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1], [0 0 0; 0 0 0])

        @test_throws ArgumentError SecondQuantizedAlgebra.bogoliubov_matrix(
            [1 0 0; 0 1 0], [0 0; 0 0], 2,
        )
        @test_throws ArgumentError SecondQuantizedAlgebra.bogoliubov_matrix(
            [1 0; 0 1], [0 0 0; 0 0 0], 2,
        )

        identity_map = Bogoliubov(a, [1 0; 0 1])
        @test iszero(simplify(conjugate(a, identity_map) - a))
        @test iszero(simplify(conjugate(a, inv(identity_map)) - a))
    end
end
