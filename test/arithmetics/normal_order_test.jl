using SecondQuantizedAlgebra
using QuantumOpticsBase
using Random: MersenneTwister, rand
using Test

@testset "Normal ordering" begin
    @testset "Fock relations" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        ad = adjoint(a)

        @test iszero(normal_order(ad * a) - ad * a)
        @test iszero(normal_order(a * ad) - (ad * a + 1))
        @test iszero(normal_order(a * ad + ad * a) - (2 * ad * a + 1))
        @test iszero(
            normal_order(a * a * ad * ad) -
                (ad * ad * a * a + 4 * ad * a + 2),
        )
        @test iszero(normal_order(a) - a)
    end

    @testset "same-site operator families" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        σ31 = Transition(hn, :σ, 3, 1)
        @test iszero(normal_order(σ12 * σ23) - Transition(hn, :σ, 1, 3))
        @test iszero(normal_order(σ12 * σ31))

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test iszero(normal_order(σx * σx) - 1)
        @test iszero(normal_order(σx * σy * σz) - im)

        hs = SpinSpace(:spin)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test iszero(normal_order(Sy * Sx) - (Sx * Sy - im * Sz))
        @test iszero(normal_order(Sx * Sy) - Sx * Sy)

        hq = PhaseSpace(:phase)
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        @test iszero(normal_order(p * x) - (x * p - im))
    end

    @testset "mixed-space canonicalization" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 3, 2)
        a = Destroy(h, :a, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        σ23 = Transition(h, :σ, 2, 3, 2)
        σ13 = Transition(h, :σ, 1, 3, 2)

        expression = σ12 * σ23 * a * a'
        expected = σ13 * (a' * a + 1)

        @test iszero(expression - expected)
        @test iszero(normal_order(expression) - expected)
        @test isequal(normal_order(normal_order(expression)), normal_order(expression))
        @test isequal(simplify(simplify(expression)), simplify(expression))
    end

    @testset "completeness and independent sites" begin
        hn = NLevelSpace(:atom, 3, 2)
        ground = Transition(hn, :σ, 2, 2)
        @test iszero(normal_order(ground) - ground)
        @test iszero(
            simplify(
                expand_completeness(normal_order(ground)) -
                    (1 - Transition(hn, :σ, 1, 1) - Transition(hn, :σ, 3, 3))
            ),
        )

        h = FockSpace(:left) ⊗ FockSpace(:right)
        left = Destroy(h, :a, 1)
        right = Destroy(h, :b, 2)
        @test iszero(normal_order(left * right') - left * right')
        @test iszero(normal_order(left * left') - (left' * left + 1))
    end
end

@testset "Weyl ordering round trips" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)

    @testset "single-mode and mixed expressions" begin
        for expression in (a, a' * a, a * a', (a * a')^2, (a * a')^3)
            symmetric = normal_to_symmetric(expression)
            @test iszero(symmetric_to_normal(symmetric) - expression)
            @test iszero(normal_to_symmetric(symmetric_to_normal(symmetric)) - symmetric)
        end

        h2 = FockSpace(:left) ⊗ FockSpace(:right)
        left = Destroy(h2, :a, 1)
        right = Destroy(h2, :b, 2)
        expression = left * right' * left' * right
        @test iszero(symmetric_to_normal(normal_to_symmetric(expression)) - expression)

        hm = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        cavity = Destroy(hm, :a, 1)
        σ = Transition(hm, :σ, 1, 2, 2)
        expression = (cavity' * σ + cavity * σ')^2
        @test iszero(symmetric_to_normal(normal_to_symmetric(expression)) - expression)
    end

    @testset "generated words" begin
        h2 = FockSpace(:left) ⊗ FockSpace(:right)
        left, right = Destroy(h2, :a, 1), Destroy(h2, :b, 2)
        pool = [left, left', right, right']
        rng = MersenneTwister(3)
        for _ in 1:40
            expression = prod(rand(rng, pool, rand(rng, 1:4)))
            @test iszero(symmetric_to_normal(normal_to_symmetric(expression)) - expression)
        end
    end

    @testset "non-Heisenberg operators are unchanged" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ = Transition(hn, :σ, 1, 2)
        @test iszero(normal_to_symmetric(σ) - σ)
        @test iszero(symmetric_to_normal(σ) - σ)

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        @test iszero(normal_to_symmetric(σx) - σx)
        @test iszero(symmetric_to_normal(σx) - σx)
    end

    @testset "public ordering boundaries are inferable" begin
        @inferred normal_order(a * a')
        @inferred normal_to_symmetric(a' * a)
        @inferred symmetric_to_normal(a' * a)
        @inferred symmetric_to_normal(normal_to_symmetric(a' * a))

        basis = FockBasis(4)
        @test to_numeric(normal_order(a * a'), basis) isa AbstractOperator
    end
end
