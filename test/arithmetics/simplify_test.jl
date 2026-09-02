using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test

@testset "Simplification contracts" begin
    @testset "symbolic coefficients are simplified" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @variables g h₁

        expression = ((g + h₁)^2 - (g^2 + 2g * h₁ + h₁^2)) * a
        @test iszero(simplify(expression))
        @test iszero(simplify((g - g) * a))
        @test isequal(
            Symbolics.expand((g + h₁)^2 * a),
            (g^2 + 2g * h₁ + h₁^2) * a
        )
        @test isequal(SymbolicUtils.simplify(2 * g * a), 2 * g * a)
    end

    @testset "operator identities are stable" begin
        hf = FockSpace(:cavity)
        a = Destroy(hf, :a)
        @test iszero(simplify(a * a' - (a' * a + 1)))
        @test isequal(simplify(a * a'), a * a')
        @test isequal(normal_order(simplify(a * a')), simplify(a * a'))

        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        σ31 = Transition(hn, :σ, 3, 1)
        σ11 = Transition(hn, :σ, 1, 1)
        σ21 = Transition(hn, :σ, 2, 1)
        @test iszero(simplify(σ12 * σ23 - Transition(hn, :σ, 1, 3)))
        @test iszero(simplify(σ12 * σ31))
        @test iszero(simplify(σ12 * σ21 - σ11))
        @test iszero(simplify(σ11 * σ12 - σ12))

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test iszero(simplify(σx * σx - 1))
        @test iszero(simplify(σx * σy - im * σz))

        hs = SpinSpace(:spin)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test iszero(simplify(Sy * Sx - (Sx * Sy - im * Sz)))

        hq = PhaseSpace(:phase)
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        @test iszero(simplify(p * x - (x * p - im)))
    end

    @testset "explicit completeness remains explicit" begin
        hn = NLevelSpace(:atom, 3, 2)
        ground = Transition(hn, :σ, 2, 2)
        @test iszero(simplify(ground) - ground)
        @test iszero(
            simplify(
                expand_completeness(ground) -
                    (1 - Transition(hn, :σ, 1, 1) - Transition(hn, :σ, 3, 3))
            ),
        )
        @test iszero(
            simplify(
                expand_completeness(ground) -
                    expand_completeness(simplify(ground))
            )
        )
    end

    @testset "simplification is idempotent across workflows" begin
        hf = FockSpace(:fock1) ⊗ FockSpace(:fock2)
        a = Destroy(hf, :a, 1)
        b = Destroy(hf, :b, 2)
        expression = (a + a') * (b + b') + commutator(a' * a, b' * b)
        @test isequal(simplify(simplify(expression)), simplify(expression))
        @test isequal(normal_order(normal_order(expression)), normal_order(expression))

        hi = FockSpace(:sites)
        site = Destroy(hi, :a)
        i = Index(hi, :i, 10, hi)
        indexed = Σ(IndexedOperator(site', i) * IndexedOperator(site, i), i)
        @test isequal(simplify(simplify(indexed)), simplify(indexed))
        @test isequal(normal_order(normal_order(indexed)), normal_order(indexed))
    end

    @testset "mixed expressions retain their public semantics" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)
        hamiltonian = a' * σ + a * σ'
        @test iszero(simplify(hamiltonian - normal_order(hamiltonian)))
        @test iszero(simplify(hamiltonian' - qadjoint(hamiltonian)))
    end

    @testset "public simplification boundaries are inferable" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test iszero((@inferred simplify(a)) - a)
        @inferred simplify(a * a')
        @inferred simplify(a + a')
        @inferred normal_order(a * a')
        @inferred Symbolics.expand(a * a')
    end
end
