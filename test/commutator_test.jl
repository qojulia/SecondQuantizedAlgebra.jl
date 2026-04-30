using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QTermDict, _to_cnum, Index, sorted_arguments

@testset "commutator" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Scalar short-circuits" begin
        @test iszero(commutator(1, 2))
        @test iszero(commutator(3, a))
        @test iszero(commutator(a, 5))
        @test iszero(commutator(2.0, ad))
    end

    @testset "Fock: [a, a†] = 1" begin
        result = commutator(a, ad)
        @test result isa QAdd
        # simplify(a*a' - a'*a) = 1 (scalar identity)
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))
        @test prefactor(only(sorted_arguments(result))) == 1
    end

    @testset "Fock: [a†, a] = -1" begin
        result = commutator(ad, a)
        @test result isa QAdd
        @test length(result) == 1
        @test prefactor(only(sorted_arguments(result))) == -1
    end

    @testset "Identical operators: [a, a] = 0" begin
        @test iszero(commutator(a, a))
    end

    @testset "Different spaces: [a, b] = 0" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1, a2))
        @test iszero(commutator(a1, a2'))
    end

    @testset "QAdd with QSym — shared space" begin
        result = commutator(ad * a, a)
        @test result isa QAdd
    end

    @testset "QAdd with QSym — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1' * a1, a2))
        @test iszero(commutator(a2, a1' * a1))
    end

    @testset "QAdd with QAdd — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1' * a1, a2' * a2))
    end

    @testset "QAdd bilinearity: [a + a†, a] = [a,a] + [a†,a]" begin
        expr = a + ad
        result = commutator(expr, a)
        @test result isa QAdd
        # [a,a] = 0, [a†,a] = -1 → total = -1
        simplified = simplify(result)
        @test length(simplified) == 1
        @test prefactor(only(sorted_arguments(simplified))) == -1
    end

    @testset "QAdd, QAdd bilinearity" begin
        result = commutator(a + ad, a + ad)
        @test result isa QAdd
        # [a+a†, a+a†] = [a,a]+[a,a†]+[a†,a]+[a†,a†] = 0+1-1+0 = 0
        simplified = simplify(result)
        @test iszero(simplified)
    end

    @testset "Sum-sum commutator preserves indices" begin
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(ad, j)

        s1 = Σ(ai, i)
        s2 = Σ(adj, j)
        result = commutator(s1, s2)
        @test result isa QAdd
        @test Set(result.indices) == Set([i, j])
    end

    @testset "Nested commutator" begin
        # [a, [a, a†]] = [a, 1] = 0
        inner = commutator(a, ad)
        result = commutator(a, inner)
        @test result isa QAdd
        @test iszero(simplify(result))
    end

    @testset "Spin: [Sx, Sy] = iSz" begin
        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        result = commutator(Sx, Sy)
        @test result isa QAdd
    end

    @testset "PhaseSpace: [X, P] = i" begin
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        result = commutator(x, p)
        @test result isa QAdd
    end

    @testset "Return type is always QAdd" begin
        @test commutator(a, ad) isa QAdd
        @test commutator(a, a) isa QAdd
        @test commutator(1, a) isa QAdd
        @test commutator(a + ad, a) isa QAdd
        @test commutator(a + ad, a + ad) isa QAdd
        @test commutator(ad * a, a) isa QAdd
    end

    @testset "Type stability" begin
        @inferred commutator(a, ad)
        @inferred commutator(ad, a)
        @inferred commutator(a, a)
        @inferred commutator(1, a)
        @inferred commutator(a, 1)
        @inferred commutator(ad * a, a)
        @inferred commutator(a, ad * a)
        @inferred commutator(ad * a, a * ad)
        @inferred commutator(a + ad, a)
        @inferred commutator(a, a + ad)
        @inferred commutator(a + ad, a + ad)

        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @inferred commutator(a1, a2)
        @inferred commutator(a1, a2')
    end

    @testset "Allocations" begin
        # Warmup
        commutator(a, ad)
        commutator(a, a)
        commutator(ad * a, a)
        commutator(a + ad, a)

        # Short-circuit (zero) should be very cheap
        @test @allocations(commutator(a, a)) < 100
        @test @allocations(commutator(1, a)) < 100

        # Non-trivial but simple
        @test @allocations(commutator(a, ad)) < 20000
        @test @allocations(commutator(ad * a, a)) < 30000

        # Bilinearity over QAdd
        @test @allocations(commutator(a + ad, a)) < 40000
    end

    @testset "anticommutator" begin
        # Fock: {a, a†} = a a† + a† a = (a†a + 1) + a†a = 1 + 2 a†a
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'
        r = anticommutator(a, ad)
        @test r isa QAdd
        @test isequal(r, 1 + 2 * ad * a)

        # {a, a} = 2 a²
        @test isequal(anticommutator(a, a), 2 * a * a)

        # Pauli: {σⱼ, σₖ} = 2 δⱼₖ
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test iszero(anticommutator(σx, σy))
        @test iszero(anticommutator(σx, σz))
        @test iszero(anticommutator(σy, σz))
        # σx² = I, so {σx, σx} = 2
        @test isequal(simplify(anticommutator(σx, σx)), QAdd(QTermDict(QSym[] => _to_cnum(2)), Index[], Tuple{Index, Index}[]))

        # Different sites: {A, B} = 2 A B
        h2 = FockSpace(:c1) ⊗ FockSpace(:c2)
        @qnumbers a1::Destroy(h2, 1) a2::Destroy(h2, 2)
        @test isequal(anticommutator(a1, a2), 2 * a1 * a2)

        # Scalars
        @test anticommutator(2, 3) == 12
        @test isequal(anticommutator(2, a), 4 * a)
        @test isequal(anticommutator(a, 3), 6 * a)
    end
end
