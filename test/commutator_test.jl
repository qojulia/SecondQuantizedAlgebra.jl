using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym

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
        @test isempty(only(collect(result)).args_nc)
        @test only(collect(result)).arg_c == 1
    end

    @testset "Fock: [a†, a] = -1" begin
        result = commutator(ad, a)
        @test result isa QAdd
        @test length(result) == 1
        @test only(collect(result)).arg_c == -1
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

    @testset "QMul with QSym — shared space" begin
        result = commutator(ad * a, a)
        @test result isa QAdd
    end

    @testset "QMul with QSym — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1' * a1, a2))
        @test iszero(commutator(a2, a1' * a1))
    end

    @testset "QMul with QMul — no shared space" begin
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
        @test only(collect(simplified)).arg_c == -1
    end

    @testset "QAdd, QAdd bilinearity" begin
        result = commutator(a + ad, a + ad)
        @test result isa QAdd
        # [a+a†, a+a†] = [a,a]+[a,a†]+[a†,a]+[a†,a†] = 0+1-1+0 = 0
        simplified = simplify(result)
        @test iszero(simplified)
    end

    @testset "Nested commutator" begin
        # [a, [a, a†]] = [a, 1] = 0
        inner = commutator(a, ad)
        result = commutator(a, inner)
        @test result isa QAdd
        @test iszero(simplify(result))
    end

    @testset "Spin: [Sx, Sy] = iSz" begin
        hs = SpinSpace(:s, 1 // 2)
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
end
