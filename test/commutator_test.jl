using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: _ZERO_QADD, simplify, QMul, QAdd, QSym

@testset "commutator" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Scalar short-circuits" begin
        @test commutator(1, 2) === _ZERO_QADD
        @test commutator(3, a) === _ZERO_QADD
        @test commutator(a, 5) === _ZERO_QADD
        @test commutator(2.0, ad) === _ZERO_QADD
    end

    @testset "Fock: [a, a†] = 1" begin
        result = commutator(a, ad)
        @test result isa QAdd
        # simplify(a*a' - a'*a) = 1 (scalar identity)
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)
        @test result.arguments[1].arg_c == 1
    end

    @testset "Fock: [a†, a] = -1" begin
        result = commutator(ad, a)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == -1
    end

    @testset "Identical operators: [a, a] = 0" begin
        result = commutator(a, a)
        @test result === _ZERO_QADD
    end

    @testset "Different spaces: [a, b] = 0" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1, a2) === _ZERO_QADD
        @test commutator(a1, a2') === _ZERO_QADD
    end

    @testset "QMul with QSym — shared space" begin
        result = commutator(ad * a, a)
        @test result isa QAdd
    end

    @testset "QMul with QSym — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1' * a1, a2) === _ZERO_QADD
        @test commutator(a2, a1' * a1) === _ZERO_QADD
    end

    @testset "QMul with QMul — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test commutator(a1' * a1, a2' * a2) === _ZERO_QADD
    end

    @testset "QAdd bilinearity: [a + a†, a] = [a,a] + [a†,a]" begin
        expr = a + ad
        result = commutator(expr, a)
        @test result isa QAdd
        # [a,a] = 0, [a†,a] = -1 → total = -1
        simplified = simplify(result)
        @test length(simplified.arguments) == 1
        @test simplified.arguments[1].arg_c == -1
    end

    @testset "QAdd, QAdd bilinearity" begin
        result = commutator(a + ad, a + ad)
        @test result isa QAdd
        # [a+a†, a+a†] = [a,a]+[a,a†]+[a†,a]+[a†,a†] = 0+1-1+0 = 0
        simplified = simplify(result)
        @test all(iszero, simplified.arguments)
    end

    @testset "Nested commutator" begin
        # [a, [a, a†]] = [a, 1] = 0
        inner = commutator(a, ad)
        result = commutator(a, inner)
        @test result isa QAdd
        @test all(iszero, simplify(result).arguments)
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
