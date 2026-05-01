using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum

@testset "normal_order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Already normal-ordered" begin
        m = ad * a
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
        @test isequal(result, m)
    end

    @testset "Single commutation: a*a† → a†a + 1" begin
        m = a * ad
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "QSym passthrough" begin
        result = normal_order(a)
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "QAdd — normal-orders each term" begin
        expr = a * ad + ad * a
        result = normal_order(expr)
        @test result isa QAdd
    end

    @testset "Different spaces don't commute" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "Return type is always QAdd" begin
        @test normal_order(a) isa QAdd
        @test normal_order(2 * a) isa QAdd
        @test normal_order(a + ad) isa QAdd
    end
end

@testset "normal_order — Transition" begin
    h = NLevelSpace(:atom, 3, 1)

    @testset "Composition: |1⟩⟨2| · |2⟩⟨3| = |1⟩⟨3|" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ23 = Transition(h, :σ, 2, 3)
        result = σ12 * σ23  # eager composition
        @test result isa QAdd
        @test length(result) == 1
        ops = first(keys(result.arguments))
        @test ops[1] isa Transition
        @test ops[1].i == 1 && ops[1].j == 3
    end

    @testset "Orthogonality: |1⟩⟨2| · |3⟩⟨1| = 0" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        result = σ12 * σ31
        @test iszero(result)
    end

    @testset "Orthogonality in longer product" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        σ23 = Transition(h, :σ, 2, 3)
        result = σ12 * σ31 * σ23
        @test iszero(result)
    end

    @testset "Ground state rewriting — single space" begin
        σ11 = Transition(h, :σ, 1, 1)
        result = normal_order(σ11, h)
        @test result isa QAdd
        @test length(result) == 3
        simplified = simplify(result)
        @test length(simplified) == 3
    end

    @testset "Ground state rewriting — ProductSpace" begin
        hc = FockSpace(:c)
        ha = NLevelSpace(:atom, 2, 1)
        hp = hc ⊗ ha
        σ11 = Transition(hp, :σ, 1, 1, 2)
        result = normal_order(σ11, hp)
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "Ground state rewriting — recursive (double projection)" begin
        h2 = NLevelSpace(:atom, 2, 1)
        σ11 = Transition(h2, :σ, 1, 1)
        m = σ11 * σ11
        result = normal_order(m, h2)
        @test result isa QAdd
    end
end

@testset "normal_order — Pauli" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)

    @testset "σx·σx = 1" begin
        result = σx * σx  # eager: identity
        @test result isa QAdd
        @test length(result) == 1
        ops, c = only(collect(result))
        @test isempty(ops)
        @test c == 1
    end

    @testset "σx·σy = iσz" begin
        result = σx * σy  # eager: iσz
        @test result isa QAdd
        @test length(result) == 1
        ops, c = only(collect(result))
        @test c == im
        @test length(ops) == 1
        @test ops[1] isa Pauli
        @test ops[1].axis == 3
    end

    @testset "σy·σx = -iσz" begin
        result = σy * σx
        @test length(result) == 1
        ops, c = only(collect(result))
        @test c == -im
        @test ops[1].axis == 3
    end

    @testset "Full Pauli cycle" begin
        # σy·σz = iσx
        ops, c = only(collect(σy * σz))
        @test c == im
        @test ops[1].axis == 1

        # σz·σx = iσy
        ops, c = only(collect(σz * σx))
        @test c == im
        @test ops[1].axis == 2

        # σz·σy = -iσx
        ops, c = only(collect(σz * σy))
        @test c == -im
        @test ops[1].axis == 1

        # σx·σz = -iσy
        ops, c = only(collect(σx * σz))
        @test c == -im
        @test ops[1].axis == 2
    end

    @testset "σy·σy = 1, σz·σz = 1" begin
        @test only(collect(σy * σy)).second == 1
        @test only(collect(σz * σz)).second == 1
    end
end

@testset "normal_order — Spin" begin
    h = SpinSpace(:s)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)

    @testset "Sy·Sx gets reordered (swap + commutator)" begin
        result = Sy * Sx  # eager ordering
        @test result isa QAdd
        @test length(result) >= 2
    end

    @testset "Same axis — already ordered" begin
        result = Sx * Sx
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "Sx·Sy already in order" begin
        result = Sx * Sy
        @test result isa QAdd
        @test length(result) == 1
    end
end

@testset "normal_order — Type stability" begin
    hf = FockSpace(:c)
    a = Destroy(hf, :a)

    @inferred normal_order(a)
    @inferred normal_order(a * a')
    @inferred normal_order(a' * a)
    @inferred normal_order(a * a' + a' * a)

    hn = NLevelSpace(:atom, 3, 1)
    σ12 = Transition(hn, :σ, 1, 2)
    σ23 = Transition(hn, :σ, 2, 3)
    @inferred normal_order(σ12 * σ23)
    @inferred normal_order(σ12, hn)

    hp = PauliSpace(:p)
    σx = Pauli(hp, :σ, 1)
    σy = Pauli(hp, :σ, 2)
    @inferred normal_order(σx * σy)

    hs = SpinSpace(:s)
    Sy = Spin(hs, :S, 2)
    Sx = Spin(hs, :S, 1)
    @inferred normal_order(Sy * Sx)

    hps = PhaseSpace(:q)
    x = Position(hps, :x)
    p = Momentum(hps, :p)
    @inferred normal_order(p * x)
end

@testset "normal_order — Allocations" begin
    hf = FockSpace(:c)
    a = Destroy(hf, :a)

    normal_order(a * a')
    normal_order(a' * a)
    normal_order(a * a' + a' * a)

    @test @allocations(normal_order(a * a')) < 1000
    @test @allocations(normal_order(a' * a)) < 1000
    @test @allocations(normal_order(a * a' + a' * a)) < 2000

    hn = NLevelSpace(:atom, 3, 1)
    σ11 = Transition(hn, :σ, 1, 1)
    normal_order(σ11, hn)
    @test @allocations(normal_order(σ11, hn)) < 5000

    hps = PhaseSpace(:q)
    x = Position(hps, :x)
    p = Momentum(hps, :p)
    normal_order(p * x)
    @test @allocations(normal_order(p * x)) < 1000
end
