using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym

@testset "normal_order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Already normal-ordered" begin
        m = ad * a
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isequal(result.arguments[1], m)
    end

    @testset "Single commutation: a*a† → a†a + 1" begin
        m = a * ad
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "QSym passthrough" begin
        result = normal_order(a)
        @test result isa QAdd
        @test length(result.arguments) == 1
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
        @test length(result.arguments) == 1
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
        m = σ12 * σ23
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        op = result.arguments[1].args_nc[1]
        @test op isa Transition
        @test op.i == 1 && op.j == 3
    end

    @testset "Orthogonality: |1⟩⟨2| · |3⟩⟨1| = 0" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        m = σ12 * σ31
        result = normal_order(m)
        @test result isa QAdd
        @test all(iszero, result.arguments)
    end

    @testset "Orthogonality in longer product" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        σ23 = Transition(h, :σ, 2, 3)
        # σ12 * σ31 = 0, so σ12 * σ31 * σ23 = 0
        m = σ12 * σ31 * σ23
        result = normal_order(m)
        @test result isa QAdd
        @test all(iszero, result.arguments)
    end

    @testset "Ground state rewriting — single space" begin
        σ11 = Transition(h, :σ, 1, 1)
        m = QMul(1, QSym[σ11])
        result = normal_order(m, h)
        # |1⟩⟨1| = 1 - |2⟩⟨2| - |3⟩⟨3| for ground_state=1, n=3
        @test result isa QAdd
        @test length(result.arguments) == 3
        # Verify content: one scalar term (1), two diagonal transitions with -1
        simplified = simplify(result)
        @test length(simplified.arguments) == 3
    end

    @testset "Ground state rewriting — ProductSpace" begin
        hc = FockSpace(:c)
        ha = NLevelSpace(:atom, 2, 1)
        hp = hc ⊗ ha
        σ11 = Transition(hp, :σ, 1, 1, 2)
        m = QMul(1, QSym[σ11])
        result = normal_order(m, hp)
        # |1⟩⟨1| = 1 - |2⟩⟨2| for ground_state=1, n=2
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "Ground state rewriting — recursive (double projection)" begin
        h2 = NLevelSpace(:atom, 2, 1)
        σ11 = Transition(h2, :σ, 1, 1)
        # σ11 * σ11 — both get ground-state expanded
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
        m = σx * σx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)
        @test result.arguments[1].arg_c == 1
    end

    @testset "σx·σy = iσz" begin
        m = σx * σy
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        t = result.arguments[1]
        @test t.arg_c == im
        @test length(t.args_nc) == 1
        @test t.args_nc[1] isa Pauli
        @test t.args_nc[1].axis == 3
    end

    @testset "σy·σx = -iσz" begin
        m = σy * σx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        t = result.arguments[1]
        @test t.arg_c == -im
        @test t.args_nc[1].axis == 3
    end

    @testset "Full Pauli cycle" begin
        # σy·σz = iσx
        result = normal_order(σy * σz)
        @test result.arguments[1].arg_c == im
        @test result.arguments[1].args_nc[1].axis == 1

        # σz·σx = iσy
        result = normal_order(σz * σx)
        @test result.arguments[1].arg_c == im
        @test result.arguments[1].args_nc[1].axis == 2

        # σz·σy = -iσx
        result = normal_order(σz * σy)
        @test result.arguments[1].arg_c == -im
        @test result.arguments[1].args_nc[1].axis == 1

        # σx·σz = -iσy
        result = normal_order(σx * σz)
        @test result.arguments[1].arg_c == -im
        @test result.arguments[1].args_nc[1].axis == 2
    end

    @testset "σy·σy = 1, σz·σz = 1" begin
        @test normal_order(σy * σy).arguments[1].arg_c == 1
        @test normal_order(σz * σz).arguments[1].arg_c == 1
    end
end

@testset "normal_order — Spin" begin
    h = SpinSpace(:s, 1 // 2)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)

    @testset "Sy·Sx gets reordered (swap + commutator)" begin
        m = Sy * Sx  # out of axis order (2 > 1)
        result = normal_order(m)
        @test result isa QAdd
        # SySx = SxSy - iSz → two terms
        @test length(result.arguments) >= 2
    end

    @testset "Same axis — already ordered" begin
        m = Sx * Sx
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
    end

    @testset "Sx·Sy already in order" begin
        m = Sx * Sy  # axes 1,2 — already in order
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
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

    hs = SpinSpace(:s, 1 // 2)
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

    # Warmup
    normal_order(a * a')
    normal_order(a' * a)
    normal_order(a * a' + a' * a)

    @test @allocations(normal_order(a * a')) < 1000
    @test @allocations(normal_order(a' * a)) < 1000
    @test @allocations(normal_order(a * a' + a' * a)) < 2000

    # Transition with ground state rewriting
    hn = NLevelSpace(:atom, 3, 1)
    σ11 = Transition(hn, :σ, 1, 1)
    normal_order(QMul(1, QSym[σ11]), hn)
    @test @allocations(normal_order(QMul(1, QSym[σ11]), hn)) < 5000

    # PhaseSpace
    hps = PhaseSpace(:q)
    x = Position(hps, :x)
    p = Momentum(hps, :p)
    normal_order(p * x)
    @test @allocations(normal_order(p * x)) < 1000
end
