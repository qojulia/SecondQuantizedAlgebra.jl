using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym, QField

@testset "simplify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Collect like terms" begin
        s = QAdd(QMul[QMul(2, QSym[ad, a]), QMul(3, QSym[ad, a])])
        result = simplify(s)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 5
    end

    @testset "Remove zero terms" begin
        s = QAdd(QMul[QMul(0, QSym[ad, a]), QMul(3, QSym[ad, a])])
        result = simplify(s)
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 3
    end

    @testset "a + a = 2a" begin
        s = a + a
        result = simplify(s)
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 2
    end

    @testset "Symbolic prefactors" begin
        @variables g h_sym
        s = QAdd([QMul(g, QSym[ad, a]), QMul(h_sym, QSym[ad, a])])
        result = simplify(s)
        @test length(result.arguments) == 1
    end

    @testset "Fock: simplify applies [a, a†] = 1" begin
        m = a * ad  # a·a† stored lazily
        result = simplify(m)
        @test result isa QAdd
        @test length(result.arguments) == 2  # a†a + 1
    end

    @testset "Transition: simplify applies composition" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        result = simplify(σ12 * σ23)
        @test length(result.arguments) == 1
        @test result.arguments[1].args_nc[1].i == 1
        @test result.arguments[1].args_nc[1].j == 3
    end

    @testset "Transition: simplify applies orthogonality" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ31 = Transition(hn, :σ, 3, 1)
        result = simplify(σ12 * σ31)
        @test all(iszero, result.arguments)
    end

    @testset "Pauli: simplify applies σⱼ²=I and product rule" begin
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)

        # σx² = I
        result = simplify(σx * σx)
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)
        @test result.arguments[1].arg_c == 1

        # σx·σy = iσz
        result = simplify(σx * σy)
        @test result.arguments[1].arg_c == im
        @test result.arguments[1].args_nc[1].axis == 3
    end

    @testset "simplify with ordering argument" begin
        result = simplify(a * ad)
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "SymbolicUtils.simplify on QField" begin
        @variables g
        s = QAdd([QMul(g, QSym[ad, a]), QMul(g, QSym[ad, a])])
        result = SymbolicUtils.simplify(s)
        @test result isa QAdd
    end

    @testset "Symbolics.expand on QField" begin
        expr = (a + ad) * (a + ad)
        result = Symbolics.expand(expr)
        @test result isa QAdd
        @test length(result.arguments) == 4
    end

    @testset "Return type is QAdd" begin
        @test simplify(a + ad) isa QAdd
        @test simplify(a) isa QAdd
        @test simplify(2 * a) isa QAdd
    end

    @testset "Type stability" begin
        @inferred simplify(a)
        @inferred simplify(a * ad)
        @inferred simplify(a + ad)
        @inferred simplify(a * ad)
        @inferred simplify(ad * a)

        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @inferred simplify(σ12 * σ23)
        @inferred simplify(σ12 + σ23)
        @inferred simplify(σ12, hn)

        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        @inferred simplify(σx * σy)
        @inferred simplify(σx * σx)
        @inferred simplify(σx + σy)

        hs = SpinSpace(:s, 1 // 2)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @inferred simplify(Sz * Sy)
        @inferred simplify(Sx + Sy)

        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @inferred simplify(p * x)
        @inferred simplify(x * p)
        @inferred simplify(x + p)
    end

    @testset "Allocations — Fock" begin
        # Warmup
        simplify(a * ad)
        simplify(ad * a)
        simplify(a + ad)

        allocs_swap = @allocations simplify(a * ad)
        allocs_ordered = @allocations simplify(ad * a)
        allocs_add = @allocations simplify(a + ad)
        allocs_sym = @allocations simplify(a)

        @test allocs_swap < 1000
        @test allocs_ordered < 1000
        @test allocs_add < 1000
        @test allocs_sym < 200

        # Scaling: (a*a')^n allocations should grow, but not explosively
        simplify((a * ad)^2)
        simplify((a * ad)^3)
        allocs_2 = @allocations simplify((a * ad)^2)
        allocs_3 = @allocations simplify((a * ad)^3)
        @test allocs_2 < 3000
        @test allocs_3 < 10000
    end

    @testset "Allocations — Transition" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        σ31 = Transition(hn, :σ, 3, 1)

        # Warmup
        simplify(σ12 * σ23)
        simplify(σ12 * σ31)
        simplify(σ12 + σ23)

        @test @allocations(simplify(σ12 * σ23)) < 1000  # composition
        @test @allocations(simplify(σ12 * σ31)) < 1000  # orthogonality → 0
        @test @allocations(simplify(σ12 + σ23)) < 1000
    end

    @testset "Allocations — Pauli" begin
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)

        # Warmup
        simplify(σx * σx)
        simplify(σx * σy)
        simplify(σx * σy * σz)

        @test @allocations(simplify(σx * σx)) < 1000   # σ² = I
        @test @allocations(simplify(σx * σy)) < 1000   # product rule
        @test @allocations(simplify(σx * σy * σz)) < 1000
    end

    @testset "Allocations — Spin" begin
        hs = SpinSpace(:s, 1 // 2)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)

        # Warmup
        simplify(Sz * Sy)
        simplify(Sz * Sy * Sx)

        @test @allocations(simplify(Sz * Sy)) < 3000
        @test @allocations(simplify(Sz * Sy * Sx)) < 3000
    end

    @testset "Allocations — PhaseSpace" begin
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)

        # Warmup
        simplify(p * x)
        simplify(x * p)

        @test @allocations(simplify(p * x)) < 1000   # P·X → X·P - i
        @test @allocations(simplify(x * p)) < 1000   # already ordered
    end

    @testset "Allocations — collect like terms" begin
        # Many duplicate terms: Dict-based collection should handle efficiently
        many = sum(a' * a for _ in 1:20) + sum(a * a' for _ in 1:20)
        simplify(many)  # warmup
        @test @allocations(simplify(many)) < 50000
    end
end
