using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test
import SecondQuantizedAlgebra: simplify

@testset "simplify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Collect like terms" begin
        s = QAdd(QMul{Int}[QMul(2, QSym[ad, a]), QMul(3, QSym[ad, a])])
        result = simplify(s)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test result.arguments[1].arg_c == 5
    end

    @testset "Remove zero terms" begin
        s = QAdd(QMul{Int}[QMul(0, QSym[ad, a]), QMul(3, QSym[ad, a])])
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
        result = simplify(a * ad, NormalOrder())
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
end
