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
    end
end
