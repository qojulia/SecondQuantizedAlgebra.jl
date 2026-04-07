using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QAdd, QSym, _CNUM_ONE
using Test

@testset "Multiplication (eager ordering)" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Eager normal ordering" begin
        # a * a† = a†a + 1 (normal ordered)
        result = a * ad
        @test result[QSym[ad, a]] == _CNUM_ONE
        @test result[QSym[]] == _CNUM_ONE
        @test length(result) == 2

        # a† * a is already normal ordered
        result2 = ad * a
        @test result2[QSym[ad, a]] == _CNUM_ONE
        @test length(result2) == 1
    end

    @testset "Cross-site operators" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2
        @test m isa QAdd
        @test length(m) == 1
        ops = first(keys(m.arguments))
        @test ops[1].space_index <= ops[2].space_index
    end

    @testset "Division" begin
        m = a / 2
        @test m isa QAdd
        @test m[QSym[a]] == 1 // 2
    end

    @testset "Power" begin
        m = a^3
        @test m isa QAdd
        @test length(first(keys(m.arguments))) == 3

        @test length(a^0) == 1  # scalar 1
        @test (a^0)[QSym[]] == _CNUM_ONE
    end

    @testset "Negation" begin
        m = -a
        @test m isa QAdd
        @test m[QSym[a]] == -1
    end

    @testset "Adjoint" begin
        m = 2 * ad * a
        ma = m'
        @test ma isa QAdd
        # (2 a†a)† = 2 a†a
        @test isequal(ma, m)
    end

    @testset "Type stability" begin
        @inferred a * ad
        @inferred 3 * a
        @inferred a * 2.0
        @inferred a / 2
        @inferred a^3
        @inferred -a
    end
end
