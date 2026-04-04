using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "printing" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = a'

    @testset "HilbertSpace" begin
        @test repr(h) == "ℋ(cavity)"
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        @test repr(h2) == "ℋ(a) ⊗ ℋ(b)"
    end

    @testset "Operators" begin
        @test repr(a) == "a"
        @test repr(ad) == "a†"
        @test sprint(show, a) == "a"
        @test sprint(show, ad) == "a†"
    end

    @testset "QMul" begin
        @test repr(1 * ad * a) == "a† * a"
        @test repr(3 * ad * a) == "3 * a† * a"
        @test repr(-1 * a) == "-a"
        @test repr(-3 * a) == "-3 * a"
        @test repr(QMul(5, QSym[])) == "5"
        @test repr(0.5 * a) == "0.5 * a"
    end

    @testset "QAdd" begin
        @test repr(ad * a + 1) == "a† * a + 1"
        @test repr(a + ad) == "a + a†"
        @test repr(a + (-1 * ad)) == "a - a†"
        @test repr(2 * a + 3 * ad) == "2 * a + 3 * a†"
    end

    @testset "Symbolic prefactors" begin
        @variables g ω
        @test repr(g * a) isa String
        @test occursin("g", repr(g * a))
        @test repr(g * a + ω * ad) isa String

        # Compound symbolic prefactors get parenthesized
        s = repr(QMul(g + ω, QSym[a]))
        @test occursin("(", s)
    end

    @testset "Type inference" begin
        s = IOBuffer(sizehint=0)
        @inferred show(s, a)
        @inferred show(s, ad)
    end
end
