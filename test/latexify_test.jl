using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Test

@testset "latexify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Operators" begin
        @test latexify(a) == L"a"
        @test occursin("dagger", latexify(ad))
        @test repr(MIME"text/latex"(), a) == latexify(a)
        @test repr(MIME"text/latex"(), ad) == latexify(ad)
    end

    @testset "QMul" begin
        s = latexify(3 * ad * a)
        @test occursin("3", s)
        @test occursin("dagger", s)

        # Unit prefactor omitted
        s1 = latexify(1 * ad * a)
        @test !startswith(s1, "\$1")

        # Scalar QMul
        @test latexify(QMul(5, QSym[])) == L"5"

        # Negative prefactor
        s_neg = latexify(-1 * a)
        @test occursin("-", s_neg)
    end

    @testset "QAdd" begin
        s = latexify(ad * a + 1)
        @test occursin("dagger", s)
        @test occursin("+", s)

        # Multiple terms
        s2 = latexify(2 * a + 3 * ad)
        @test occursin("2", s2)
        @test occursin("3", s2)
    end

    @testset "MIME text/latex" begin
        @test repr(MIME"text/latex"(), a * ad) == latexify(a * ad)
        @test repr(MIME"text/latex"(), a + ad) == latexify(a + ad)
    end
end
