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
        la = latexify(ad)
        @test occursin("dagger", la)
    end

    @testset "QMul" begin
        s = latexify(3 * ad * a)
        @test occursin("3", s)
        @test occursin("dagger", s)

        @test latexify(QMul(5, QSym[])) == L"5"
    end

    @testset "QAdd" begin
        s = latexify(ad * a + 1)
        @test occursin("dagger", s)
        @test occursin("+", s)
    end
end
