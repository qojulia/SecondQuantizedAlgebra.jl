using SecondQuantizedAlgebra
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
    end

    @testset "QMul" begin
        @test repr(1 * ad * a) == "a† * a"
        @test repr(-1 * a) == "-a"
        @test repr(QMul(5, QSym[])) == "5"
    end

    @testset "QAdd" begin
        s = ad * a + 1
        r = repr(s)
        @test occursin("a†", r)
        @test occursin("+", r) || occursin("-", r)
    end
end
