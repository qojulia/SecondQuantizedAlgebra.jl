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
        @test occursin("-", repr(a + (-1 * ad)))  # order may vary
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

    @testset "NLevel printing" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        @test repr(σ12) == "σ₁₂"
        @test repr(hn) == "ℋ(atom)"
    end

    @testset "Pauli printing" begin
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test repr(σx) == "σx"
        @test repr(σy) == "σy"
        @test repr(σz) == "σz"
        @test repr(hp) == "ℋ(p)"
    end

    @testset "Spin printing" begin
        hs = SpinSpace(:s, 1 // 2)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test repr(Sx) == "Sx"
        @test repr(Sy) == "Sy"
        @test repr(Sz) == "Sz"
        @test repr(hs) == "ℋ(s)"
    end
end
