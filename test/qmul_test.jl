using SecondQuantizedAlgebra
using Test

@testset "QMul" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Construction" begin
        m = a * ad
        @test m isa QMul{Int}
        @test m.arg_c == 1
        @test length(m.args_nc) == 2
        # Same space: order preserved (a then a†)
        @test m.args_nc[1] isa Destroy
        @test m.args_nc[2] isa Create
    end

    @testset "Scalar multiplication" begin
        m1 = 3 * a
        @test m1 isa QMul{Int}
        @test m1.arg_c == 3
        @test m1.args_nc == [a]

        m2 = a * 2.0
        @test m2 isa QMul{Float64}
        @test m2.arg_c == 2.0

        m3 = 0 * a
        @test m3 isa QMul{Int}
        @test m3.arg_c == 0
    end

    @testset "QSym * QSym" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2
        @test m isa QMul{Int}
        @test m.args_nc[1].space_index <= m.args_nc[2].space_index
    end

    @testset "QMul * QSym and QSym * QMul" begin
        m1 = (2 * a) * ad
        @test m1 isa QMul{Int}
        @test m1.arg_c == 2
        @test length(m1.args_nc) == 2

        m2 = ad * (3 * a)
        @test m2 isa QMul{Int}
        @test m2.arg_c == 3
        @test length(m2.args_nc) == 2
    end

    @testset "QMul * QMul" begin
        m1 = 2 * a
        m2 = 3 * ad
        m3 = m1 * m2
        @test m3 isa QMul{Int}
        @test m3.arg_c == 6
        @test length(m3.args_nc) == 2
    end

    @testset "QMul * Number" begin
        m = (2 * a) * 3
        @test m isa QMul{Int}
        @test m.arg_c == 6
    end

    @testset "Division" begin
        m = a / 2
        @test m isa QMul
        @test m.arg_c == 1 // 2
    end

    @testset "Power" begin
        m = a^3
        @test m isa QMul
        @test length(m.args_nc) == 3
    end

    @testset "Negation" begin
        m = -a
        @test m isa QMul{Int}
        @test m.arg_c == -1
    end

    @testset "Lazy — no commutation" begin
        m = a * ad
        @test m isa QMul
        @test length(m.args_nc) == 2
    end

    @testset "Equality and hashing" begin
        m1 = 2 * a * ad
        m2 = 2 * a * ad
        @test isequal(m1, m2)
        @test hash(m1) == hash(m2)
    end

    @testset "iszero" begin
        m = 0 * a
        @test iszero(m)
        @test !iszero(2 * a)
    end

    @testset "Adjoint" begin
        m = 2 * ad * a
        ma = m'
        @test ma isa QMul
        @test ma.arg_c == 2
        @test ma.args_nc[1] isa Create
        @test ma.args_nc[2] isa Destroy
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(QMul)
    end

    @testset "Type stability" begin
        @inferred a * ad
        @inferred 3 * a
        @inferred a * 2.0
        @inferred (2 * a) * ad
        @inferred (2 * a) * (3 * ad)
    end
end
