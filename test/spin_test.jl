using SecondQuantizedAlgebra
using Test

@testset "spin" begin
    @testset "SpinSpace" begin
        h = SpinSpace(:s, 1 // 2)
        @test h isa HilbertSpace
        @test h.name == :s
        @test h.spin == 1 // 2
        @test SpinSpace(:s, 1 // 2) == SpinSpace(:s, 1 // 2)
        @test SpinSpace(:s, 1 // 2) != SpinSpace(:s, 1)
    end

    @testset "Spin construction — single space" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        @test Sx isa Spin
        @test Sx isa QSym
        @test Sx.axis == 1
        @test Sy.axis == 2
        @test Sz.axis == 3
        @test Sx.space_index == 1
        @test_throws ArgumentError Spin(h, :S, 4)
    end

    @testset "Spin construction — product space" begin
        h = FockSpace(:c) ⊗ SpinSpace(:s, 1)
        Sx = Spin(h, :S, 1, 2)
        @test Sx.space_index == 2
        @test_throws ArgumentError Spin(h, :S, 1, 1)
    end

    @testset "Adjoint — Hermitian" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        @test Sx' == Sx
    end

    @testset "Equality and hashing" begin
        h = SpinSpace(:s, 1 // 2)
        S1 = Spin(h, :S, 1)
        S2 = Spin(h, :S, 1)
        S3 = Spin(h, :S, 2)
        @test isequal(S1, S2)
        @test !isequal(S1, S3)
        @test hash(S1) == hash(S2)
        @test hash(S1) != hash(S3)
    end

    @testset "Arithmetic" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        m = Sx * Sy
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = Sx + Sy
        @test s isa QAdd{Int}
    end

    @testset "@qnumbers" begin
        h = SpinSpace(:s, 1 // 2)
        @qnumbers Sx::Spin(h, 1)
        @test Sx isa Spin
        @test Sx.name == :Sx
        @test Sx.axis == 1
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(SpinSpace)
        all_concrete(Spin)
    end

    @testset "Type stability" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        @inferred Spin(:S, 1, 1)
        @inferred adjoint(Sx)
        @inferred isequal(Sx, Sy)
        @inferred hash(Sx, UInt(0))
        @inferred Sx * Sy
        @inferred Sx + Sy
    end

    @testset "Allocations" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sx2 = Spin(h, :S, 1)

        @test @allocations(Spin(:S, 1, 1)) == 0
        @test @allocations(adjoint(Sx)) == 0
        @test @allocations(isequal(Sx, Sx2)) == 0
        @test @allocations(hash(Sx, UInt(0))) == 0
    end
end
