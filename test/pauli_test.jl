using SecondQuantizedAlgebra
using Test

@testset "pauli" begin
    @testset "PauliSpace" begin
        h = PauliSpace(:p)
        @test h isa HilbertSpace
        @test h.name == :p
        @test PauliSpace(:p) == PauliSpace(:p)
        @test PauliSpace(:p) != PauliSpace(:q)
    end

    @testset "Pauli construction — single space" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)
        @test σx isa Pauli
        @test σx isa QSym
        @test σx.axis == 1
        @test σy.axis == 2
        @test σz.axis == 3
        @test σx.space_index == 1
        @test_throws ArgumentError Pauli(h, :σ, 4)
    end

    @testset "Pauli construction — product space" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p)
        σx = Pauli(h, :σ, 1, 2)
        @test σx.space_index == 2
        @test_throws ArgumentError Pauli(h, :σ, 1, 1)
    end

    @testset "Adjoint — Hermitian" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        @test σx' == σx
    end

    @testset "Equality and hashing" begin
        h = PauliSpace(:p)
        σ1 = Pauli(h, :σ, 1)
        σ2 = Pauli(h, :σ, 1)
        σ3 = Pauli(h, :σ, 2)
        @test isequal(σ1, σ2)
        @test !isequal(σ1, σ3)
        @test hash(σ1) == hash(σ2)
        @test hash(σ1) != hash(σ3)
    end

    @testset "Arithmetic" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)

        m = σx * σy
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = σx + σy
        @test s isa QAdd{Int}
    end

    @testset "@qnumbers" begin
        h = PauliSpace(:p)
        @qnumbers σx::Pauli(h, 1)
        @test σx isa Pauli
        @test σx.name == :σx
        @test σx.axis == 1
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(PauliSpace)
        all_concrete(Pauli)
    end

    @testset "Type stability" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)

        @inferred Pauli(:σ, 1, 1)
        @inferred adjoint(σx)
        @inferred isequal(σx, σy)
        @inferred hash(σx, UInt(0))
        @inferred σx * σy
        @inferred σx + σy
    end

    @testset "Allocations" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σx2 = Pauli(h, :σ, 1)

        @test @allocations(Pauli(:σ, 1, 1)) == 0
        @test @allocations(adjoint(σx)) == 0
        @test @allocations(isequal(σx, σx2)) == 0
        @test @allocations(hash(σx, UInt(0))) == 0
    end
end
