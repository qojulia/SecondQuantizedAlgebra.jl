using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QMul, QAdd, QSym, HilbertSpace
using Test

@testset "nlevel" begin
    @testset "NLevelSpace" begin
        h = NLevelSpace(:atom, 3, 1)
        @test h isa HilbertSpace
        @test h.name == :atom
        @test h.n == 3
        @test h.ground_state == 1
        @test NLevelSpace(:a, 3, 1) == NLevelSpace(:a, 3, 1)
        @test NLevelSpace(:a, 3, 1) != NLevelSpace(:b, 3, 1)
    end

    @testset "Transition construction — single space" begin
        h = NLevelSpace(:atom, 3, 1)
        σ = Transition(h, :σ, 1, 2)
        @test σ isa Transition
        @test σ isa QSym
        @test σ.name == :σ
        @test σ.i == 1
        @test σ.j == 2
        @test σ.space_index == 1
    end

    @testset "Transition construction — product space" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2, 2)
        @test σ.space_index == 2
        @test_throws ArgumentError Transition(h, :σ, 1, 2, 1)
    end

    @testset "Adjoint" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = σ12'
        @test σ21 isa Transition
        @test σ21.i == 2
        @test σ21.j == 1
        @test σ12'' == σ12
    end

    @testset "Equality and hashing" begin
        h = NLevelSpace(:atom, 3, 1)
        σ1 = Transition(h, :σ, 1, 2)
        σ2 = Transition(h, :σ, 1, 2)
        σ3 = Transition(h, :σ, 2, 1)
        @test isequal(σ1, σ2)
        @test !isequal(σ1, σ3)
        @test hash(σ1) == hash(σ2)
        @test hash(σ1) != hash(σ3)
    end

    @testset "Arithmetic" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)

        m = σ12 * σ21
        @test m isa QMul
        @test length(m.args_nc) == 2

        s = σ12 + σ21
        @test s isa QAdd
    end

    @testset "@qnumbers" begin
        h = NLevelSpace(:atom, 2, 1)
        @qnumbers σ::Transition(h, 1, 2)
        @test σ isa Transition
        @test σ.name == :σ
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(NLevelSpace)
        all_concrete(Transition)
    end

    @testset "Type stability" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)

        @inferred Transition(:σ, 1, 2, 1)
        @inferred adjoint(σ12)
        @inferred isequal(σ12, σ21)
        @inferred hash(σ12, UInt(0))
        @inferred σ12 * σ21
        @inferred σ12 + σ21
    end

    @testset "Allocations" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)

        # Construction and adjoint
        @test @allocations(Transition(:σ, 1, 2, 1)) == 0
        @test @allocations(adjoint(σ12)) == 0

        # Equality / hashing
        σ12b = Transition(h, :σ, 1, 2)
        @test @allocations(isequal(σ12, σ12b)) == 0
        @test @allocations(hash(σ12, UInt(0))) == 0
    end
end
