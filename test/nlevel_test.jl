using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, HilbertSpace
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
        @test m isa QAdd

        s = σ12 + σ21
        @test s isa QAdd
    end

    @testset "@qnumbers" begin
        h = NLevelSpace(:atom, 2, 1)
        @qnumbers σ::Transition(h, 1, 2)
        @test σ isa Transition
        @test σ.name == :σ
    end

    @testset "Ground state projector" begin
        h = NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2)
        σee = Transition(h, :σ, 2, 2)
        σgg = Transition(h, :σ, 1, 1)

        # σ*σ' = σ_gg (ground state projector)
        @test isequal(simplify(normal_order(σ * σ')), simplify(σgg))
        # σ_gg = 1 - σ_ee (completeness relation, requires Hilbert space context)
        @test isequal(simplify(σgg, h), simplify(1 - σee, h))
    end

    @testset "Algebraic relations" begin
        h = NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2)
        σee = Transition(h, :σ, 2, 2)
        σgg = Transition(h, :σ, 1, 1)

        # normal_order applies transition product rules
        @test isequal(simplify(normal_order(σ' * σ)), simplify(σee))
        # σ*σ' = σgg (ground state projector), which is 1 - σee
        no_result = simplify(normal_order(σ * σ'))
        @test isequal(no_result, simplify(σgg))
    end

    @testset "Product space operations" begin
        ha1 = NLevelSpace(:atom1, 2, 1)
        ha2 = NLevelSpace(:atom2, 2, 1)
        hprod = ha1 ⊗ ha2

        σ1 = Transition(hprod, :σ1, 1, 2, 1)
        σ2 = Transition(hprod, :σ2, 1, 2, 2)

        @test isequal(
            simplify(normal_order(σ1' * σ1)), simplify(Transition(hprod, :σ1, 2, 2, 1))
        )
        no_result = simplify(normal_order(σ2 * σ2'))
        @test isequal(no_result, simplify(Transition(hprod, :σ2, 1, 1, 2)))
        # Different subspaces don't interact
        @test isequal(simplify(σ1 * σ2), simplify(σ1 * σ2))
    end

    @testset "Symbolic levels" begin
        levels = (:g, :e, :a)
        h = NLevelSpace(:atom, levels)
        @test h.n == 3
        @test h.levels == [:g, :e, :a]
        @test h.ground_state == 1

        # Transition with symbol levels resolves to integer indices
        σge = Transition(h, :σ, :g, :e)
        @test σge.i == 1
        @test σge.j == 2
        σea = Transition(h, :σ, :e, :a)
        @test σea.i == 2
        @test σea.j == 3

        # Unknown level throws
        @test_throws ArgumentError Transition(h, :σ, :x, :g)

        # Integer construction still works
        σ12 = Transition(h, :σ, 1, 2)
        @test isequal(σ12, σge)

        # ProductSpace with symbolic levels
        hf = FockSpace(:c)
        hp = hf ⊗ h
        σge_p = Transition(hp, :σ, :g, :e, 2)
        @test σge_p.i == 1
        @test σge_p.j == 2
        @test σge_p.space_index == 2

        # Equality: spaces with different levels are not equal
        h_int = NLevelSpace(:atom, 3, 1)
        @test h != h_int
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
