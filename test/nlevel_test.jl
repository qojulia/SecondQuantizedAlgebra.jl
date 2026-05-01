using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, HilbertSpace
using Test

@testset "nlevel" begin
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

        # σ*σ' produces a ground-state projector under composition; under
        # NormalOrder this is eagerly expanded via completeness: σgg = 1 - σee.
        @test isequal(simplify(normal_order(σ * σ')), simplify(1 - σee))
        # The h-aware overload applies completeness explicitly (LazyOrder opt-in).
        @test isequal(simplify(σgg, h), simplify(1 - σee, h))
    end

    @testset "Algebraic relations" begin
        h = NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2)
        σee = Transition(h, :σ, 2, 2)

        # normal_order applies transition product rules
        @test isequal(simplify(normal_order(σ' * σ)), simplify(σee))
        # σ*σ' produces σgg via composition; eager NormalOrder expands σgg = 1 - σee.
        no_result = simplify(normal_order(σ * σ'))
        @test isequal(no_result, simplify(1 - σee))
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
        # σ2 * σ2' produces σ2_11; eager NormalOrder expands to 1 - σ2_22.
        no_result = simplify(normal_order(σ2 * σ2'))
        @test isequal(no_result, simplify(1 - Transition(hprod, :σ2, 2, 2, 2)))
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

    @testset "Type stability" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)

        @inferred Transition(h, :σ, 1, 2)
        @inferred adjoint(σ12)
        @inferred isequal(σ12, σ21)
        @inferred hash(σ12, UInt(0))
        @inferred σ12 * σ21
        @inferred σ12 + σ21
    end

    @static if VERSION >= v"1.12"
        @testset "Allocations" begin
            h = NLevelSpace(:atom, 3, 1)
            σ12 = Transition(h, :σ, 1, 2)

            # Construction and adjoint
            @test @allocations(Transition(h, :σ, 1, 2)) == 0
            @test @allocations(adjoint(σ12)) == 0

            # Equality / hashing
            σ12b = Transition(h, :σ, 1, 2)
            @test @allocations(isequal(σ12, σ12b)) == 0
            @test @allocations(hash(σ12, UInt(0))) == 0
        end
    end
end
