using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify
using Test

@testset "nlevel" begin
    @testset "Transition construction — product space" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2, 2)
        @test acts_on(σ) == [2]
        @test_throws ArgumentError Transition(h, :σ, 1, 2, 1)
    end

    @testset "Adjoint" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = σ12'
        @test is_transition(σ21)
        @test σ21 == Transition(h, :σ, 2, 1)
        @test σ12'' == σ12
    end

    @testset "Arithmetic" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)

        @test iszero(σ12 * σ21 - Transition(h, :σ, 1, 1))
        @test iszero((σ12 + σ21) - (σ21 + σ12))
    end

    @testset "@qnumbers" begin
        h = NLevelSpace(:atom, 2, 1)
        @qnumbers σ::Transition(h, 1, 2)
        @test is_transition(σ)
        @test operator_name(σ) == :σ
    end

    @testset "Ground state projector" begin
        h = NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2)
        σee = Transition(h, :σ, 2, 2)
        σgg = Transition(h, :σ, 1, 1)

        # σ*σ' produces a ground-state projector under composition; expand
        # completeness explicitly: σgg = 1 - σee.
        @test isequal(simplify(expand_completeness(σ * σ')), simplify(1 - σee))
        # The h-aware overload applies completeness explicitly to user-constructed σgg.
        @test isequal(expand_completeness(σgg), simplify(1 - σee))
    end

    @testset "Ground state projector — non-default state" begin
        h = NLevelSpace(:atom, 4, 2)
        σ(i, j) = Transition(h, :σ, i, j)

        # The product rules remain local to the transition basis: σ¹²σ²¹ = σ¹¹,
        # while σ²¹σ¹² is the ground-state projector and expands on request.
        @test isequal(σ(1, 2) * σ(2, 1), simplify(σ(1, 1)))
        @test isequal(
            expand_completeness(σ(2, 1) * σ(1, 2)),
            simplify(1 - σ(1, 1) - σ(3, 3) - σ(4, 4)),
        )
    end

    @testset "Ground state projector — indexed" begin
        h = NLevelSpace(:atom, 3, 2)
        i = Index(h, :i, 5, h)
        σ(α, β) = IndexedOperator(Transition(h, :σ, α, β), i)

        # Completeness expansion preserves the site index on the surviving
        # excited-state projectors.
        @test isequal(
            expand_completeness(σ(2, 1) * adjoint(σ(2, 1))),
            simplify(1 - σ(1, 1) - σ(3, 3)),
        )
    end

    @testset "Completeness is explicit and site-local" begin
        h = NLevelSpace(:atom, 3, 1)
        σ(i, j) = Transition(h, :σ, i, j)

        cycle = σ(1, 2) * σ(2, 3) * σ(3, 1)
        @test iszero(cycle - σ(1, 1))
        @test isequal(
            expand_completeness(cycle),
            simplify(1 - σ(2, 2) - σ(3, 3)),
        )

        ha = NLevelSpace(:atom_a, 2, 1)
        hb = NLevelSpace(:atom_b, 3, 2)
        hproduct = ha ⊗ hb
        σa(i, j) = Transition(hproduct, :σa, i, j, 1)
        σb(i, j) = Transition(hproduct, :σb, i, j, 2)

        @test isequal(
            expand_completeness(σa(1, 2) * σa(2, 1)),
            simplify(1 - σa(2, 2)),
        )
        @test isequal(
            expand_completeness(σb(2, 1) * σb(1, 2)),
            simplify(1 - σb(1, 1) - σb(3, 3)),
        )
        @test iszero(commutator(σa(1, 2), σb(2, 1)))
    end

    @testset "Algebraic relations" begin
        h = NLevelSpace(:atom, 2, 1)
        σ = Transition(h, :σ, 1, 2)
        σee = Transition(h, :σ, 2, 2)

        # normal_order applies transition product rules
        @test isequal(simplify(normal_order(σ' * σ)), simplify(σee))
        # σ*σ' produces σgg via composition; expand explicitly: σgg = 1 - σee.
        no_result = simplify(expand_completeness(σ * σ'))
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
        # σ2 * σ2' produces σ2_11; expand explicitly to 1 - σ2_22.
        no_result = simplify(expand_completeness(σ2 * σ2'))
        @test isequal(no_result, simplify(1 - Transition(hprod, :σ2, 2, 2, 2)))
        # Different subspaces don't interact
        @test iszero(commutator(σ1, σ2))
        @test iszero(simplify(σ1 * σ2 - σ2 * σ1))
    end

    @testset "Symbolic levels" begin
        levels = (:g, :e, :a)
        h = NLevelSpace(:atom, levels)
        @test repr(h) == "ℋ(atom)"

        # Transition with symbol levels resolves to integer indices
        σge = Transition(h, :σ, :g, :e)
        @test σge == Transition(h, :σ, 1, 2)
        σea = Transition(h, :σ, :e, :a)
        @test σea == Transition(h, :σ, 2, 3)

        # Unknown level throws
        @test_throws ArgumentError Transition(h, :σ, :x, :g)

        # Integer construction still works
        σ12 = Transition(h, :σ, 1, 2)
        @test isequal(σ12, σge)

        # ProductSpace with symbolic levels
        hf = FockSpace(:c)
        hp = hf ⊗ h
        σge_p = Transition(hp, :σ, :g, :e, 2)
        @test acts_on(σge_p) == [2]
        @test σge_p == Transition(hp, :σ, 1, 2, 2)

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
        @inferred σ21 * σ12
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
