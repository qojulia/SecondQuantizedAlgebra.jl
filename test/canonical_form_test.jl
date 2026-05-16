using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, Transition, NO_INDEX

# ============================================================================
# Canonical form invariants
# ============================================================================
# The canonical basis for an NLevelSpace is {σⁱʲ : (i,j) ≠ (g,g)} ∪ {1}.
# Same-site composition that would produce σᵍᵍ is eagerly rewritten via
# completeness σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ. These tests pin that invariant across
# constructions.

# Helper: is this op a ground-state projector of an NLevelSpace?
function _is_gs_projector(op)
    op isa Transition || return false
    return op.i == op.ground_state && op.j == op.ground_state
end

# Helper: every dict key in `expr` is free of ground-state projectors
function _no_gs_projectors(expr::QAdd)
    for term in keys(expr.arguments)
        any(_is_gs_projector, term.ops) && return false
    end
    return true
end

@testset "Canonical form" begin

    # ============================================================================
    # Field-level invariants on Transition
    # ============================================================================
    @testset "Transition carries ground_state and n_levels" begin
        h = NLevelSpace(:atom, 3, 2)        # 3 levels, ground state = 2
        σ = Transition(h, :σ, 1, 3)
        @test σ.ground_state == 2
        @test σ.n_levels == 3

        # Defaults: ground_state == 1
        h2 = NLevelSpace(:atom, 2)
        σ2 = Transition(h2, :σ, 1, 2)
        @test σ2.ground_state == 1
        @test σ2.n_levels == 2

        # ProductSpace: GS comes from the relevant subspace
        hf = FockSpace(:c)
        hp = hf ⊗ h
        σp = Transition(hp, :σ, 1, 3, 2)
        @test σp.ground_state == 2
        @test σp.n_levels == 3
    end

    @testset "Field preservation through ops" begin
        h = NLevelSpace(:atom, 4, 3)        # ground state = 3
        σ = Transition(h, :σ, 1, 2)

        # adjoint preserves
        σadj = adjoint(σ)
        @test σadj.ground_state == 3
        @test σadj.n_levels == 4

        # IndexedOperator preserves
        i = Index(h, :i, 10, h)
        σi = IndexedOperator(σ, i)
        @test σi.ground_state == 3
        @test σi.n_levels == 4

        # change_index preserves
        j = Index(h, :j, 10, h)
        σj = change_index(σi, i, j)
        @test σj.ground_state == 3
        @test σj.n_levels == 4
    end

    @testset "Equality: different ground_state distinguishes operators" begin
        h_g1 = NLevelSpace(:atom, 3, 1)
        h_g2 = NLevelSpace(:atom, 3, 2)
        σ_g1 = Transition(h_g1, :σ, 1, 2)
        σ_g2 = Transition(h_g2, :σ, 1, 2)
        @test σ_g1 != σ_g2
        @test hash(σ_g1) != hash(σ_g2)
    end

    @testset "Equality: different n_levels distinguishes operators" begin
        h_n2 = NLevelSpace(:atom, 2, 1)
        h_n3 = NLevelSpace(:atom, 3, 1)
        σ_n2 = Transition(h_n2, :σ, 1, 2)
        σ_n3 = Transition(h_n3, :σ, 1, 2)
        @test σ_n2 != σ_n3
        @test hash(σ_n2) != hash(σ_n3)
    end

    # ============================================================================
    # Eager completeness: same-site composition producing σᵍᵍ
    # ============================================================================
    @testset "Composition+expand: σ¹²·σ²¹ → 1 - σ²² (2-level, g=1)" begin
        h = NLevelSpace(:atom, 2, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ22 = Transition(h, :σ, 2, 2)

        result = expand_completeness(σ12 * σ21)
        @test result isa QAdd
        @test isequal(result, simplify(1 - σ22))
        @test _no_gs_projectors(result)
    end

    @testset "Composition+expand: σ¹²·σ²¹ → 1 - σ²² - σ³³ (3-level, g=1)" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ22 = Transition(h, :σ, 2, 2)
        σ33 = Transition(h, :σ, 3, 3)

        result = expand_completeness(σ12 * σ21)
        @test isequal(result, simplify(1 - σ22 - σ33))
        @test _no_gs_projectors(result)
    end

    @testset "Composition+expand: ground state ≠ 1 (g=2)" begin
        h = NLevelSpace(:atom, 3, 2)        # ground state = 2
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ11 = Transition(h, :σ, 1, 1)
        σ33 = Transition(h, :σ, 3, 3)

        # σ²¹·σ¹² composes to σ²² which IS the ground state for g=2 → expand
        result = expand_completeness(σ21 * σ12)
        @test isequal(result, simplify(1 - σ11 - σ33))
        @test _no_gs_projectors(result)

        # σ¹²·σ²¹ composes to σ¹¹ (NOT ground state when g=2) → stays atomic
        result2 = σ12 * σ21
        @test length(result2) == 1
        (term, c) = only(collect(result2)); ops = term.ops
        @test length(ops) == 1
        @test ops[1] == σ11   # σ¹¹ is not the ground state, so it's kept
    end

    @testset "Composition+expand: indexed σᵢ¹²·σᵢ²¹ preserves index" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        h = hf ⊗ ha

        @variables N
        i = Index(h, :i, N, ha)
        σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

        result = expand_completeness(σ(1, 2) * σ(2, 1))
        @test result isa QAdd
        # Should be 1 - σᵢ²², so two dict entries: identity + -σᵢ²²
        @test length(result) == 2
        @test _no_gs_projectors(result)
        # The σ²² entry should preserve index `i`
        for term in keys(result.arguments)
            for op in term.ops
                if op isa Transition
                    @test op.index == i
                    @test op.ground_state == 1
                    @test op.n_levels == 2
                end
            end
        end
    end

    @testset "Composition+expand: σᵍᵍ via longer chain σ¹²·σ²³·σ³¹ (3-level, g=1)" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ23 = Transition(h, :σ, 2, 3)
        σ31 = Transition(h, :σ, 3, 1)
        σ22 = Transition(h, :σ, 2, 2)
        σ33 = Transition(h, :σ, 3, 3)

        # σ¹²·σ²³·σ³¹ → σ¹³·σ³¹ → σ¹¹ → 1 - σ²² - σ³³ (after expand_completeness)
        result = expand_completeness(σ12 * σ23 * σ31)
        @test isequal(result, simplify(1 - σ22 - σ33))
        @test _no_gs_projectors(result)
    end

    # ============================================================================
    # Eager * canonicalizes any product containing σᵍᵍ
    # ============================================================================
    @testset "Different-site σᵍᵍ_i · σᵍᵍ_j (atomic; expand via expand_completeness)" begin
        h = NLevelSpace(:atom, 3, 2)
        i = SecondQuantizedAlgebra.Index(h, :i, 10, 1)
        j = SecondQuantizedAlgebra.Index(h, :j, 10, 1)
        σi = SecondQuantizedAlgebra.IndexedOperator(Transition(h, :σ, 2, 2), i)
        σj = SecondQuantizedAlgebra.IndexedOperator(Transition(h, :σ, 2, 2), j)
        expanded = expand_completeness(σi * σj)
        @test _no_gs_projectors(expanded)
        @test isequal(expanded, expand_completeness(normal_order(σi * σj)))
    end

    @testset "σᵍᵍ * σⁱʲ (expand opt-in)" begin
        h = NLevelSpace(:atom, 3, 1)
        σgg = Transition(h, :σ, 1, 1)
        σ12 = Transition(h, :σ, 1, 2)
        @test isequal(expand_completeness(σgg * σ12), expand_completeness(normal_order(σgg * σ12)))
    end

    # ============================================================================
    # User-constructed σᵍᵍ stays atomic until expand_completeness fires
    # ============================================================================
    @testset "User-constructed σᵍᵍ stays atomic (no composition fired)" begin
        h = NLevelSpace(:atom, 2, 1)
        σgg = Transition(h, :σ, 1, 1)
        # Direct construction: no arithmetic, no expansion
        @test σgg isa Transition
        @test σgg.i == 1 && σgg.j == 1

        # simplify without h: still atomic
        @test isequal(simplify(σgg), simplify(σgg))

        # simplify wrapped in expand_completeness: expands
        σ22 = Transition(h, :σ, 2, 2)
        @test isequal(expand_completeness(simplify(σgg)), expand_completeness(simplify(1 - σ22)))

        # normal_order wrapped in expand_completeness: expands
        @test isequal(expand_completeness(normal_order(σgg)), expand_completeness(simplify(1 - σ22)))
    end

    # ============================================================================
    # σᵍᵍ-never-appears invariant for representative expressions
    # ============================================================================
    @testset "After expand_completeness: σᵍᵍ removed from commutators" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

        @qnumbers a::Destroy(h, 1)
        σ(α, β) = Transition(h, :σ, α, β, 2)

        H_jc = a' * σ(1, 2) + a * σ(2, 1)
        for op in (σ(1, 2), σ(2, 1), σ(2, 2), a' * σ(1, 2), a * σ(2, 1))
            result = expand_completeness(commutator(H_jc, op))
            result isa QAdd || continue
            @test _no_gs_projectors(result)
        end
    end

    @testset "After expand_completeness: σᵍᵍ removed for indexed sums" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

        @variables N Δ
        @qnumbers a::Destroy(h, 1)
        σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
        g(idx) = IndexedVariable(:g, idx)

        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        k = Index(h, :k, N, ha)

        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

        for op in (
                a' * σ(1, 2, j),
                a * σ(2, 1, j),
                σ(1, 2, j) * σ(2, 1, k),
                σ(2, 2, j),
            )
            result = expand_completeness(commutator(H, op))
            result isa QAdd || continue
            @test _no_gs_projectors(result)
        end
    end

    @testset "After expand_completeness: σᵍᵍ removed for higher-level systems" begin
        h = NLevelSpace(:atom, 4, 2)
        σ(α, β) = Transition(h, :σ, α, β)

        # σ¹² · σ²¹ → σ¹¹ (atomic, NOT ground state)
        r1 = σ(1, 2) * σ(2, 1)
        @test _no_gs_projectors(r1)

        # σ²¹ · σ¹² → σ²² (ground state). After expand_completeness: 1 - σ¹¹ - σ³³ - σ⁴⁴
        r2 = expand_completeness(σ(2, 1) * σ(1, 2))
        @test _no_gs_projectors(r2)
        @test length(r2) == 4

        # σ³² · σ²³ → σ³³ (NOT ground state, g=2) → kept atomic
        r3 = σ(3, 2) * σ(2, 3)
        (term, _) = only(collect(r3))
        @test term.ops == [Transition(h, :σ, 3, 3)]
    end

    @testset "After expand_completeness: ProductSpace with multiple NLevelSpaces" begin
        ha = NLevelSpace(:atomA, 2, 1)
        hb = NLevelSpace(:atomB, 3, 2)
        h = ha ⊗ hb

        σA(α, β) = Transition(h, :σA, α, β, 1)
        σB(α, β) = Transition(h, :σB, α, β, 2)

        rA = expand_completeness(σA(1, 2) * σA(2, 1))
        @test _no_gs_projectors(rA)

        rB = expand_completeness(σB(2, 1) * σB(1, 2))
        @test _no_gs_projectors(rB)

        rB2 = σB(1, 2) * σB(2, 1)
        term, _ = only(collect(rB2))
        @test length(term.ops) == 1
        @test term.ops[1].i == 1 && term.ops[1].j == 1
        @test term.ops[1].ground_state == 2
    end
end
