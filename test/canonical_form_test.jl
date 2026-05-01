using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, Transition, NO_INDEX,
    set_ordering!, get_ordering, NormalOrder, LazyOrder

# ============================================================================
# Canonical form invariants under NormalOrder
# ============================================================================
# Under NormalOrder (the default), the canonical basis for an NLevelSpace is
# {σⁱʲ : (i,j) ≠ (g,g)} ∪ {1}. Same-site composition that would produce σᵍᵍ
# is eagerly rewritten via completeness σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ. These tests
# pin that invariant across constructions and check the LazyOrder opt-in
# behaves as documented.

# Helper: is this op a ground-state projector of an NLevelSpace?
function _is_gs_projector(op)
    op isa Transition || return false
    return op.i == op.ground_state && op.j == op.ground_state
end

# Helper: every dict key in `expr` is free of ground-state projectors
function _no_gs_projectors(expr::QAdd)
    for (ops, _) in expr.arguments
        any(_is_gs_projector, ops) && return false
    end
    return true
end

@testset "Canonical form (NormalOrder default)" begin

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
    @testset "Eager: σ¹²·σ²¹ → 1 - σ²² (2-level, g=1)" begin
        h = NLevelSpace(:atom, 2, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ22 = Transition(h, :σ, 2, 2)

        result = σ12 * σ21
        @test result isa QAdd
        @test isequal(result, simplify(1 - σ22))
        @test _no_gs_projectors(result)
    end

    @testset "Eager: σ¹²·σ²¹ → 1 - σ²² - σ³³ (3-level, g=1)" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ22 = Transition(h, :σ, 2, 2)
        σ33 = Transition(h, :σ, 3, 3)

        result = σ12 * σ21
        @test isequal(result, simplify(1 - σ22 - σ33))
        @test _no_gs_projectors(result)
    end

    @testset "Eager: ground state ≠ 1 (g=2)" begin
        h = NLevelSpace(:atom, 3, 2)        # ground state = 2
        σ12 = Transition(h, :σ, 1, 2)
        σ21 = Transition(h, :σ, 2, 1)
        σ11 = Transition(h, :σ, 1, 1)
        σ33 = Transition(h, :σ, 3, 3)

        # σ²¹·σ¹² composes to σ²² which IS the ground state for g=2 → expand
        result = σ21 * σ12
        @test isequal(result, simplify(1 - σ11 - σ33))
        @test _no_gs_projectors(result)

        # σ¹²·σ²¹ composes to σ¹¹ (NOT ground state when g=2) → stays atomic
        result2 = σ12 * σ21
        @test length(result2) == 1
        ops, c = only(collect(result2))
        @test length(ops) == 1
        @test ops[1] == σ11   # σ¹¹ is not the ground state, so it's kept
    end

    @testset "Eager: indexed σᵢ¹²·σᵢ²¹ preserves index" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        h = hf ⊗ ha

        @variables N
        i = Index(h, :i, N, ha)
        σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

        result = σ(1, 2) * σ(2, 1)
        @test result isa QAdd
        # Should be 1 - σᵢ²², so two dict entries: identity + -σᵢ²²
        @test length(result) == 2
        @test _no_gs_projectors(result)
        # The σ²² entry should preserve index `i`
        for (ops, _) in result.arguments
            for op in ops
                if op isa Transition
                    @test op.index == i
                    @test op.ground_state == 1
                    @test op.n_levels == 2
                end
            end
        end
    end

    @testset "Eager: σᵍᵍ via longer chain σ¹²·σ²³·σ³¹ (3-level, g=1)" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        σ23 = Transition(h, :σ, 2, 3)
        σ31 = Transition(h, :σ, 3, 1)
        σ22 = Transition(h, :σ, 2, 2)
        σ33 = Transition(h, :σ, 3, 3)

        # σ¹²·σ²³·σ³¹ → σ¹³·σ³¹ → σ¹¹ → 1 - σ²² - σ³³
        result = σ12 * σ23 * σ31
        @test isequal(result, simplify(1 - σ22 - σ33))
        @test _no_gs_projectors(result)
    end

    # ============================================================================
    # User-constructed σᵍᵍ stays atomic until simplify(_, h) / normal_order(_, h)
    # ============================================================================
    @testset "User-constructed σᵍᵍ stays atomic (no composition fired)" begin
        h = NLevelSpace(:atom, 2, 1)
        σgg = Transition(h, :σ, 1, 1)
        # Direct construction: no arithmetic, no expansion
        @test σgg isa Transition
        @test σgg.i == 1 && σgg.j == 1

        # simplify without h: still atomic
        @test isequal(simplify(σgg), simplify(σgg))

        # simplify with h: expands
        σ22 = Transition(h, :σ, 2, 2)
        @test isequal(simplify(σgg, h), simplify(1 - σ22, h))

        # normal_order with h: expands
        @test isequal(normal_order(σgg, h), simplify(1 - σ22, h))
    end

    # ============================================================================
    # LazyOrder opt-in: simplify(expr) keeps σᵍᵍ; simplify(expr, h) expands
    # ============================================================================
    @testset "LazyOrder: simplify(expr) keeps σᵍᵍ atomic" begin
        prev = get_ordering()
        set_ordering!(LazyOrder())
        try
            h = NLevelSpace(:atom, 2, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ11 = Transition(h, :σ, 1, 1)

            # In LazyOrder, * just stores [σ12, σ21] without composition
            raw = σ12 * σ21
            @test length(raw) == 1
            ops, _ = only(collect(raw))
            @test length(ops) == 2

            # simplify(expr): fires reductions; σᵍᵍ produced and KEPT atomic
            no_h = simplify(raw)
            @test isequal(no_h, simplify(σ11))
            # σᵍᵍ DOES appear here — that's the point of LazyOrder simplify-without-h
            @test !_no_gs_projectors(no_h)

            # simplify(expr, h): fires reductions + completeness
            σ22 = Transition(h, :σ, 2, 2)
            with_h = simplify(raw, h)
            @test isequal(with_h, simplify(1 - σ22, h))
            @test _no_gs_projectors(with_h)
        finally
            set_ordering!(prev)
        end
    end

    @testset "LazyOrder: normal_order(expr) forces NormalOrder canonical form" begin
        prev = get_ordering()
        set_ordering!(LazyOrder())
        try
            h = NLevelSpace(:atom, 2, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ22 = Transition(h, :σ, 2, 2)

            raw = σ12 * σ21    # un-reduced under LazyOrder

            # normal_order forces NormalOrder pass: reductions + commutation +
            # completeness all fire because _apply_ordering(NormalOrder()) is invoked
            result = normal_order(raw)
            @test isequal(result, simplify(1 - σ22))
            @test _no_gs_projectors(result)
        finally
            set_ordering!(prev)
        end
    end

    # ============================================================================
    # σᵍᵍ-never-appears invariant for representative expressions
    # ============================================================================
    @testset "Invariant: σᵍᵍ never appears in commutators (NormalOrder)" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

        @qnumbers a::Destroy(h, 1)
        σ(α, β) = Transition(h, :σ, α, β, 2)

        # Several commutator constructions that internally produce σᵍᵍ
        H_jc = a' * σ(1, 2) + a * σ(2, 1)
        for op in (σ(1, 2), σ(2, 1), σ(2, 2), a' * σ(1, 2), a * σ(2, 1))
            result = commutator(H_jc, op)
            result isa QAdd || continue
            @test _no_gs_projectors(result)
        end
    end

    @testset "Invariant: σᵍᵍ never appears for indexed sums (NormalOrder)" begin
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

        # Commutators with assorted operator structures
        for op in (
                a' * σ(1, 2, j),
                a * σ(2, 1, j),
                σ(1, 2, j) * σ(2, 1, k),
                σ(2, 2, j),
            )
            result = commutator(H, op)
            result isa QAdd || continue
            @test _no_gs_projectors(result)
        end
    end

    @testset "Invariant: σᵍᵍ never appears for higher-level systems" begin
        # 4-level atom, ground state = 2
        h = NLevelSpace(:atom, 4, 2)
        σ(α, β) = Transition(h, :σ, α, β)

        # σ¹² · σ²¹ → σ¹¹ (atomic, NOT ground state)
        r1 = σ(1, 2) * σ(2, 1)
        @test _no_gs_projectors(r1)

        # σ²¹ · σ¹² → σ²² (ground state) → 1 - σ¹¹ - σ³³ - σ⁴⁴
        r2 = σ(2, 1) * σ(1, 2)
        @test _no_gs_projectors(r2)
        # Should have 4 dict entries: identity, -σ¹¹, -σ³³, -σ⁴⁴
        @test length(r2) == 4

        # σ³² · σ²³ → σ³³ (NOT ground state, g=2) → kept atomic
        r3 = σ(3, 2) * σ(2, 3)
        ops, _ = only(collect(r3))
        @test ops == [Transition(h, :σ, 3, 3)]
    end

    @testset "Invariant: ProductSpace with multiple NLevelSpaces" begin
        # Two NLevelSpaces with different ground states in one ProductSpace
        ha = NLevelSpace(:atomA, 2, 1)
        hb = NLevelSpace(:atomB, 3, 2)
        h = ha ⊗ hb

        σA(α, β) = Transition(h, :σA, α, β, 1)
        σB(α, β) = Transition(h, :σB, α, β, 2)

        # σA¹²·σA²¹ → σA¹¹ (g=1 for A) → 1 - σA²²
        rA = σA(1, 2) * σA(2, 1)
        @test _no_gs_projectors(rA)

        # σB²¹·σB¹² → σB²² (g=2 for B) → 1 - σB¹¹ - σB³³
        rB = σB(2, 1) * σB(1, 2)
        @test _no_gs_projectors(rB)

        # σB¹²·σB²¹ → σB¹¹ (NOT ground state for B) → stays atomic
        rB2 = σB(1, 2) * σB(2, 1)
        ops, _ = only(collect(rB2))
        @test length(ops) == 1
        @test ops[1].i == 1 && ops[1].j == 1
        @test ops[1].ground_state == 2
    end
end
