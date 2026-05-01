using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CollectiveTransition, NO_INDEX

@testset "CollectiveTransition" begin
    h = NLevelSpace(:atom, 3, 1)
    S(i, j) = CollectiveTransition(h, :S, i, j)

    # ============================================================================
    # Field-level invariants
    # ============================================================================
    @testset "Struct fields and basic properties" begin
        op = S(1, 2)
        @test op.name == :S
        @test op.i == 1
        @test op.j == 2
        @test op.space_index == 1
        @test op.index == NO_INDEX        # never a user-supplied index

        # Adjoint flips i and j
        @test S(1, 2)' == S(2, 1)
        @test S(2, 2)' == S(2, 2)
        @test S(1, 3)' == S(3, 1)

        # Equality discriminates on every field
        @test S(1, 2) == S(1, 2)
        @test S(1, 2) != S(2, 1)
        @test S(1, 2) != S(1, 3)

        # Different name → different operator
        @test CollectiveTransition(h, :S, 1, 2) != CollectiveTransition(h, :T, 1, 2)
    end

    @testset "Bounds checking" begin
        @test_throws ArgumentError CollectiveTransition(h, :S, 0, 2)
        @test_throws ArgumentError CollectiveTransition(h, :S, 1, 4)
        @test_throws ArgumentError CollectiveTransition(h, :S, -1, 2)
    end

    @testset "Symbolic-level constructor" begin
        h_sym = NLevelSpace(:atom, (:g, :e, :a))
        Sge = CollectiveTransition(h_sym, :S, :g, :e)
        @test Sge.i == 1
        @test Sge.j == 2
    end

    @testset "ProductSpace constructors" begin
        hf = FockSpace(:cavity)
        hp = hf ⊗ h
        S_p(i, j) = CollectiveTransition(hp, :S, i, j, 2)
        @test S_p(1, 2).space_index == 2
        @test S_p(1, 2).i == 1
        @test S_p(1, 2).j == 2

        # Auto-detect (only one NLevelSpace in hp)
        @test CollectiveTransition(hp, :S, 1, 2) == S_p(1, 2)

        # Symbolic levels
        h_sym = NLevelSpace(:atom, (:g, :e))
        hp_sym = hf ⊗ h_sym
        Sge = CollectiveTransition(hp_sym, :S, :g, :e, 2)
        @test Sge.i == 1
        @test Sge.j == 2
    end

    @testset "Indexing is rejected" begin
        i = Index(h, :i, 10, h)
        @test_throws ArgumentError IndexedOperator(S(1, 2), i)
    end

    # ============================================================================
    # Algebra under NormalOrder (default)
    # ============================================================================
    @testset "Ordered products stay unchanged (single dict entry)" begin
        # Descending lex on (i, j): largest left
        for (a, b) in [
                (S(2, 1), S(1, 3)),    # a.i=2 > b.i=1
                (S(2, 3), S(2, 1)),    # a.i==b.i, a.j=3 > b.j=1
                (S(2, 2), S(2, 2)),    # equal
                (S(3, 1), S(1, 2)),    # a.i=3 > b.i=1
            ]
            r = a * b
            @test r isa QAdd
            @test length(r) == 1
            ops, c = only(collect(r))
            @test ops == [a, b]
            @test c == 1
        end
    end

    @testset "Unordered products fire the commutator (su(N) algebra)" begin
        # S(1,2) * S(2,3): a.j == b.i = 2 → +S(1,3); b.j == a.i?? 3 != 1, no second term.
        # Result: S(1,3) + S(2,3) * S(1,2)
        @test isequal(S(1, 2) * S(2, 3), S(1, 3) + S(2, 3) * S(1, 2))

        # S(1,2) * S(1,3): a.j != b.i, b.j != a.i → just a swap.
        # Result: S(1,3) * S(1,2)
        @test isequal(S(1, 2) * S(1, 3), S(1, 3) * S(1, 2))

        # S(1,2) * S(3,1): a.j != b.i (2 != 3); b.j == a.i = 1 → -S(3,2).
        # Result: -S(3,2) + S(3,1) * S(1,2)
        @test isequal(S(1, 2) * S(3, 1), -S(3, 2) + S(3, 1) * S(1, 2))

        # Both delta terms fire simultaneously: S(1,2) * S(2,1).
        # a.j == b.i = 2 → +S(a.i, b.j) = +S(1,1).
        # b.j == a.i = 1 → -S(b.i, a.j) = -S(2,2).
        # Result: S(1,1) - S(2,2) + S(2,1) * S(1,2). This is exactly the su(2)
        # commutator [S(1,2), S(2,1)] = S(1,1) - S(2,2) plus the swapped pair.
        @test isequal(S(1, 2) * S(2, 1), S(1, 1) - S(2, 2) + S(2, 1) * S(1, 2))

        # S(1,2) * S(1,1): _isordered((1,2),(1,1))? a.i==b.i, a.j=2 >= b.j=1 → ordered.
        # No swap fires; the product stays as [S(1,2), S(1,1)].
        r = S(1, 2) * S(1, 1)
        @test length(r) == 1
        ops, _ = only(collect(r))
        @test ops == [S(1, 2), S(1, 1)]

        # S(1,1) * S(1,2): unordered. a.j == b.i = 1 → +S(1,2). b.j != a.i → no second.
        # Result: S(1,2) + S(1,2) * S(1,1)
        @test isequal(S(1, 1) * S(1, 2), S(1, 2) + S(1, 2) * S(1, 1))
    end

    @testset "Self-product is non-zero (kept as 2-op product)" begin
        # S(1,2) * S(1,2): _isordered((1,2),(1,2))? equal → ordered. Stays as a
        # 2-op product. Note: collectively S^{12}·S^{12} = Σ_{k≠l} σ_k^{12} σ_l^{12}
        # is generally non-zero (only the k=l diagonal vanishes via σ²=0), so
        # leaving it un-collapsed is correct.
        r = S(1, 2) * S(1, 2)
        @test length(r) == 1
        ops, _ = only(collect(r))
        @test ops == [S(1, 2), S(1, 2)]
    end

    @testset "Adjoint involution" begin
        # `adjoint` reverses operator order and daggers each but does *not* re-run
        # the eager swap rule, so multi-op adjoints may not be in `_isordered_ct`
        # canonical form (apply `normal_order`/`simplify` if you need to compare).
        # The involution property holds independently:
        @test S(1, 2)'' == S(1, 2)
        @test S(2, 1)'' == S(2, 1)
        @test isequal((S(1, 2) * S(2, 3))'', S(1, 2) * S(2, 3))
        @test isequal((S(1, 2) * S(2, 1))'', S(1, 2) * S(2, 1))
    end

    @testset "Different subspaces do not interact" begin
        h2 = NLevelSpace(:atom2, 3, 1)
        hp = h ⊗ h2
        S1(i, j) = CollectiveTransition(hp, :S1, i, j, 1)
        S2(i, j) = CollectiveTransition(hp, :S2, i, j, 2)

        # Cross-subspace ops commute (different sites, site-sort handles order)
        @test isequal(S2(2, 1) * S1(1, 3), S1(1, 3) * S2(2, 1))
    end

    @testset "Different names on same subspace do not interact" begin
        # Two CT operators on the same NLevelSpace but with different names.
        # The swap rule requires a.name == b.name; otherwise nothing fires.
        T(i, j) = CollectiveTransition(h, :T, i, j)
        r = S(1, 2) * T(2, 1)
        @test length(r) == 1
        ops, _ = only(collect(r))
        @test length(ops) == 2
    end

    @testset "Mixing with Fock operators (Tavis-Cummings-like)" begin
        hc = FockSpace(:cavity)
        hp = hc ⊗ h
        @qnumbers a::Destroy(hp, 1)
        Sp(i, j) = CollectiveTransition(hp, :S, i, j, 2)

        # Commuting transitions: S(1,2) and S(2,3) → S(1,3) on the diagonal.
        # Multiplied by 'a' from the right (Fock commutes through atoms):
        @test isequal(Sp(1, 2) * Sp(2, 3) * a, a * Sp(1, 3) + Sp(2, 3) * Sp(1, 2) * a)

        # Two transitions that yield zero commutator: [S(1,2), S(1,3)] = 0 by su(N) rule.
        @test isequal(Sp(1, 2) * Sp(1, 3) * a, Sp(1, 3) * Sp(1, 2) * a)
        @test iszero(simplify(Sp(1, 2) * Sp(1, 3) * a - a * Sp(1, 3) * Sp(1, 2)))
    end

    @testset "CollectiveTransition and Transition are unrelated operators" begin
        # Both can sit on the same NLevelSpace, but they are NOT the same physics.
        # CT's index is NO_INDEX; Transition's user-supplied index differs.
        # _same_site requires matching index, so Transition·CT never fires reductions.
        σ(i, j) = Transition(h, :σ, i, j)
        @test S(1, 2) != σ(1, 2)
        @test hash(S(1, 2)) != hash(σ(1, 2))
        # Product: stored as 2-op product (no algebra fires)
        r = S(1, 2) * σ(2, 1)
        @test length(r) == 1
        ops, _ = only(collect(r))
        @test length(ops) == 2
    end

    # ============================================================================
    # Numeric conversion through ManyBodyBasis (symmetric subspace)
    # ============================================================================
    @testset "Numeric conversion: ManyBodyBasis" begin
        n_levels = 4
        n_atoms = 3
        h4 = NLevelSpace(:atoms, n_levels, 1)
        S4(i, j) = CollectiveTransition(h4, :S, i, j)

        b1 = NLevelBasis(n_levels)
        b = ManyBodyBasis(b1, bosonstates(b1, n_atoms))

        # Symmetric subspace dim = binomial(n_atoms + n_levels - 1, n_levels - 1)
        @test length(b) == binomial(n_atoms + n_levels - 1, n_levels - 1)

        for i in 1:n_levels, j in 1:n_levels
            op = S4(i, j)
            num = to_numeric(op, b)
            expected = QuantumOpticsBase.manybodyoperator(b, QuantumOpticsBase.transition(b1, i, j))
            @test num == expected
        end

        # Sum of population operators = N · I (symmetric subspace identity)
        pop_sum = sum(to_numeric(S4(k, k), b) for k in 1:n_levels)
        @test pop_sum ≈ Float64(n_atoms) * one(b)

        # Hermitian conjugation: S^{ij}† == S^{ji}
        @test to_numeric(S4(1, 2)', b) ≈ dagger(to_numeric(S4(1, 2), b))

        # Algebra: S^{12} S^{21} - S^{21} S^{12} = S^{11} - S^{22}
        lhs = to_numeric(S4(1, 2), b) * to_numeric(S4(2, 1), b) -
            to_numeric(S4(2, 1), b) * to_numeric(S4(1, 2), b)
        rhs = to_numeric(S4(1, 1), b) - to_numeric(S4(2, 2), b)
        @test lhs ≈ rhs
    end

    @testset "Numeric conversion: rejects non-NLevel one-body basis" begin
        n_atoms = 2
        b1 = FockBasis(3)
        b = ManyBodyBasis(b1, bosonstates(b1, n_atoms))
        Sx = CollectiveTransition(NLevelSpace(:atom, 3, 1), :S, 1, 2)
        @test_throws ArgumentError to_numeric(Sx, b)
    end

    # ============================================================================
    # Type stability
    # ============================================================================
    @testset "Type stability" begin
        @inferred CollectiveTransition(h, :S, 1, 2)
        @inferred adjoint(S(1, 2))
        @inferred isequal(S(1, 2), S(2, 1))
        @inferred hash(S(1, 2), UInt(0))
    end
end
