using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, _to_cnum, NO_INDEX,
    create_index_arrays, operators, prefactor, sorted_arguments

@testset "Indexing" begin

    # ========== Setup ==========
    hf = FockSpace(:f)
    hn = NLevelSpace(:n, 2, 1)
    hp = PauliSpace(:p)
    hs = SpinSpace(:s)
    hq = PhaseSpace(:q)
    h_prod = hf ⊗ hn

    a = Destroy(hf, :a)
    ad = a'
    sigma = Transition(hn, :σ, 1, 2)

    # ========== Index construction ==========
    @testset "Index basics" begin
        i = Index(hf, :i, 10, hf)
        @test i.name == :i
        @test i.range == 10
        @test i.space_index == 1
        @test has_index(i)
        @test !has_index(NO_INDEX)

        # Index equality ignores sym (only name, range, space_index)
        j = Index(hf, :i, 10, hf)
        @test i == j

        # Different name
        k = Index(hf, :k, 10, hf)
        @test i != k

        # Different range
        m = Index(hf, :i, 5, hf)
        @test i != m
    end

    @testset "Index in ProductSpace" begin
        i = Index(h_prod, :i, 10, hf)
        @test i.space_index == 1

        j = Index(h_prod, :j, 10, hn)
        @test j.space_index == 2

        # Error for space not in ProductSpace
        hother = FockSpace(:other)
        @test_throws ArgumentError Index(h_prod, :i, 10, hother)
    end

    @testset "Index from integer space_index" begin
        i = Index(hf, :i, 10, 1)
        @test i.space_index == 1
        @test has_index(i)
    end

    @testset "Index hashing" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :i, 10, hf)
        @test hash(i) == hash(j)

        k = Index(hf, :k, 10, hf)
        @test hash(i) != hash(k)
    end

    # ========== IndexedOperator ==========
    @testset "IndexedOperator — all operator types" begin
        i = Index(hf, :i, 10, hf)

        # Fock
        ai = IndexedOperator(a, i)
        @test ai isa Destroy
        @test ai.index == i
        @test has_index(ai.index)

        adi = IndexedOperator(ad, i)
        @test adi isa Create
        @test adi.index == i

        # NLevel
        j = Index(hn, :j, 5, hn)
        sj = IndexedOperator(sigma, j)
        @test sj isa Transition
        @test sj.index == j

        # Pauli
        sx = Pauli(hp, :σ, 1)
        ip = Index(hp, :i, 10, hp)
        sxi = IndexedOperator(sx, ip)
        @test sxi isa Pauli
        @test sxi.index == ip

        # Spin
        Sx = Spin(hs, :S, 1)
        is_ = Index(hs, :i, 10, hs)
        Sxi = IndexedOperator(Sx, is_)
        @test Sxi isa Spin
        @test Sxi.index == is_

        # PhaseSpace
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        iq = Index(hq, :i, 10, hq)
        xi = IndexedOperator(x, iq)
        pi = IndexedOperator(p, iq)
        @test xi isa Position
        @test pi isa Momentum
        @test xi.index == iq
    end

    @testset "IndexedOperator preserves adjoint" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'
        @test adi isa Create
        @test adi.index == i
    end

    # ========== Same-site with indices ==========
    @testset "Indexed operators on different indices are different sites" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        # Different index → different site → commute
        @test commutator(ai, aj') isa QAdd
        result = commutator(ai, aj')
        @test iszero(result)
    end

    @testset "Indexed operators on same index interact" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # Same index → same site → [a_i, a†_i] = 1
        result = commutator(ai, adi)
        @test result isa QAdd
        @test length(result) == 1
        @test isempty(only(collect(result)).first)
        @test only(collect(result)).second == 1
    end

    @testset "Indexed vs non-indexed are different sites" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        # Non-indexed a and indexed a_i have different indices (NO_INDEX vs i)
        result = commutator(a, ai')
        @test iszero(result)
    end

    # ========== IndexedVariable ==========
    @testset "IndexedVariable" begin
        i = Index(hf, :i, 10, hf)
        gi = IndexedVariable(:g, i)
        @test gi isa Symbolics.Num

        # Different index gives different variable
        j = Index(hf, :j, 10, hf)
        gj = IndexedVariable(:g, j)
        @test !isequal(gi, gj)
    end

    @testset "DoubleIndexedVariable" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)

        Jij = DoubleIndexedVariable(:J, i, j)
        @test Jij isa Symbolics.Num

        # identical=false → 0 when i==j
        Jii_nodiag = DoubleIndexedVariable(:J, i, i; identical = false)
        @test isequal(Jii_nodiag, Symbolics.Num(0))

        # identical=true (default) → nonzero when i==i
        Jii = DoubleIndexedVariable(:J, i, i)
        @test !isequal(Jii, Symbolics.Num(0))
    end

    @testset "IndexedVariable in prefactor" begin
        i = Index(hf, :i, 10, hf)
        gi = IndexedVariable(:g, i)
        ai = IndexedOperator(a, i)

        # Can multiply: g_i * a_i
        expr = gi * ai
        @test expr isa QAdd
    end

    # ========== Σ (symbolic sums) ==========
    @testset "Sigma single index" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        s = Σ(ai, i)
        @test s isa QAdd
        @test length(s.indices) == 1
        @test s.indices[1] == i
        @test isempty(s.non_equal)
    end

    @testset "Sigma with product" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        s = Σ(adi * ai, i)
        @test s isa QAdd
        @test length(s.indices) == 1
        @test length(s) == 1
        @test length(operators(only(sorted_arguments(s)))) == 2
    end

    @testset "Sigma with non_equal" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)

        s = Σ(ai, i, [j])
        @test length(s.non_equal) == 1
        @test s.non_equal[1] == (i, j)
    end

    @testset "Sigma multi-index" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)

        # Σ(a_i, i, j) — j is independent → folds to 10 * Σ(a_i, i)
        s = Σ(ai, i, j)
        @test s isa QAdd
        @test length(s.indices) == 1
        @test i in s.indices
        @test isequal(s, Σ(ai, i) * 10)

        # Multi-index where both are needed (different spaces, both indexed)
        i_n = Index(h_prod, :i, 10, hn)
        m_f = Index(h_prod, :m, 5, hf)
        σi = IndexedOperator(sigma, i_n)
        am = IndexedOperator(a, m_f)
        s2 = Σ(σi * am, i_n, m_f)
        @test length(s2.indices) == 2
    end

    @testset "Sigma with scalar" begin
        i = Index(hf, :i, 10, hf)
        s = Σ(1, i)
        @test s isa QAdd
        # Scalar 1 doesn't depend on i → simplified to 10 * 1
        @test isempty(get_indices(s))
    end

    @testset "Sigma with QSym" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        s = Σ(ai, i)
        @test s isa QAdd
        @test length(s) == 1
    end

    @testset "Sigma alias ∑" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        @test ∑(ai, i) isa QAdd
    end

    # ========== change_index ==========
    @testset "change_index — operators" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)

        # Destroy
        ai = IndexedOperator(a, i)
        aj = change_index(ai, i, j)
        @test aj isa Destroy
        @test aj.index == j

        # Create
        adi = ai'
        adj = change_index(adi, i, j)
        @test adj isa Create
        @test adj.index == j

        # Transition
        k = Index(hn, :k, 5, hn)
        l = Index(hn, :l, 5, hn)
        sk = IndexedOperator(sigma, k)
        sl = change_index(sk, k, l)
        @test sl isa Transition
        @test sl.index == l

        # Pauli
        sx = Pauli(hp, :σ, 1)
        ip = Index(hp, :i, 10, hp)
        jp = Index(hp, :j, 10, hp)
        sxi = IndexedOperator(sx, ip)
        sxj = change_index(sxi, ip, jp)
        @test sxj isa Pauli
        @test sxj.index == jp

        # Spin
        Sx = Spin(hs, :S, 1)
        is_ = Index(hs, :i, 10, hs)
        js_ = Index(hs, :j, 10, hs)
        Sxi = IndexedOperator(Sx, is_)
        Sxj = change_index(Sxi, is_, js_)
        @test Sxj isa Spin
        @test Sxj.index == js_

        # Position / Momentum
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        iq = Index(hq, :i, 10, hq)
        jq = Index(hq, :j, 10, hq)
        xi = IndexedOperator(x, iq)
        pi = IndexedOperator(p, iq)
        xj = change_index(xi, iq, jq)
        pj = change_index(pi, iq, jq)
        @test xj.index == jq
        @test pj.index == jq
    end

    @testset "change_index — no-op when index doesn't match" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        k = Index(hf, :k, 10, hf)

        ai = IndexedOperator(a, i)
        result = change_index(ai, j, k)  # i != j, so no change
        @test result.index == i
    end

    @testset "change_index — QAdd (product)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        m = adi * ai  # a†_i * a_i
        mj = change_index(m, i, j)
        @test mj isa QAdd
        for op in operators(mj)
            @test op.index == j
        end
    end

    @testset "change_index — QAdd prefactor with IndexedVariable" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        gi = IndexedVariable(:g, i)
        ai = IndexedOperator(a, i)

        m = gi * ai  # g(i) * a_i
        mj = change_index(m, i, j)
        @test operators(mj)[1].index == j
        # Prefactor should also be substituted: g(i) → g(j)
        gj = IndexedVariable(:g, j)
        @test isequal(real(prefactor(mj)), gj)
    end

    @testset "change_index — QAdd" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        s = ai + adi  # a_i + a†_i
        sj = change_index(s, i, j)
        @test sj isa QAdd
        for (ops, _) in sj
            for op in ops
                @test op.index == j
            end
        end
    end

    @testset "change_index — QAdd with indices and non_equal" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        k = Index(hf, :k, 10, hf)
        ai = IndexedOperator(a, i)

        s = Σ(ai, i, [j])  # Σ_i (i≠j) a_i
        sk = change_index(s, j, k)
        @test sk isa QAdd
        # non_equal should be updated: (i, j) → (i, k)
        @test any(p -> p == (i, k), sk.non_equal)
    end

    @testset "change_index — Number passthrough" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        @test change_index(42, i, j) == 42
    end

    @testset "change_index — CNum" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        gi = IndexedVariable(:g, i)
        c = _to_cnum(gi)
        cj = change_index(c, i, j)
        gj = IndexedVariable(:g, j)
        @test isequal(real(cj), gj)
    end

    # ========== get_indices ==========
    @testset "get_indices" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)

        # Number → empty
        @test isempty(get_indices(42))

        # Non-indexed operator → empty
        @test isempty(get_indices(a))

        # Indexed operator → [i]
        ai = IndexedOperator(a, i)
        @test get_indices(ai) == [i]

        # QAdd with indexed operators
        aj = IndexedOperator(a, j)
        m = ai' * aj
        inds = get_indices(m)
        @test length(inds) == 2
        @test i in inds
        @test j in inds

        # QAdd with repeated index → unique
        m2 = ai' * ai
        @test length(get_indices(m2)) == 1

        # QAdd
        s = ai + aj
        inds = get_indices(s)
        @test i in inds
        @test j in inds

        # QAdd with summation indices
        s2 = Σ(ai, i)
        inds = get_indices(s2)
        @test i in inds
    end

    # ========== Eager diagonal splitting ==========
    @testset "Σ construction — diagonal split on product" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ_i(a†_j * a_i) — splits into: Σ_{i≠j}(a_i * a†_j) + a†_j * a_j + 1
        # Eager ordering rewrites a†_j * a_i → a_i * a†_j (off-diagonal),
        # and the diagonal (i→j) a_j * a†_j → a†_j * a_j + 1 (commutation).
        expr = Σ(adj * ai, i)
        @test expr isa QAdd
        # Should have 3 terms: off-diagonal product + diagonal normal-ordered + scalar 1
        @test length(expr) == 3
        # Should have i≠j constraint
        @test any(p -> p == (i, j) || p == (j, i), expr.non_equal)
        # Diagonal term: a†_j * a_j (eager ordering composed a_j * a†_j → a†_j * a_j + 1)
        aj = IndexedOperator(a, j)
        diag_key = QSym[adj, aj]
        @test haskey(expr.arguments, diag_key)
        # Scalar term from commutator
        @test haskey(expr.arguments, QSym[])
    end

    @testset "Σ construction — already non_equal → no split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ_{i≠j}(a†_j * a_i) — already constrained, no diagonal extraction
        expr = Σ(adj * ai, i, [j])
        @test length(expr) == 1
    end

    @testset "Σ construction — independent index folds to range" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)

        # Σ(Σ(a_i, i), j) — j is independent → 10 * Σ(a_i, i)
        inner = Σ(ai, i)
        outer = Σ(inner, j)
        expected = inner * 10
        @test isequal(outer, expected)
    end

    @testset "Σ construction — symbolic N independent index" begin
        @variables N_sym::Real
        i = Index(hf, :i, N_sym, hf)
        j = Index(hf, :j, N_sym, hf)
        ai = IndexedOperator(a, i)

        # Σ(-a_i, i, j) = N * Σ(-a_i, i) since j is independent
        result = Σ(-ai, i, j)
        expected = Σ(-ai, i) * N_sym
        @test isequal(simplify(result), simplify(expected))
    end

    @testset "Σ construction — nested same-space double sum" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ(Σ(a†_j * a_i, i), j) — nested sum, same space
        # Should split: Σ_{i≠j} Σ_j(a†_j a_i) + Σ_j(a†_j a_j)
        result = Σ(Σ(adj * ai, i), j)
        @test result isa QAdd
        @test length(result) > 1
    end

    @testset "Sum * QSym — diagonal split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ_i(a_i) * a†_j — splits diagonal
        # Eager ordering: a_i * a†_j (off-diag), a_j * a†_j → a†_j * a_j + 1 (diag)
        sum_expr = Σ(ai, i)
        result = sum_expr * adj
        @test result isa QAdd
        @test length(result) == 3
        @test any(p -> p == (i, j) || p == (j, i), result.non_equal)
    end

    @testset "Sum * product — diagonal split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        adj = IndexedOperator(a', j)

        sum_expr = Σ(ai, i)
        product = adj * aj
        result = sum_expr * product
        @test result isa QAdd
        @test any(p -> p == (i, j) || p == (j, i), result.non_equal)
    end

    @testset "QSym * Sum — diagonal split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # a†_j * Σ_i(a_i) — splits diagonal
        # Off-diag i ≠ j: a_i * a†_j (commute, then normal-order canonical).
        # Diag i = j: a†_j * a_j (already normal-ordered, no commutator remainder).
        sum_expr = Σ(ai, i)
        result = adj * sum_expr
        @test result isa QAdd
        @test length(result) == 2
        @test any(p -> p == (i, j) || p == (j, i), result.non_equal)
    end

    @testset "Sum * QSym — different space, no split" begin
        i = Index(h_prod, :i, 10, hn)
        σi = IndexedOperator(sigma, i)
        sum_expr = Σ(σi, i)

        # Σ_i(σ_i) * a — different spaces, no diagonal split
        result = sum_expr * a
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "QAdd * QAdd — clashing index error" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        s1 = Σ(ai, i)
        s2 = Σ(ai', i)
        @test_throws ArgumentError s1 * s2
    end

    @testset "QAdd * QAdd — different indices, same space splits" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        s1 = Σ(ai, i)
        s2 = Σ(adj, j)
        result = s1 * s2
        @test result isa QAdd
        @test length(result) > 1
        @test any(p -> p == (i, j) || p == (j, i), result.non_equal)
    end

    @testset "Legacy: NLevel diagonal splitting (DoubleSum equivalence)" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        i2 = Index(h2, :i, 10, ha2)
        j2 = Index(h2, :j, 10, ha2)
        σ2(α, β, k) = IndexedOperator(Transition(h2, :σ, α, β, 2), k)

        # Σ(Σ(σ_21_i * σ_12_j, i), j) splits into off-diagonal + diagonal
        inner = Σ(σ2(2, 1, i2) * σ2(1, 2, j2), i2)
        double = Σ(inner, j2)
        @test double isa QAdd

        # The diagonal term: σ(2,1,j) * σ(1,2,j) eagerly composes to σ(2,2,j)
        diag_ops = QSym[σ2(2, 2, j2)]
        @test haskey(double.arguments, diag_ops)
    end

    @testset "Legacy: Issue #221 — symbolic N double sum" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        @variables N_s::Real
        i2 = Index(h2, :i, N_s, ha2)
        j2 = Index(h2, :j, N_s, ha2)
        σ2(α, β, k) = IndexedOperator(Transition(h2, :σ, α, β, 2), k)

        # Σ(-σ_22_i, i, j) where j is independent → N * Σ(-σ_22_i, i)
        @test isequal(
            simplify(Σ(-σ2(2, 2, i2), i2, j2)),
            simplify(Σ(-σ2(2, 2, i2), i2) * N_s),
        )
    end

    @testset "Legacy: Issue #223 — commutators with sums" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        @variables N_c::Real
        i2 = Index(h2, :i, N_c, ha2)
        j2 = Index(h2, :j, N_c, ha2)
        k2 = Index(h2, :k, N_c, ha2)
        σ2(α, β, m) = IndexedOperator(Transition(h2, :σ, α, β, 2), m)

        H = Σ(2 * σ2(2, 2, j2), i2, j2)
        dict_N = Dict(SymbolicUtils.unwrap(N_c) => 10)
        sub_dict(x) = simplify(substitute(x, dict_N))

        @test isequal(sub_dict(simplify(commutator(H, σ2(2, 1, k2)))), simplify(20 * σ2(2, 1, k2)))
        @test isequal(sub_dict(simplify(commutator(H, σ2(1, 2, k2)))), simplify(-20 * σ2(1, 2, k2)))
    end

    @testset "Legacy: Issue #256 — cross-subspace sum * operator" begin
        @variables N_a::Real M_b::Real
        hA = NLevelSpace(:atomA, (:g, :r))
        hB = NLevelSpace(:atomB, (:g, :r))
        h2 = hA ⊗ hB

        σA(α, β, m) = IndexedOperator(Transition(h2, :σA, α, β, 1), m)
        σB(α, β, m) = IndexedOperator(Transition(h2, :σB, α, β, 2), m)

        i = Index(h2, :i, N_a, hA)
        j = Index(h2, :j, M_b, hB)
        k = Index(h2, :k, M_b, hB)

        op = σA(:g, :r, i)
        doublesum = Σ(σB(:r, :r, j) * σB(:r, :r, k), j, k)

        # Different subspaces: no diagonal split, no error
        result = op * doublesum
        @test result isa QAdd
        @test !iszero(result)

        # Commutator of different subspaces = 0
        @test iszero(simplify(commutator(op, doublesum)))
    end

    @testset "Legacy: Nested sum equivalences" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        i2 = Index(h2, :i, 4, ha2)
        j2 = Index(h2, :j, 4, ha2)
        σ2(α, β, k) = IndexedOperator(Transition(h2, :σ, α, β, 2), k)

        # Σ(Σ(σ_21_i * σ_12_j, i, [j]), j, [i]) vs Σ(Σ(σ_21_i * σ_12_j, i, [j]), j)
        # Adding redundant [i] to outer produces extra (j,i) constraint
        # but same terms and indices
        lhs = Σ(Σ(σ2(2, 1, i2) * σ2(1, 2, j2), i2, [j2]), j2, [i2])
        rhs = Σ(Σ(σ2(2, 1, i2) * σ2(1, 2, j2), i2, [j2]), j2)
        @test isequal(lhs.arguments, rhs.arguments)
        @test isequal(Set(lhs.indices), Set(rhs.indices))

        # Σ(Σ(σ_12_i * σ_21_j, i), j) splits into off-diagonal + diagonal
        inner = Σ(σ2(1, 2, i2) * σ2(2, 1, j2), i2)
        dsum = Σ(inner, j2)
        @test isequal(
            Σ(Σ(σ2(1, 2, i2) * σ2(2, 1, j2), i2, [j2]), j2) +
                Σ(σ2(1, 2, j2) * σ2(2, 1, j2), j2),
            dsum,
        )
    end

    @testset "Legacy: Issue #221 — symbolic N variants" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        @variables N1_s::Real c1_s::Real
        i2 = Index(h2, :i, N1_s, ha2)
        j2 = Index(h2, :j, N1_s, ha2)
        σ2(α, β, k) = IndexedOperator(Transition(h2, :σ, α, β, 2), k)

        # Σ(3*σ_22_i, i, j) where j independent → N * Σ(3*σ_22_i, i)
        @test isequal(
            simplify(Σ(3 * σ2(2, 2, i2), i2, j2)),
            simplify(Σ(3 * σ2(2, 2, i2), i2) * N1_s),
        )

        # Σ(c1*σ_22_i, i, j) where j independent → N * Σ(c1*σ_22_i, i)
        @test isequal(
            simplify(Σ(c1_s * σ2(2, 2, i2), i2, j2)),
            simplify(Σ(c1_s * σ2(2, 2, i2), i2) * N1_s),
        )
    end

    @testset "Legacy: Issue #223 — commutators with symbolic g" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        @variables N_g::Real g_s::Real
        i2 = Index(h2, :i, N_g, ha2)
        j2 = Index(h2, :j, N_g, ha2)
        k2 = Index(h2, :k, N_g, ha2)
        σ2(α, β, m) = IndexedOperator(Transition(h2, :σ, α, β, 2), m)

        dict_N = Dict(SymbolicUtils.unwrap(N_g) => 10)
        sub_dict(x) = simplify(substitute(x, dict_N))

        # H with symbolic prefactor g
        H_g = Σ(g_s * σ2(2, 2, j2), i2, j2)
        @test isequal(
            sub_dict(simplify(commutator(H_g, σ2(2, 1, k2)))),
            simplify(10 * g_s * σ2(2, 1, k2)),
        )
        @test isequal(
            sub_dict(simplify(commutator(H_g, σ2(1, 2, k2)))),
            simplify(-10 * g_s * σ2(1, 2, k2)),
        )

        # H with reversed index order (j first, then i)
        H_ji_g = Σ(g_s * σ2(2, 2, j2), j2, i2)
        @test isequal(
            sub_dict(simplify(commutator(H_ji_g, σ2(2, 1, k2)))),
            simplify(10 * g_s * σ2(2, 1, k2)),
        )
    end

    @testset "Legacy: Multi-mode multi-atom" begin
        ha2 = NLevelSpace(:atom, 2, 1)
        hf2 = FockSpace(:cavity)
        h2 = hf2 ⊗ ha2

        i_a = Index(h2, :i, 4, ha2)
        j_a = Index(h2, :j, 4, ha2)
        k_f = Index(h2, :k, 2, hf2)
        l_f = Index(h2, :l, 2, hf2)

        g_ik = DoubleIndexedVariable(:g, i_a, k_f)
        a2(k) = IndexedOperator(Destroy(h2, :a, 1), k)
        σ2(α, β, k) = IndexedOperator(Transition(h2, :σ, α, β, 2), k)

        Ssum1 = Σ(g_ik * a2(k_f) * σ2(2, 1, i_a), i_a)
        Ssum2 = Σ(conj(g_ik) * a2(k_f)' * σ2(1, 2, i_a), i_a)

        # adjoint distributes over sums
        @test isequal(adjoint(Ssum1), Ssum2)

        # Double sum nesting works across different spaces
        H = Σ(Ssum1 + Ssum2, k_f)
        @test H isa QAdd
        @test k_f in H.indices
        @test i_a in H.indices
    end

    # ========== Arithmetic with indexed operators ==========
    @testset "Indexed operator arithmetic" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # Multiplication
        m = adi * ai
        @test m isa QAdd
        @test length(operators(m)) == 2

        # Addition
        s = ai + adi
        @test s isa QAdd
        @test length(s) == 2

        # Scalar multiplication
        m2 = 2 * ai
        @test m2 isa QAdd
        @test prefactor(m2) == 2

        # Mixed indexed + non-indexed product
        m3 = a' * ai  # different sites
        @test m3 isa QAdd
        @test length(operators(m3)) == 2
    end

    @testset "Simplify indexed expressions" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # simplify(a_i * a†_i) should apply [a_i, a†_i] = 1
        result = simplify(ai * adi)
        @test result isa QAdd
        # Should be a†_i * a_i + 1
        @test length(result) == 2
    end

    @testset "Normal order indexed" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        result = normal_order(ai * adi)
        @test result isa QAdd
    end

    # ========== QAdd with indices through operations ==========
    @testset "QAdd sum algebra preserves indices" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        s1 = Σ(ai, i)
        s2 = s1 + a  # Add non-indexed term to sum
        @test s2 isa QAdd
        @test i in s2.indices

        # Multiply sum by scalar
        s3 = 2 * s1
        @test s3 isa QAdd
        @test i in s3.indices
    end

    @testset "QAdd * QSym preserves indices" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        s = Σ(ai, i)
        result = s * a'  # (Σ_i a_i) * a†
        @test result isa QAdd
        @test i in result.indices
    end

    @testset "QAdd * QAdd merges indices" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        s1 = Σ(ai, i)
        s2 = Σ(aj', j)
        result = s1 * s2  # (Σ_i a_i)(Σ_j a†_j)
        @test result isa QAdd
        @test i in result.indices
        @test j in result.indices
        # Same space → eager diagonal split produces off-diagonal + diagonal
        @test length(result) > 1
        @test any(p -> p == (i, j) || p == (j, i), result.non_equal)
    end

    @testset "Subtraction preserves indices" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        s = Σ(ai, i)
        neg = -s
        @test neg isa QAdd
        @test i in neg.indices
    end

    # ========== Commutator with indexed sums ==========
    @testset "Commutator — indexed QSym, QSym" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [a_i, a†_j] = 0 (different indices → different sites)
        result = commutator(ai, adj)
        @test iszero(result)
    end

    @testset "Commutator — sum with external operator (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [Σ_i a_i, a†_j] — sum collapses: only i=j term survives
        sum_expr = Σ(ai, i)
        result = commutator(sum_expr, adj)
        @test result isa QAdd
        # Should yield [a_j, a†_j] = 1
        simplified = simplify(result)
        @test length(simplified) == 1
        @test isempty(operators(only(sorted_arguments(simplified))))
        @test prefactor(only(sorted_arguments(simplified))) == 1
    end

    @testset "Commutator — external operator with sum (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [a_j, Σ_i a†_i] — sum collapses: only i=j term survives
        sum_expr = Σ(ai', i)
        result = commutator(adj', sum_expr)
        @test result isa QAdd
        # Should yield [a_j, a†_j] = 1
        simplified = simplify(result)
        @test length(simplified) == 1
        @test isempty(operators(only(sorted_arguments(simplified))))
        @test prefactor(only(sorted_arguments(simplified))) == 1
    end

    @testset "Commutator — sum with product (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [Σ_i a_i, a†_j * a_j] — sum over i, product has j-indexed ops
        sum_expr = Σ(ai, i)
        prod_expr = adj * IndexedOperator(a, j)
        result = commutator(sum_expr, prod_expr)
        @test result isa QAdd
    end

    @testset "Commutator — product with sum (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [a†_j * a_j, Σ_i a_i] — product has j-indexed ops, sum over i
        prod_expr = adj * IndexedOperator(a, j)
        sum_expr = Σ(ai, i)
        result = commutator(prod_expr, sum_expr)
        @test result isa QAdd
    end

    # ========== Single Hilbert space (non-ProductSpace) ==========
    @testset "Index on single HilbertSpace" begin
        N = 10
        h = NLevelSpace(:atom, 2, 1)
        i = Index(h, :i, N, h)

        σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
        @test σ(1, 2, i) == (σ(1, 2, i)')'
    end

    # ========== Mixed operator types (NLevel × NLevel) ==========
    @testset "Mixed NLevel operator types across subspaces" begin
        hc = NLevelSpace(:cavity, 3, 1)
        ha = NLevelSpace(:atom, 2, 1)
        h = hc ⊗ ha

        @variables N::Real
        i = Index(h, :i, N, ha)

        S(x, y) = Transition(h, :S, x, y, 1)
        σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)

        @test S(2, 1) * σ(1, 2, i) isa QAdd
        @test σ(1, 2, i) * S(2, 1) isa QAdd
        @test isequal(S(2, 1) * σ(1, 2, i), σ(1, 2, i) * S(2, 1))
    end

    # ========== Indexed variable arithmetic ==========
    @testset "IndexedVariable in prefactor" begin
        i = Index(h_prod, :i, 10, hn)
        gi = IndexedVariable(:g, i)

        # g(i)*a puts indexed variable in prefactor, operator in operators()
        m = gi * a
        @test m isa QAdd
        @test isequal(operators(m), [a])

        # a*g(i) same
        m2 = a * gi
        @test m2 isa QAdd
        @test isequal(operators(m2), [a])

        # g(i)*a†
        m3 = gi * a'
        @test m3 isa QAdd
        @test isequal(operators(m3), [a'])

        # a†*g(i)
        m4 = a' * gi
        @test m4 isa QAdd
        @test isequal(operators(m4), [a'])

        # σ*g(i)
        σi = IndexedOperator(sigma, i)
        m5 = σi * gi
        @test m5 isa QAdd
        @test isequal(operators(m5), [σi])

        # g(i)*σ
        m6 = gi * σi
        @test m6 isa QAdd
        @test isequal(operators(m6), [σi])

        # g(i)*QAdd preserves operators
        qmul = a' * a
        m7 = gi * qmul
        @test length(operators(m7)) == 2
        m8 = qmul * gi
        @test length(operators(m8)) == 2
    end

    @testset "IndexedVariable addition commutativity" begin
        i = Index(h_prod, :i, 10, hn)
        j = Index(h_prod, :j, 10, hn)
        gi = IndexedVariable(:g, i)
        σi = IndexedOperator(sigma, i)
        σj = IndexedOperator(sigma, j)

        @test isequal(gi + a, a + gi)
        @test isequal(a + σi, σi + a)
        @test isequal(σj + σi, σi + σj)

        qadd = a + a'
        @test isequal(gi + qadd, qadd + gi)
        @test length(qadd + gi) == 3
        @test length(qadd + σi) == 3
        @test isequal(σi + qadd, qadd + σi)

        qmul = a' * a
        @test isequal(gi + qmul, qmul + gi)
        @test isequal(gi + σj, σj + gi)
    end

    @testset "Negation of indexed operators" begin
        i = Index(h_prod, :i, 10, hn)
        gi = IndexedVariable(:g, i)
        σi = IndexedOperator(sigma, i)

        @test isequal(-σi, -1 * σi)
        @test isequal(-gi, -1 * gi)
    end

    # ========== Indexed commutators across spaces ==========
    @testset "Indexed commutator — different spaces vanish" begin
        i = Index(h_prod, :i, 10, hn)
        σi = IndexedOperator(sigma, i)
        qadd = a + a'
        qmul = a' * a

        # σ on NLevel space, a on Fock space → different spaces → commutator = 0
        @test iszero(simplify(commutator(σi, qadd)))
        @test iszero(simplify(commutator(σi, qmul)))
    end

    # ========== Same-index Fock normal ordering ==========
    @testset "Same-index Fock: a_m · a_m† = a_m† · a_m + 1" begin
        m = Index(hf, :m, 10, hf)
        ai = IndexedOperator(a, m)
        result = simplify(normal_order(ai * ai'))
        expected = simplify(ai' * ai + 1)
        @test isequal(result, expected)
    end

    # ========== Same-site transition product = 0 ==========
    @testset "Same-site σ₁₂·σ₁₂ = 0 via normal_order" begin
        i = Index(h_prod, :i, 10, hn)
        σi = IndexedOperator(sigma, i)
        result = simplify(normal_order(σi * σi))
        @test iszero(result)
    end

    # ========== change_index producing zero ==========
    @testset "change_index — DoubleIndexedVariable identical=false → 0" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        Ωij = DoubleIndexedVariable(:Ω, i, j; identical = false)
        # Construction with equal indices gives 0
        @test isequal(DoubleIndexedVariable(:Ω, i, i; identical = false), Symbolics.Num(0))
        # change_index collapsing to equal indices gives 0
        @test isequal(change_index(Ωij, i, j), Symbolics.Num(0))
        @test isequal(change_index(Ωij, j, i), Symbolics.Num(0))
        # change_index on identical=true (default) keeps the value
        Ωij_diag = DoubleIndexedVariable(:Ω, i, j)
        @test !isequal(change_index(Ωij_diag, i, j), Symbolics.Num(0))
        # change_index through QAdd with identical=false prefactor
        ai = IndexedOperator(a, i)
        m = Ωij * ai
        mj = change_index(m, i, j)
        @test iszero(mj)
    end

    # ========== Single sums ==========
    @testset "Single sum operations" begin
        N = 10
        i = Index(h_prod, :i, N, hn)
        j = Index(h_prod, :j, N, hn)
        gi = IndexedVariable(:g, i)
        gj = IndexedVariable(:g, j)
        σi(x, y) = IndexedOperator(Transition(h_prod, :σ, x, y, 2), i)
        Γij = DoubleIndexedVariable(:Γ, i, j)

        s1 = Σ(σi(1, 2) * a', i)
        s2 = Σ(IndexedOperator(Transition(h_prod, :σ, 2, 1, 2), i) * a, i)
        s3 = Σ(a' * σi(1, 2) + a * IndexedOperator(Transition(h_prod, :σ, 2, 1, 2), i), i)

        # adjoint distributes over sums
        @test isequal(adjoint(s1), s2)
        # sum of sum = sum of terms
        @test isequal(s3, s1 + s2)

        # Σ(0, i) stays as a QAdd (lazy)
        @test Σ(0, i) isa QAdd

        # Sum of variable over unrelated index → N * variable
        r_gj = Σ(gj, i)
        @test isempty(get_indices(r_gj))
        @test isequal(simplify(r_gj * a), simplify(N * gj * a))
        r_Γij = Σ(Γij, Index(h_prod, :k, N, hn))
        @test isempty(get_indices(r_Γij))
        @test isequal(simplify(r_Γij * a), simplify(N * Γij * a))

        # ∑ and Σ equivalence
        @test isequal(∑(σi(1, 2), i), Σ(σi(1, 2), i))
        @test isequal(
            ∑(σi(1, 2) * IndexedOperator(Transition(h_prod, :σ, 2, 1, 2), j), i, j),
            Σ(σi(1, 2) * IndexedOperator(Transition(h_prod, :σ, 2, 1, 2), j), i, j),
        )

        # change_index on ∑
        @test isequal(
            change_index(∑(2gj, j), j, i),
            ∑(2gi, i),
        )
    end

    @testset "Sum addition" begin
        N = 10
        i = Index(h_prod, :i, N, hn)
        σi = IndexedOperator(sigma, i)
        s1 = Σ(σi * a', i)

        # Sum + operator
        @test (s1 + a') isa QAdd
        @test (s1 + σi) isa QAdd

        # Sum + QAdd
        qadd = a + a'
        @test length(qadd + s1) == 3
        @test isequal(s1 + qadd, qadd + s1)

        # Sum + QAdd (product)
        qmul = a' * a
        @test s1 + qmul isa QAdd
    end

    # ========== Sum of a constant ==========
    @testset "Sum of constant over index" begin
        @variables α_sum::Real
        N = 10
        i = Index(h_prod, :i, N, hn)
        # Constant (no dependence on i) → simplified to N * α
        s = Σ(Symbolics.Num(α_sum), i)
        @test s isa QAdd
        @test isempty(get_indices(s))
        @test isequal(simplify(s * a), simplify(N * α_sum * a))
    end

    # ========== Average addition with indexed sums (PR 28) ==========
    @testset "Average addition with indexed sums (PR 28)" begin
        ha = NLevelSpace(:atom1, 2, 1)
        hb = NLevelSpace(:atom2, 2, 1)
        h = ha ⊗ hb
        @variables N
        i1 = Index(h, :i1, N, ha)
        i2 = Index(h, :i2, N, hb)
        s1(α, β, i) = IndexedOperator(Transition(h, :S1, α, β, 1), i)
        s2(α, β, i) = IndexedOperator(Transition(h, :S2, α, β, 2), i)

        term = average(Σ(s1(2, 1, i1) * s2(1, 2, i2), i1, i2))
        @test (term + term) isa SymbolicUtils.BasicSymbolic
    end

    # ========== conj on indexed variables ==========
    @testset "conj on indexed variables" begin
        i = Index(hf, :i, 10, hf)
        gi = IndexedVariable(:g, i)  # Real-typed
        # conj of real indexed variable is identity
        @test isequal(conj(gi), gi)

        # conj distributes through products with indexed variables
        ai = IndexedOperator(a, i)
        expr = gi * ai
        adj_expr = adjoint(expr)
        # adjoint of g_i * a_i should be g_i * a_i† (g_i is real)
        @test adj_expr isa QAdd
        @test operators(adj_expr)[1] == ai'
    end

    # ========== create_index_arrays ==========
    @testset "create_index_arrays" begin
        i = Index(h_prod, :i, 10, hn)
        j = Index(h_prod, :j, 5, hn)

        # Single index → returns a collected vector
        arr1 = create_index_arrays([i], [1:10])
        @test arr1 isa Vector{Int}
        @test arr1 == collect(1:10)

        # Multi-index → flat Cartesian product
        ranges = [1:10, 1:5]
        arr2 = create_index_arrays([i, j], ranges)
        @test isequal(arr2, vec(collect(Iterators.product(ranges...))))
    end

    # ========== Type stability ==========
    @testset "Type stability" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        @inferred Σ(ai, i)
        @inferred change_index(ai, i, Index(hf, :j, 10, hf))
        @inferred get_indices(ai)
        @inferred commutator(ai, adi)
    end
end
