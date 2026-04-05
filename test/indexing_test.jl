using SecondQuantizedAlgebra
using Test
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: QMul, QAdd, QSym, CNum, _to_cnum, NO_INDEX

@testset "Indexing" begin

    # ========== Setup ==========
    hf = FockSpace(:f)
    hn = NLevelSpace(:n, 2, 1)
    hp = PauliSpace(:p)
    hs = SpinSpace(:s, 1 // 2)
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
        @test all(iszero, result.arguments)
    end

    @testset "Indexed operators on same index interact" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # Same index → same site → [a_i, a†_i] = 1
        result = commutator(ai, adi)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isempty(result.arguments[1].args_nc)
        @test result.arguments[1].arg_c == 1
    end

    @testset "Indexed vs non-indexed are different sites" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)

        # Non-indexed a and indexed a_i have different indices (NO_INDEX vs i)
        result = commutator(a, ai')
        @test all(iszero, result.arguments)
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
        @test expr isa QMul
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

    @testset "Sigma with QMul" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        s = Σ(adi * ai, i)
        @test s isa QAdd
        @test length(s.indices) == 1
        @test length(s.arguments) == 1
        @test length(s.arguments[1].args_nc) == 2
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

        s = Σ(ai, i, j)
        @test s isa QAdd
        @test length(s.indices) == 2
        @test i in s.indices
        @test j in s.indices
    end

    @testset "Sigma with scalar" begin
        i = Index(hf, :i, 10, hf)
        s = Σ(1, i)
        @test s isa QAdd
        @test length(s.indices) == 1
    end

    @testset "Sigma with QSym" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        s = Σ(ai, i)
        @test s isa QAdd
        @test length(s.arguments) == 1
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

    @testset "change_index — QMul" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        m = adi * ai  # a†_i * a_i
        mj = change_index(m, i, j)
        @test mj isa QMul
        for op in mj.args_nc
            @test op.index == j
        end
    end

    @testset "change_index — QMul prefactor with IndexedVariable" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        gi = IndexedVariable(:g, i)
        ai = IndexedOperator(a, i)

        m = gi * ai  # g(i) * a_i
        mj = change_index(m, i, j)
        @test mj.args_nc[1].index == j
        # Prefactor should also be substituted: g(i) → g(j)
        gj = IndexedVariable(:g, j)
        @test isequal(real(mj.arg_c), gj)
    end

    @testset "change_index — QAdd" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        s = ai + adi  # a_i + a†_i
        sj = change_index(s, i, j)
        @test sj isa QAdd
        for term in sj.arguments
            for op in term.args_nc
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

        # QMul with indexed operators
        aj = IndexedOperator(a, j)
        m = ai' * aj
        inds = get_indices(m)
        @test length(inds) == 2
        @test i in inds
        @test j in inds

        # QMul with repeated index → unique
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

    # ========== expand_sums ==========
    @testset "expand_sums — passthrough for non-sums" begin
        @test expand_sums(a) === a
        m = a' * a
        @test expand_sums(m) === m
    end

    @testset "expand_sums — no indices → passthrough" begin
        s = a + a'
        @test expand_sums(s) === s
    end

    @testset "expand_sums — diagonal split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ_i(a†_j * a_i) — sum over i, j is external
        # Should split into: a†_j * a_j (diagonal, i→j) + Σ_{i≠j}(a†_j * a_i)
        expr = Σ(adj * ai, i)
        expanded = expand_sums(expr)
        @test expanded isa QAdd
        # Should have more terms than original
        @test length(expanded.arguments) > length(expr.arguments)
        # Should have i≠j constraint
        @test any(p -> p == (i, j) || p == (j, i), expanded.non_equal)
    end

    @testset "expand_sums — already non_equal → no split" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # Σ_{i≠j}(a†_j * a_i) — already constrained
        expr = Σ(adj * ai, i, [j])
        expanded = expand_sums(expr)
        # No new terms should be added (already non_equal)
        @test length(expanded.arguments) == length(expr.arguments)
    end

    # ========== Arithmetic with indexed operators ==========
    @testset "Indexed operator arithmetic" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # Multiplication
        m = adi * ai
        @test m isa QMul
        @test length(m.args_nc) == 2

        # Addition
        s = ai + adi
        @test s isa QAdd
        @test length(s.arguments) == 2

        # Scalar multiplication
        m2 = 2 * ai
        @test m2 isa QMul
        @test m2.arg_c == 2

        # Mixed indexed + non-indexed product
        m3 = a' * ai  # different sites
        @test m3 isa QMul
        @test length(m3.args_nc) == 2
    end

    @testset "Simplify indexed expressions" begin
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        adi = ai'

        # simplify(a_i * a†_i) should apply [a_i, a†_i] = 1
        result = simplify(ai * adi)
        @test result isa QAdd
        # Should be a†_i * a_i + 1
        @test length(result.arguments) == 2
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
        @test all(iszero, result.arguments)
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
        @test length(simplified.arguments) == 1
        @test isempty(simplified.arguments[1].args_nc)
        @test simplified.arguments[1].arg_c == 1
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
        @test length(simplified.arguments) == 1
        @test isempty(simplified.arguments[1].args_nc)
        @test simplified.arguments[1].arg_c == 1
    end

    @testset "Commutator — sum with QMul (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [Σ_i a_i, a†_j * a_j] — sum over i, QMul has j-indexed ops
        sum_expr = Σ(ai, i)
        qmul_expr = adj * IndexedOperator(a, j)
        result = commutator(sum_expr, qmul_expr)
        @test result isa QAdd
    end

    @testset "Commutator — QMul with sum (diagonal collapse)" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(a', j)

        # [a†_j * a_j, Σ_i a_i] — QMul has j-indexed ops, sum over i
        qmul_expr = adj * IndexedOperator(a, j)
        sum_expr = Σ(ai, i)
        result = commutator(qmul_expr, sum_expr)
        @test result isa QAdd
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
