using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym, CNum, _to_cnum, NO_INDEX

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
        # Scalar 1 doesn't depend on i → simplified to 10 * 1
        @test isempty(get_indices(s))
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

        @test S(2, 1) * σ(1, 2, i) isa QMul
        @test σ(1, 2, i) * S(2, 1) isa QMul
        @test isequal(S(2, 1) * σ(1, 2, i), σ(1, 2, i) * S(2, 1))
    end

    # ========== Indexed variable arithmetic ==========
    @testset "IndexedVariable in prefactor" begin
        i = Index(h_prod, :i, 10, hn)
        gi = IndexedVariable(:g, i)

        # g(i)*a puts indexed variable in prefactor, operator in args_nc
        m = gi * a
        @test m isa QMul
        @test isequal(m.args_nc, [a])

        # a*g(i) same
        m2 = a * gi
        @test m2 isa QMul
        @test isequal(m2.args_nc, [a])

        # g(i)*a†
        m3 = gi * a'
        @test m3 isa QMul
        @test isequal(m3.args_nc, [a'])

        # a†*g(i)
        m4 = a' * gi
        @test m4 isa QMul
        @test isequal(m4.args_nc, [a'])

        # σ*g(i)
        σi = IndexedOperator(sigma, i)
        m5 = σi * gi
        @test m5 isa QMul
        @test isequal(m5.args_nc, [σi])

        # g(i)*σ
        m6 = gi * σi
        @test m6 isa QMul
        @test isequal(m6.args_nc, [σi])

        # g(i)*QMul preserves operators
        qmul = a' * a
        m7 = gi * qmul
        @test length(m7.args_nc) == 2
        m8 = qmul * gi
        @test length(m8.args_nc) == 2
    end

    @testset "IndexedVariable addition commutativity" begin
        i = Index(h_prod, :i, 10, hn)
        j = Index(h_prod, :j, 10, hn)
        gi = IndexedVariable(:g, i)
        σi = IndexedOperator(sigma, i)
        σj = IndexedOperator(sigma, j)

        @test isequal(gi + a, a + gi)
        # TODO: QAdd isequal is not order-independent for indexed operators
        @test_broken isequal(a + σi, σi + a)
        @test_broken isequal(σj + σi, σi + σj)

        qadd = a + a'
        @test isequal(gi + qadd, qadd + gi)
        @test length((qadd + gi).arguments) == 3
        @test length((qadd + σi).arguments) == 3
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
        @test all(iszero, simplify(commutator(σi, qadd)).arguments)
        @test all(iszero, simplify(commutator(σi, qmul)).arguments)
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
        @test all(iszero, result.arguments)
    end

    # ========== change_index producing zero ==========
    @testset "change_index — DoubleIndexedVariable identical=false → 0" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        Ωij = DoubleIndexedVariable(:Ω, i, j; identical = false)
        # When both indices become the same, identical=false should give 0
        @test isequal(DoubleIndexedVariable(:Ω, i, i; identical = false), Symbolics.Num(0))
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
        @test length((qadd + s1).arguments) == 3
        # TODO: QAdd isequal is not order-independent for sums
        @test_broken isequal(s1 + qadd, qadd + s1)

        # Sum + QMul
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
        @test adj_expr isa QMul
        @test adj_expr.args_nc[1] == ai'
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
