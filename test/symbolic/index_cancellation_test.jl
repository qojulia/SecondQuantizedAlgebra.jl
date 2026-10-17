using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: is_indexed_sum
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: constraint_pairs, get_sum_indices, get_sum_non_equal,
    has_sum_metadata
import SecondQuantizedAlgebra as SQA

@testset "Index scope and cancellation" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
    a = Destroy(h, :a, 1)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
    @variables N g Δ
    i = Index(h, :i, N, 2)
    j = Index(h, :j, N, 2)
    k = Index(h, :k, N, 2)

    @testset "cancellation removes dead summation scope" begin
        s = Σ(g * σ(1, 2, i), i)
        @test iszero(s - s)
        @test isempty(get_indices(s - s))
        @test isempty(get_indices(s + (-s)))

        s_i = Σ(g * σ(1, 2, i), i)
        s_j = Σ(g * σ(1, 2, j), j)
        @test Set(get_indices(s_i - s_j)) == Set([i, j])
        @test get_indices((s_i + s_j) - s_i) == [j]
    end

    @testset "sum residuals carry their range" begin
        hf = FockSpace(:modes)
        b = Destroy(hf, :b)
        m = Index(hf, :m, N, hf)
        n = Index(hf, :n, N, hf)
        bm = IndexedOperator(b, m)
        bn = IndexedOperator(b, n)
        @test iszero(simplify(Σ(bm + 1, m) - (Σ(bm, m) + N)))
        @test iszero(simplify(Σ(bm * bm', m) - (Σ(bm' * bm, m) + N)))
        @test iszero(simplify(commutator(Σ(bm, m), Σ(bn', n)) - N))
    end

    @testset "nested constrained sums count every binder" begin
        hf = FockSpace(:modes)
        b = Destroy(hf, :b)
        n = 5
        m = Index(hf, :m, n, hf)
        p = Index(hf, :p, n, hf)
        q = Index(hf, :q, n, hf)
        local k = Index(hf, :k, n, hf)
        one_op = 1 + 0 * b' * b
        @test iszero(Σ(Σ(Σ(one_op, m, [p, q]), p), q) - 80)
        c2 = Σ(Σ(one_op, m, [p, q]), p)
        @test iszero(c2 - 16)

        bk = IndexedOperator(b, k)
        bp = IndexedOperator(b, p)
        kept = Σ(Σ(bk' * bk, k, [p, q]), p)
        @test get_indices(kept) == [k]
        @test iszero(kept - Σ(4 * bk' * bk, k, [q]))

        split = Σ(Σ(bk' * bp, k, [q]), p)
        @test iszero(split - (Σ(Σ(bk' * bp, k, [p, q]), p) + Σ(bp' * bp, p, [q])))
        @test all(t -> !isempty(get_indices(t)), SecondQuantizedAlgebra.sorted_arguments(split))
        @test occursin("p≠q", string(split))
    end

    @testset "diagonal and off-diagonal terms remain distinct" begin
        # The public constraint view exposes the split without depending on how terms are stored.
        off = Σ(σ(1, 2, i) * σ(2, 1, j), i)
        @test any(pair -> pair == (i, j) || pair == (j, i), constraint_pairs(off))
        @test iszero(
            simplify(
                assume_distinct_index(σ(2, 1, i) * σ(1, 2, j), [(i, j)]) -
                    assume_distinct_index(σ(2, 1, i) * σ(1, 2, j), [(j, i)]),
            )
        )
        # Keep a genuine indexed Pauli product in the public suite as well.
        hp = PauliSpace(:pauli_scope)
        ip = Index(hp, :ip, N, hp)
        jp = Index(hp, :jp, N, hp)
        σxi = IndexedOperator(Pauli(hp, :σ, 1), ip)
        σyj = IndexedOperator(Pauli(hp, :σ, 2), jp)
        σzi = IndexedOperator(Pauli(hp, :σ, 3), ip)
        actual = Σ(σxi, ip) * Σ(σyj, jp)
        diagonal = im * Σ(σzi, ip)
        off_diagonal = Σ(Σ(σxi * σyj, ip, [jp]), jp)
        @test iszero(simplify(actual - diagonal - off_diagonal))
        @test any(pair -> pair == (ip, jp) || pair == (jp, ip), constraint_pairs(actual))
    end

    @testset "constrained sums count excluded sites" begin
        hf = FockSpace(:counted)
        b = Destroy(hf, :b)
        M = 5
        ii = Index(hf, :i, M, hf)
        p = Index(hf, :p, M, hf)
        q = Index(hf, :q, M, hf)
        bi = IndexedOperator(b, ii)
        bp = IndexedOperator(b, p)
        bq = IndexedOperator(b, q)

        @test iszero(simplify(Σ(bi * bi', ii) - (M + Σ(bi' * bi, ii))))
        @test iszero(simplify(Σ(bi * bi', ii, [p]) - ((M - 1) + Σ(bi' * bi, ii, [p]))))
        @test iszero(
            simplify(
                Σ(bi * bi' * bp' * bp, ii, [p]) -
                    ((M - 1) * bp' * bp + Σ(bi' * bi * bp' * bp, ii, [p])),
            ),
        )
        @test iszero(simplify(Σ(Σ(bi * bi', ii, [p]), p) - (M * (M - 1) + (M - 1) * Σ(bi' * bi, ii))))

        # p and q may coincide, so the count M - |{p, q}| is unknown until declared.
        undecided = Σ(bi * bi' * bp' * bq, ii, [p, q])
        @test ii in undecided.indices
        distinct = assume_distinct_index(undecided, [(p, q)])
        expected = assume_distinct_index(
            (M - 2) * bp' * bq + Σ(bi' * bi * bp' * bq, ii, [p, q]), [(p, q)],
        )
        @test iszero(simplify(distinct - expected))
    end

    @testset "average only retains live sum metadata" begin
        averaged = average(Σ(σ(1, 2, i), i))
        @test is_indexed_sum(averaged)
        @test has_sum_metadata(averaged)
        @test get_sum_indices(averaged) == [i]
        @test isempty(get_sum_non_equal(averaged))
        @test iszero(undo_average(averaged) - Σ(σ(1, 2, i), i))

        hf = FockSpace(:average_modes)
        b = Destroy(hf, :b)
        m = Index(hf, :m, N, hf)
        n = Index(hf, :n, N, hf)
        comm = commutator(Σ(IndexedOperator(b, m), m), Σ(IndexedOperator(b', n), n))
        @test !has_sum_metadata(average(comm))
    end

    @testset "indexed equations use public semantics" begin
        H = -Δ * a' * a + Σ(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
        @test iszero(simplify(commutator(H, a) - (Δ * a - g * Σ(σ(1, 2, i), i))))
        @test iszero(simplify(commutator(-Δ * a' * a, σ(1, 2, i))))
        @test iszero(simplify(H - adjoint(H)))

        leibniz = commutator(H, σ(1, 2, j)) * σ(2, 1, k) +
            σ(1, 2, j) * commutator(H, σ(2, 1, k))
        direct = commutator(H, σ(1, 2, j) * σ(2, 1, k))
        @test iszero(
            simplify(
                assume_distinct_index(direct - leibniz, [(j, k)]),
            )
        )
    end

    @testset "binder elimination reaches a fixpoint" begin
        # Σ_{i,k,l} [i≠k, i≠l] 1 = N(N-1)²: `i` is countable only after `k` and `l` are gone.
        h5 = NLevelSpace(:atom, 2)
        i5, k5, l5 = Index(h5, :i, 5, h5), Index(h5, :k, 5, h5), Index(h5, :l, 5, h5)
        ne = SQA.NonEqualPair[(i5, k5), (i5, l5)]
        temp = SQA.QTermDict(SQA.QTerm(SQA.Op[], ne) => SQA.CNUM_ONE)
        out = SQA.settle_scope!(SQA.QTermDict(), temp, [i5, k5, l5])
        @test only(values(out)) == SQA.to_cnum(80)
        @test isempty(only(keys(out)).ne)

        stemp = SQA.SuperTermDict(SQA.SuperTerm(SQA.Op[], SQA.Op[], ne) => SQA.CNUM_ONE)
        sout = SQA.settle_scope!(SQA.SuperTermDict(), stemp, [i5, k5, l5])
        @test only(values(sout)) == SQA.to_cnum(80)

        J = SQA.single_qadd(SQA.CNUM_ONE, SQA.Op[Transition(h5, :σ, 1, 2)])
        kept, cne, rate = SQA.count_unmentioned_binders(J, J, SQA.CNUM_ONE, [i5, k5, l5], ne)
        @test isempty(kept) && isempty(cne) && rate == SQA.to_cnum(80)
    end

    @testset "diagonal split counts coinciding externals once" begin
        # Σᵢ over i ∈ {k, l}: setting l = k afterwards must give the sum built at k = l.
        s22(q) = σ(2, 2, q)
        l = Index(h, :l, N, 2)
        X = Σ(s22(i) * s22(k) * s22(l), i)
        @test iszero(simplify(change_index(X, l, k) - Σ(s22(i) * s22(k), i)))
        P = Σ(s22(i), i) * (s22(k) * s22(l))
        @test iszero(simplify(change_index(P, l, k) - Σ(s22(i), i) * s22(k)))
        @test iszero(simplify(P - X))
    end
end
