using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: constraint_pairs, get_sum_indices, get_sum_non_equal,
    has_sum_metadata

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
end
