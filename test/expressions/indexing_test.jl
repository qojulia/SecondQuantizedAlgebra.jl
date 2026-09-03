using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: constraint_pairs

@testset "Indexed public API" begin
    hf = FockSpace(:fock)
    hn = NLevelSpace(:atom, 2, 1)
    hp = PauliSpace(:pauli)
    hs = SpinSpace(:spin)
    hq = PhaseSpace(:phase)

    @testset "index values and scope" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        @test index_name(i) == :i
        @test index_range(i) == 10
        @test isequal(index_sym(i), index_sym(Index(hf, :i, 10, hf)))
        @test acts_on(i) == [1]
        @test has_index(i)
        @test i == Index(hf, :i, 10, hf)
        @test i != j
        @test i != Index(hf, :i, 5, hf)
        @test hash(i) == hash(Index(hf, :i, 10, hf))

        integer_space = Index(hf, :m, 10, 1)
        @test acts_on(integer_space) == [1]
        @test has_index(integer_space)

        unindexed = Destroy(hf, :unindexed)
        @test !has_index(operator_index(unindexed))

        concrete = i(3)
        @test has_index(concrete)
        @test index_slot(concrete) == 3
        @test index_name(concrete) == :i_3
        @test index_range(concrete) == index_range(i)
        @test i(3) == i(3)
        @test i(3) != i(4)

        product = hf ⊗ hn
        @test acts_on(Index(product, :m, 4, hf)) == [1]
        @test acts_on(Index(product, :n, 4, hn)) == [2]
        @test_throws ArgumentError Index(product, :x, 4, FockSpace(:other))
        @test_throws ArgumentError Index(hf, "i", 4, hf)
    end

    @testset "index-aware operators cover every family" begin
        i = Index(hf, :i, 10, hf)
        @test is_destroy(IndexedOperator(Destroy(hf, :a), i))
        @test is_create(IndexedOperator(Create(hf, :a), i))

        j = Index(hn, :j, 10, hn)
        @test is_transition(IndexedOperator(Transition(hn, :σ, 1, 2), j))

        pidx = Index(hp, :p, 10, hp)
        sidx = Index(hs, :s, 10, hs)
        qidx = Index(hq, :q, 10, hq)
        @test is_pauli(IndexedOperator(Pauli(hp, :σ, 1), pidx))
        @test is_spin(IndexedOperator(Spin(hs, :S, 1), sidx))
        @test is_position(IndexedOperator(Position(hq, :x), qidx))
        @test is_momentum(IndexedOperator(Momentum(hq, :p), qidx))

        indexed = IndexedOperator(Destroy(hf, :a), i)
        @test operator_index(indexed) == i
        @test operator_index(indexed') == i
        @test is_create(indexed')
        @test acts_on(indexed) == [1]
    end

    @testset "change_index covers every operator family" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)

        @test operator_index(change_index(IndexedOperator(Destroy(hf, :a), i), i, j)) == j
        @test operator_index(change_index(IndexedOperator(Create(hf, :a), i), i, j)) == j

        ni = Index(hn, :ni, 10, hn)
        nj = Index(hn, :nj, 10, hn)
        @test operator_index(
            change_index(IndexedOperator(Transition(hn, :σ, 1, 2), ni), ni, nj),
        ) == nj

        pi = Index(hp, :pi, 10, hp)
        pj = Index(hp, :pj, 10, hp)
        @test operator_index(change_index(IndexedOperator(Pauli(hp, :σ, 1), pi), pi, pj)) == pj

        si = Index(hs, :si, 10, hs)
        sj = Index(hs, :sj, 10, hs)
        @test operator_index(change_index(IndexedOperator(Spin(hs, :S, 1), si), si, sj)) == sj

        qi = Index(hq, :qi, 10, hq)
        qj = Index(hq, :qj, 10, hq)
        @test operator_index(change_index(IndexedOperator(Position(hq, :x), qi), qi, qj)) == qj
        @test operator_index(change_index(IndexedOperator(Momentum(hq, :p), qi), qi, qj)) == qj
    end

    @testset "indexed parameters" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        gi = IndexedVariable(:g, i)
        gj = IndexedVariable(:g, j)
        @test !isequal(gi, gj)
        @test !iszero(gi)

        Jij = DoubleIndexedVariable(:J, i, j)
        @test !iszero(Jij)
        @test iszero(DoubleIndexedVariable(:J, i, i; identical = false))
        @test !iszero(DoubleIndexedVariable(:J, i, i))
    end

    @testset "get_indices reports public expression scope" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(Destroy(hf, :a), i)
        aj = IndexedOperator(Destroy(hf, :a), j)

        @test isempty(get_indices(42))
        @test isempty(get_indices(Destroy(hf, :a)))
        @test get_indices(ai) == [i]
        @test Set(get_indices(ai' * aj)) == Set([i, j])
        @test get_indices(ai' * ai) == [i]
        @test Set(get_indices(ai + aj)) == Set([i, j])
        @test get_indices(Σ(ai, i)) == [i]
    end

    @testset "same-index and distinct-index algebra" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        a = Destroy(hf, :a)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        @test iszero(commutator(ai, aj'))
        @test iszero(commutator(ai, ai') - 1)
        @test iszero(commutator(a, ai'))
        @test iszero(normal_order(ai * ai') - (ai' * ai + 1))

        distinct = assume_distinct_index(ai * aj', [(i, j)])
        @test any(pair -> pair == (i, j) || pair == (j, i), constraint_pairs(distinct))
        @test operators(distinct) == operators(ai * aj')
        @test !isequal(distinct, ai * aj')

        hproduct = hf ⊗ hn
        cavity = Destroy(hproduct, :a, 1)
        σ = IndexedOperator(Transition(hproduct, :σ, 1, 2, 2), Index(hproduct, :k, 10, 2))
        @test iszero(commutator(cavity, σ))
    end

    @testset "index constraints control same-site reductions" begin
        h = NLevelSpace(:atom, 2, 1)
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        σ12(k) = IndexedOperator(Transition(h, :σ, 1, 2), k)
        σ11(k) = IndexedOperator(Transition(h, :σ, 1, 1), k)

        same_site = σ12(i) * σ12(i)'
        different_sites = assume_distinct_index(
            σ12(i) * σ12(j)', [(i, j)],
        )

        @test iszero(same_site - σ11(i))
        @test !iszero(different_sites - σ11(i))
        @test any(pair -> pair == (i, j) || pair == (j, i), constraint_pairs(different_sites))
    end

    @testset "summation scope and diagonal splitting" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        single = Σ(ai, i)
        @test get_indices(single) == [i]
        @test isempty(constraint_pairs(single))
        @test isequal(∑(ai, i), single)
        @test iszero(Σ(1, i) - 10)
        @test isequal(Σ(ai, i, j), 10 * single)

        off_diagonal = Σ(aj' * ai, i)
        @test i in get_indices(off_diagonal)
        @test j in get_indices(off_diagonal)
        @test any(pair -> pair == (i, j) || pair == (j, i), constraint_pairs(off_diagonal))

        already_distinct = Σ(aj' * ai, i, [j])
        @test any(pair -> pair == (i, j) || pair == (j, i), constraint_pairs(already_distinct))

        @test iszero(commutator(single, aj') - 1)
        @test iszero(commutator(Σ(ai, i), Σ(aj', j)) - 10)

        hmix = hf ⊗ hn
        atom_index = Index(hmix, :k, 10, hn)
        other_index = Index(hmix, :m, 10, hf)
        atom = IndexedOperator(Transition(hmix, :σ, 1, 2, 2), atom_index)
        @test isequal(Σ(atom, atom_index, other_index), Σ(atom, atom_index) * 10)
    end

    @testset "public sum products split only physical diagonals" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        adj = adjoint(aj)

        product_sum = Σ(adj * ai, i)
        @test length(product_sum) == 2
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(product_sum),
        )
        @test haskey(product_sum, operators(adj * aj))

        constrained = Σ(adj * ai, i, [j])
        @test length(constrained) == 1
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(constrained),
        )

        right_product = Σ(ai, i) * adj
        left_product = adj * Σ(ai, i)
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(right_product),
        )
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(left_product),
        )

        other = FockSpace(:other)
        other_sum = Σ(
            IndexedOperator(Destroy(other, :b), Index(other, :k, 10, other)),
            Index(other, :k, 10, other)
        )
        @test length(other_sum * a) == 1

        @test_throws ArgumentError Σ(ai, i) * Σ(adjoint(ai), i)
    end

    @testset "sum scope survives public arithmetic" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        summed = Σ(ai, i)
        @test isequal(adjoint(summed), Σ(ai', i))
        @test iszero(Σ(0, i))
        @test iszero(simplify(Σ(IndexedVariable(:g, j), i) - 10 * IndexedVariable(:g, j)))

        renamed = change_index(Σ(aj, j), j, i)
        @test isequal(renamed, summed)
        @test length(summed + Σ(aj, j)) == 2
        @test isequal(summed + renamed, 2 * summed)

        scoped = Σ(ai, i) + Σ(ai, i, [j])
        @test length(scoped) == 2
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(scoped),
        )
    end

    @testset "public substitution prunes dead constraints" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        indexed = IndexedOperator(a, i)

        scalar = substitute(Σ(indexed, i, [j]), Dict(indexed => 2))
        @test iszero(scalar - 2)
        @test isempty(constraint_pairs(scalar))
    end

    @testset "ordering distinguishes public constraint scopes" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        k = Index(hf, :k, 10, hf)
        indexed = IndexedOperator(a, i)
        constrained = Σ(indexed, i, [j]) + Σ(indexed, i, [k])
        ordered = SecondQuantizedAlgebra.sorted_arguments(constrained)
        @test length(ordered) == 2
        @test Set(constraint_pairs(ordered[1] + ordered[2])) ==
            Set(constraint_pairs(constrained))
    end

    @testset "nested indexed sums retain diagonal semantics" begin
        h = NLevelSpace(:atom, 2, 1)
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)

        double_sum = Σ(Σ(σ(2, 1, i) * σ(1, 2, j), i), j)
        @test haskey(double_sum, operators(Σ(σ(2, 2, j), j)))
        @test any(
            pair -> pair == (i, j) || pair == (j, i),
            constraint_pairs(double_sum),
        )

        independent = Σ(Σ(σ(2, 2, i), i), j)
        @test iszero(simplify(independent - 10 * Σ(σ(2, 2, i), i)))
    end

    @testset "change_index is simultaneous and semantic" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        k = Index(hf, :k, 10, hf)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        gi = IndexedVariable(:g, i)
        gj = IndexedVariable(:g, j)

        @test isequal(change_index(ai, i, j), aj)
        @test isequal(change_index(ai, k, j), ai)
        @test iszero(simplify(change_index(gi * ai, i, j) - gj * aj))
        @test iszero(simplify(change_index(Σ(ai, i, [j]), j, k) - Σ(ai, i, [k])))
        @test change_index(42, i, j) == 42

        two_sites = IndexedOperator(
            Transition(hf ⊗ hn, :σ, 1, 2, 2),
            Index(hf ⊗ hn, :i, 10, hn)
        ) *
            IndexedOperator(
            Transition(hf ⊗ hn, :σ, 2, 2, 2),
            Index(hf ⊗ hn, :j, 10, hn)
        )
        swapped = change_index(
            two_sites, Dict(
                Index(hf ⊗ hn, :i, 10, hn) => Index(hf ⊗ hn, :j, 10, hn),
                Index(hf ⊗ hn, :j, 10, hn) => Index(hf ⊗ hn, :i, 10, hn),
            )
        )
        expected = IndexedOperator(
            Transition(hf ⊗ hn, :σ, 1, 2, 2),
            Index(hf ⊗ hn, :j, 10, hn)
        ) *
            IndexedOperator(
            Transition(hf ⊗ hn, :σ, 2, 2, 2),
            Index(hf ⊗ hn, :i, 10, hn)
        )
        @test iszero(simplify(swapped - expected))
        @test change_index(two_sites, Dict{Index, Index}()) === two_sites
    end

    @testset "index substitution can remove forbidden diagonal terms" begin
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        Ωij = DoubleIndexedVariable(:Ω, i, j; identical = false)
        a = IndexedOperator(Destroy(hf, :a), i)
        @test iszero(change_index(Ωij, j, i))
        @test iszero(change_index(Ωij * a, j, i))
        @test iszero(change_index(Ωij + 3, j, i) - 3)
        @test iszero(change_index(Ωij * a, Dict(j => i)))
    end

    @testset "cross-space and nested sums" begin
        ha = NLevelSpace(:atom_a, (:g, :e))
        hb = NLevelSpace(:atom_b, (:g, :e))
        h = ha ⊗ hb
        i = Index(h, :i, 4, ha)
        j = Index(h, :j, 4, hb)
        σa = IndexedOperator(Transition(h, :σa, :g, :e, 1), i)
        σb = IndexedOperator(Transition(h, :σb, :g, :e, 2), j)
        sum_b = Σ(σb' * σb, j)
        @test iszero(commutator(σa, sum_b))
        @test get_indices(Σ(σa * σb, i, j)) == [i, j]
        @test iszero(simplify(Σ(σa, i, j) - 4 * Σ(σa, i)))
    end

    @testset "public indexed operations are inferable" begin
        a = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        ai = IndexedOperator(a, i)
        @inferred IndexedOperator(a, i)
        @inferred Σ(ai, i)
        @inferred change_index(ai, i, Index(hf, :j, 10, hf))
        @inferred get_indices(ai)
        @inferred commutator(ai, ai')
    end
end
