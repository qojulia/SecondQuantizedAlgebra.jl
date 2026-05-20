using SecondQuantizedAlgebra
using Test
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
import SecondQuantizedAlgebra: QAdd, QSym, Index, simplify, sorted_arguments,
    has_sum_metadata, get_sum_indices

# Regression tests for `.indices` hygiene under +/- and for per-term `average`
# metadata. These exercise the cancellation paths that historically left dead
# summation indices on QAdd.indices and the average() path that historically
# attached SumIndices metadata to every term regardless of whether the term
# actually depended on a bound index.

@testset "Index hygiene under +/-" begin
    hf = FockSpace(:f)
    ha = NLevelSpace(:atom, 2)
    h = hf ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    @variables N
    g(k) = IndexedVariable(:g, k)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)
    iF = Index(hf, :i, N, hf)
    jF = Index(hf, :j, N, hf)

    @testset "Σ minus itself: full cancellation drops the index" begin
        s = Σ(g(i) * σ(1, 2, i), i)
        d = s - s
        @test d isa QAdd
        @test isempty(d.arguments)
        @test isempty(d.indices)
    end

    @testset "Σ + (-Σ): full cancellation drops the index" begin
        s = Σ(g(i) * σ(1, 2, i), i)
        d = s + (-s)
        @test isempty(d.arguments)
        @test isempty(d.indices)
    end

    @testset "Disjoint-index Σ subtraction: indices stay until cancelled" begin
        s_i = Σ(g(i) * σ(1, 2, i), i)
        s_j = Σ(g(j) * σ(1, 2, j), j)
        d = s_i - s_j
        @test d isa QAdd
        # Both summands survive (different bound names), so both indices stay.
        @test Set(d.indices) == Set([i, j])
    end

    @testset "Partial cancellation prunes only the dead index" begin
        # x_i and x_j each independent sums; subtract a clone of x_i so only j
        # has any surviving term that depends on it.
        x_i = Σ(g(i) * σ(1, 2, i), i)
        x_j = Σ(g(j) * σ(1, 2, j), j)
        combined = x_i + x_j
        d = combined - x_i
        @test d isa QAdd
        @test !isempty(d.arguments)
        # The i-bound terms cancelled; only j-bound contributions remain.
        @test j in d.indices
        @test !(i in d.indices)
    end

    @testset "Two-index Σ minus itself: both indices dropped" begin
        # Σ_i Σ_j (...) is built via product of two single sums on Fock space.
        ai = IndexedOperator(a, iF)
        aj_dag = IndexedOperator(a', jF)
        s = Σ(ai, iF) * Σ(aj_dag, jF)
        d = s - s
        @test d isa QAdd
        @test isempty(d.arguments)
        @test isempty(d.indices)
    end

    @testset "Sum-sum commutator [Σ a_i, Σ a†_j] yields a constant" begin
        # The off-diagonal contributions cancel exactly between H*σ and σ*H.
        # Only the c-number commutator residual survives, so no operator
        # depends on i or j and `.indices` must be empty. Per the package's
        # per-term scope convention, the residual is a single c-number (1),
        # not 2 (which is the double-emission bug from bidirectional pinning).
        ai = IndexedOperator(a, iF)
        aj_dag = IndexedOperator(a', jF)
        c = commutator(Σ(ai, iF), Σ(aj_dag, jF))
        @test c isa QAdd
        @test isempty(c.indices)
        @test length(c.arguments) == 1
        only_term = first(keys(c.arguments))
        @test isempty(only_term.ops)
        # Residual: a single [a_k, a_k'] = 1 contribution summed over k = N.
        cnum_residual = sum(coef for (_, coef) in c.arguments; init = Complex(0))
        @test isequal(cnum_residual, Complex(N))
    end

    @testset "Σ * Σ on same Pauli space: diagonal pin not double-emitted" begin
        hP = PauliSpace(:p)
        iP = Index(hP, :ip, N, hP)
        jP = Index(hP, :jp, N, hP)
        σxi = IndexedOperator(Pauli(hP, :σ, 1), iP)
        σyj = IndexedOperator(Pauli(hP, :σ, 2), jP)
        p = Σ(σxi, iP) * Σ(σyj, jP)
        # Expected dict entries:
        #   1. diagonal contribution: i*σz on a single survivor index.
        #   2. off-diagonal σx_i * σy_j with ne = [(i, j)].
        # Before fix: 3 entries (two duplicated diagonals from i↔j pin).
        @test length(p.arguments) == 2
        # Diagonal pieces (no `ne` constraint), exactly one.
        diag_entries = [t for t in keys(p.arguments) if isempty(t.ne)]
        @test length(diag_entries) == 1
        # Off-diagonal piece carries (i, j) inequality, exactly one.
        off_entries = [t for t in keys(p.arguments) if !isempty(t.ne)]
        @test length(off_entries) == 1
    end

    @testset "Σ over expression: constants pick up range factor" begin
        # `Σ_i (a_i + 1)` should mathematically yield `N + Σ_i a_i`. The `+1`
        # term doesn't depend on `i`, but it lives inside the Σ_i scope so
        # summing it over i = 1..N gives N.
        hfL = FockSpace(:fL)
        @qnumbers aL::Destroy(hfL)
        iL = Index(hfL, :iL, N, hfL)
        aiL = IndexedOperator(aL, iL)
        r = Σ(aiL + 1, iL)
        cnums = [coef for (t, coef) in r.arguments if isempty(t.ops)]
        @test length(cnums) == 1
        @test isequal(cnums[1], Complex(N))
    end

    @testset "Σ_i (a_i a'_i) yields N + Σ_i a'_i a_i" begin
        # The eager a*a' = a'a + 1 canonicalization produces a constant `+1`
        # residual; when wrapped in `Σ_i`, the constant must pick up N.
        hfL = FockSpace(:fM)
        @qnumbers aL::Destroy(hfL)
        iL = Index(hfL, :iM, N, hfL)
        aiL = IndexedOperator(aL, iL)
        aiL_dag = IndexedOperator(aL', iL)
        r = Σ(aiL * aiL_dag, iL)
        cnums = [coef for (t, coef) in r.arguments if isempty(t.ops)]
        @test length(cnums) == 1
        @test isequal(cnums[1], Complex(N))
    end

    @testset "[Σ a_i, Σ a'_j] = N (canonical commutator)" begin
        # The diagonal pin in the product emits the +1 residual under the
        # surviving bound scope, so it must carry an N factor. Currently
        # returns the un-scaled constant 1, which is mathematically wrong.
        c = commutator(
            Σ(IndexedOperator(a, iF), iF),
            Σ(IndexedOperator(a', jF), jF)
        )
        cnums = [coef for (t, coef) in c.arguments if isempty(t.ops)]
        @test length(cnums) == 1
        @test isequal(cnums[1], Complex(N))
    end

    @testset "QAdd + QSym cancellation drops the dead index" begin
        # Build Σ(σ_i, i) + (-Σ(σ_j, j))-equivalent via QSym path: the QAdd
        # term cancels against the added QSym only when shapes match. Here we
        # build a single-term QAdd that cancels a QSym addition.
        op_j = σ(1, 2, j)
        wrapped = Σ(σ(1, 2, i), i) - op_j + op_j  # adding op_j then subtracting cancels
        @test wrapped isa QAdd
        # i is the only surviving bound index; the σ_j cancellation must not
        # have stranded anything on `.indices`.
        @test i in wrapped.indices
        for idx in wrapped.indices
            @test idx == i
        end
    end
end

@testset "average() per-term Σ metadata" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    @variables N Δ
    g(k) = IndexedVariable(:g, k)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    @testset "average() on QAdd with only sum-independent terms drops Σ" begin
        # commutator(H, σ_j) is a sum of two j-only terms even though H has a
        # bound index i. The averaged result must NOT wrap each term in Σ_i.
        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
        c = commutator(H, σ(1, 2, j))
        # The QAdd we average has no surviving dependence on i.
        @test isempty(c.indices)
        avg = average(1im * c)
        # No subterm should carry SumIndices metadata.
        for arg in SymbolicUtils.arguments(SymbolicUtils.unwrap(avg))
            @test !has_sum_metadata(arg)
        end
    end

    @testset "average() splits metadata: indep terms bare, dep terms wrapped" begin
        # Σ(g(i)*σ_i, i) * σ_j produces a diagonal (j-only) term and an
        # off-diagonal (i≠j) sum term in the same QAdd. After average(), only
        # the dep term should carry SumIndices metadata.
        prod = Σ(g(i) * σ(1, 2, i), i) * σ(2, 1, j)
        @test i in prod.indices
        avg = average(prod)
        unwrapped = SymbolicUtils.unwrap(avg)
        dep_count = 0
        indep_count = 0
        # average returns a sum of (cnum * avg-call) or bare avg-call; the avg
        # nodes carry metadata. Walk the argument tree shallowly.
        function carries_sum(x)
            x isa SymbolicUtils.BasicSymbolic || return false
            has_sum_metadata(x) && return true
            SymbolicUtils.iscall(x) || return false
            return any(carries_sum, SymbolicUtils.arguments(x))
        end
        for arg in SymbolicUtils.arguments(unwrapped)
            if carries_sum(arg)
                dep_count += 1
            else
                indep_count += 1
            end
        end
        @test dep_count >= 1   # the (i≠j) off-diagonal piece
        @test indep_count >= 1 # the pinned diagonal piece
    end

    @testset "average() of c-number (empty-ops) term carries no Σ metadata" begin
        # If a QAdd has Σ scope but a term is purely a c-number (e.g. arising
        # from a commutator residual), the averaged contribution is a bare
        # scalar and must not be wrapped in Σ.
        hf = FockSpace(:f)
        @qnumbers b::Destroy(hf)
        iF = Index(hf, :i, N, hf)
        jF = Index(hf, :j, N, hf)
        bi = IndexedOperator(b, iF)
        bj_dag = IndexedOperator(b', jF)
        c = commutator(Σ(bi, iF), Σ(bj_dag, jF))
        avg = average(c)
        # The result is a plain number / Num, with no Σ metadata anywhere.
        unwrapped = SymbolicUtils.unwrap(avg)
        @test !(unwrapped isa SymbolicUtils.BasicSymbolic) || !has_sum_metadata(unwrapped)
    end
end

@testset "Superradiant-laser commutators (QC reference form)" begin
    # Cross-check against the form reported in QuantumCumulants.jl docs:
    # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, kk) = IndexedOperator(Transition(h, :σ, α, β, 2), kk)
    @variables N Δ
    g(kk) = IndexedVariable(:g, kk)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    @testset "[H, σ_j^{12}] reduces to two j-only terms" begin
        c = commutator(H, σ(1, 2, j))
        # Closed form: i is pinned to j on both branches of the commutator,
        # and the off-diagonal (i≠j) contributions cancel exactly.
        expected = g(j) * a * σ(2, 2, j) - g(j) * a * σ(1, 1, j)
        @test c == expected
        @test isempty(c.indices)
        @test length(c.arguments) == 2
    end

    @testset "[H, a' σ_j^{12}] keeps the legitimate i≠j off-diagonal Σ" begin
        c = commutator(H, a' * σ(1, 2, j))
        # Diagonal contributions (i pinned to j) and the (i≠j) off-diagonal
        # built independently. The off-diagonal piece is the same QAdd as
        # what the diagonal-split machinery produces from
        # Σ_i g(i) σ_i^{21} σ_j^{12}, minus its own pinned i=j piece.
        full = Σ(g(i) * σ(2, 1, i) * σ(1, 2, j), i)
        off_diag = full - g(j) * σ(2, 2, j)
        expected = g(j) * σ(2, 2, j) - Δ * a' * σ(1, 2, j) +
            g(j) * a' * a * σ(2, 2, j) - g(j) * a' * a * σ(1, 1, j) +
            off_diag
        @test c == expected
        @test c.indices == [i]
        @test length(c.arguments) == 5
        # Exactly one surviving term carries the (i,j) inequality.
        ne_terms = [t for t in keys(c.arguments) if !isempty(t.ne)]
        @test length(ne_terms) == 1
    end

    @testset "[H, σ_j^{12} σ_k^{21}] matches Leibniz under j≠k" begin
        c = commutator(H, σ(1, 2, j) * σ(2, 1, k))
        @test isempty(c.indices)
        # Without j≠k the direct form retains non-canonical j/k interleavings
        # that don't reduce (free indices stay Undetermined). Under j≠k both
        # the direct and Leibniz expansions canonicalize to the same form.
        leibniz = commutator(H, σ(1, 2, j)) * σ(2, 1, k) +
            σ(1, 2, j) * commutator(H, σ(2, 1, k))
        @test assume_distinct_index(c, [(j, k)]) ==
            assume_distinct_index(leibniz, [(j, k)])
    end
end
