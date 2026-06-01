using SecondQuantizedAlgebra
using Test
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
import SecondQuantizedAlgebra: QAdd, QSym, QField, Index, simplify, normal_order,
    _depends_on_index_term, _any_depends_on_index, _iszero_cnum,
    has_sum_metadata, get_sum_indices, is_average, undo_average

# ============================================================================
# Invariant helpers
#
# These two predicates are the structural contracts that QAdd and average()
# outputs must satisfy. Bugs A (stale `.indices` after cancellation) and B
# (blanket `SumIndices` metadata on every averaged term) both violated one
# of these but slipped through the example-based suite because every
# isolated primitive looked correct on hand-rolled inputs.
# ============================================================================

"""
Every bound index advertised in `q.indices` must be carried by at least one
surviving term, and every dict entry must have a nonzero coefficient.
A violation means downstream code (printing, average) will wrap dead
indices around terms that no longer depend on them.
"""
function check_qadd_invariants(q::QAdd)
    failures = String[]
    for idx in q.indices
        _any_depends_on_index(q, idx) ||
            push!(failures, "index $(idx.name) in .indices has no surviving carrier")
    end
    for (term, c) in q.arguments
        _iszero_cnum(c) && push!(failures, "term $(term.ops) has zero prefactor")
    end
    return failures
end

"""
Walk an averaged symbolic expression and collect every AvgFunc node. For
each, if it carries `SumIndices` metadata then the underlying QField must
actually depend on at least one of those indices. Blanket metadata (Bug B)
violates this.
"""
function check_average_metadata(avg)
    failures = String[]
    _walk_avg!(failures, SymbolicUtils.unwrap(avg))
    return failures
end

function _walk_avg!(failures, x)
    x isa SymbolicUtils.BasicSymbolic || return
    if is_average(x)
        if has_sum_metadata(x)
            sum_idxs = get_sum_indices(x)
            inner = undo_average(x)
            inner isa QAdd || return
            any_dep = false
            for idx in sum_idxs
                if _any_depends_on_index(inner, idx)
                    any_dep = true
                    break
                end
            end
            any_dep || push!(
                failures,
                "average term $(inner) carries SumIndices=$(sum_idxs) but " *
                    "does not depend on any of them"
            )
        end
        return
    end
    SymbolicUtils.iscall(x) || return
    for arg in SymbolicUtils.arguments(x)
        _walk_avg!(failures, arg)
    end
    return
end

assert_qadd_ok(q::QAdd, ctx::AbstractString) = begin
    fs = check_qadd_invariants(q)
    @test isempty(fs) || (println("QAdd invariant violated [$ctx]: ", fs); false)
end

assert_avg_ok(avg, ctx::AbstractString) = begin
    fs = check_average_metadata(avg)
    @test isempty(fs) || (println("avg invariant violated [$ctx]: ", fs); false)
end

# ============================================================================
# Directed-fuzz corpus
#
# Each entry is a (label, ()->QAdd) pair. The harness applies the QAdd
# invariant to every entry, then averages the entry and applies the
# metadata invariant. The corpus exercises the composition paths that
# bugs A and B lived on: sum-product, sum-sum, commutator-of-sums,
# cancellation, and double-pinning. New paths can be appended without
# rewriting the runner.
# ============================================================================

@testset "Invariants on QAdd and average() outputs" begin

    # Operator setup shared by all corpus entries.
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    @variables N Δ κ
    g(k) = IndexedVariable(:g, k)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    # Pure-Fock setup for sum-sum commutator and Spin/Pauli laws.
    hf = FockSpace(:f)
    @qnumbers b::Destroy(hf)
    iF = Index(hf, :i, N, hf)
    jF = Index(hf, :j, N, hf)
    bi = IndexedOperator(b, iF)
    bj = IndexedOperator(b, jF)
    bj_dag = IndexedOperator(b', jF)

    H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    corpus = [
        ("Σ minus itself (Bug A)", () -> Σ(g(i) * σ(1, 2, i), i) - Σ(g(i) * σ(1, 2, i), i)),
        ("Σ product diag split", () -> Σ(g(i) * σ(1, 2, i), i) * σ(2, 1, j)),
        ("Σ * external op other space", () -> Σ(g(i) * σ(1, 2, i), i) * a),
        ("commutator(H, σ_j)", () -> commutator(H, σ(1, 2, j))),
        ("commutator(H, a' σ_j)", () -> commutator(H, a' * σ(1, 2, j))),
        ("commutator(H, σ_j σ_k)", () -> commutator(H, σ(1, 2, j) * σ(2, 1, k))),
        ("sum-sum commutator (Bug A)", () -> commutator(Σ(bi, iF), Σ(bj_dag, jF))),
        (
            "(Σ_i b_i)(Σ_j b†_j) cancels", () -> Σ(bi, iF) * Σ(bj_dag, jF) -
                Σ(bi, iF) * Σ(bj_dag, jF),
        ),
        (
            "nested commutator on sums", () -> begin
                # Inner commutator binds j; rename one side to avoid bound-name clash.
                kF = Index(hf, :k, N, hf)
                bk = IndexedOperator(b, kF)
                commutator(Σ(bi, iF), commutator(Σ(bj_dag, jF), Σ(bk, kF)))
            end,
        ),
        ("simplify(H - H)", () -> simplify(H - H)),
        ("normal_order on sum prod", () -> normal_order(Σ(g(i) * σ(2, 2, i), i) * a' * a)),
        (
            "partial cancel: i dies, j lives",
            () -> (Σ(g(i) * σ(1, 2, i), i) + Σ(g(j) * σ(2, 1, j), j)) - Σ(g(i) * σ(1, 2, i), i),
        ),
        ("double pin (i pinned twice)", () -> Σ(g(i) * σ(1, 2, i) * σ(2, 1, i), i) * σ(2, 1, j) * σ(1, 2, k)),
        ("anticommutator on sums", () -> anticommutator(Σ(bi, iF), Σ(bj_dag, jF))),
        # ----- Operations the original corpus did not exercise -----
        ("expand_completeness on Σ", () -> expand_completeness(Σ(σ(1, 1, i) * a, i))),
        ("adjoint of Σ-product", () -> (Σ(g(i) * σ(1, 2, i), i) * σ(2, 1, j))'),
        (
            "substitute on Σ", () -> substitute(
                Σ(g(i) * σ(1, 2, i), i),
                Dict(σ(1, 2, j) => 2 * σ(1, 2, j))
            ),
        ),
        (
            "change_index after cancel", () -> change_index(
                (Σ(g(i) * σ(1, 2, i), i) + Σ(g(j) * σ(2, 1, j), j)) - Σ(g(i) * σ(1, 2, i), i),
                j, k
            ),
        ),
        (
            "Pauli sum-sum commutator", () -> begin
                hP = PauliSpace(:p)
                iP = Index(hP, :ip, N, hP)
                jP = Index(hP, :jp, N, hP)
                σxi = IndexedOperator(Pauli(hP, :σ, 1), iP)
                σyj = IndexedOperator(Pauli(hP, :σ, 2), jP)
                commutator(Σ(σxi, iP), Σ(σyj, jP))
            end,
        ),
        (
            "Dicke [H_dicke, σy_j]", () -> begin
                # Cavity ⊗ Pauli-atom, atoms collectively coupled. Tests Pauli
                # commutator emission interacting with both an indexed sum over
                # atoms and a Fock-side prefactor (a + a').
                hcD = FockSpace(:cavity)
                hpD = PauliSpace(:atom)
                hD = hcD ⊗ hpD
                aD = Destroy(hD, :a, 1)
                iD = Index(hD, :id, N, hpD)
                jD = Index(hD, :jd, N, hpD)
                σxD(kk) = IndexedOperator(Pauli(hD, :σ, 1, 2), kk)
                σyD(kk) = IndexedOperator(Pauli(hD, :σ, 2, 2), kk)
                σzD(kk) = IndexedOperator(Pauli(hD, :σ, 3, 2), kk)
                @variables ω₀ ωₐ λ
                H_dicke = ω₀ * aD' * aD + ωₐ * Σ(σzD(iD), iD) +
                    λ * (aD + aD') * Σ(σxD(iD), iD)
                commutator(H_dicke, σyD(jD))
            end
        ),
    ]

    @testset "$label" for (label, builder) in corpus
        q = builder()
        @test q isa QAdd
        assert_qadd_ok(q, label)
        avg = average(q)
        assert_avg_ok(avg, label)
    end
end

# ============================================================================
# Algebraic laws on the operator zoo
#
# Property tests, not example tests. Each law holds for every (A, B) pair
# we construct, so adding a new operator type or composition shape only
# requires extending the basis; the assertions need no edits.
# ============================================================================

@testset "Algebraic laws" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, kk) = IndexedOperator(Transition(h, :σ, α, β, 2), kk)
    @variables N
    g(kk) = IndexedVariable(:g, kk)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    hs = SpinSpace(:s)
    Sx = Spin(hs, :S, 1)
    Sy = Spin(hs, :S, 2)
    Sz = Spin(hs, :S, 3)

    # Build a small operator basis spanning leaf / product / sum / sum-prod
    # across Fock, NLevel, and Spin so the law assertions hit every code
    # path the algebra uses.
    basis = [
        ("a", a),
        ("a'", a'),
        ("a'*a", a' * a),
        ("a + a'", a + a'),
        ("σ_j^{12}", σ(1, 2, j)),
        ("σ_j^{11}", σ(1, 1, j)),
        ("a * σ_j^{12}", a * σ(1, 2, j)),
        ("Σ_i g(i) σ_i^{12}", Σ(g(i) * σ(1, 2, i), i)),
        ("Σ_i g(i) a σ_i^{21}", Σ(g(i) * a * σ(2, 1, i), i)),
        ("Sx", Sx),
        ("Sx + Sy", Sx + Sy),
        ("Sx * Sy", Sx * Sy),
    ]

    @testset "commutator antisymmetry [A,B] + [B,A] = 0" begin
        for (la, A) in basis, (lb, B) in basis
            # Skip pairs that share a bound index name (product would throw).
            A isa QAdd && B isa QAdd &&
                !isempty(A.indices) && !isempty(B.indices) &&
                any(idx -> idx in B.indices, A.indices) && continue
            sum = commutator(A, B) + commutator(B, A)
            @test iszero(simplify(sum)) ||
                (println("antisymmetry failed for ($la, $lb): ", sum); false)
        end
    end

    @testset "commutator bilinearity [A, B+C] = [A,B] + [A,C]" begin
        # Pick a representative spread; the full Cartesian cube is overkill.
        A_cases = [a, σ(1, 2, j), Sx]
        B_cases = [a', σ(2, 1, j), Sy]
        C_cases = [a, σ(1, 1, j), Sz]
        for (A, B, C) in zip(A_cases, B_cases, C_cases)
            lhs = commutator(A, B + C)
            rhs = commutator(A, B) + commutator(A, C)
            @test iszero(simplify(lhs - rhs))
        end
    end

    @testset "average linearity average(A + B) = average(A) + average(B)" begin
        for (la, A) in basis, (lb, B) in basis
            A isa QAdd && B isa QAdd &&
                !isempty(A.indices) && !isempty(B.indices) &&
                any(idx -> idx in B.indices, A.indices) && continue
            lhs = average(A + B)
            rhs = average(A) + average(B)
            # Both sides go through `simplify` on the underlying QAdd so
            # symbolic prefactors match up.
            d = simplify(undo_average(lhs) - undo_average(rhs))
            @test iszero(d) ||
                (println("linearity failed for ($la, $lb)"); false)
        end
    end

    @testset "undo_average ∘ average roundtrip" begin
        for (label, A) in basis
            qa = A isa QAdd ? A : (1 * A)
            r = undo_average(average(qa))
            @test iszero(simplify(r - qa)) ||
                (println("roundtrip failed for $label: r=$r, qa=$qa"); false)
        end
    end

    @testset "Spin Jacobi identity [[Sx,Sy],Sz] + [[Sy,Sz],Sx] + [[Sz,Sx],Sy] = 0" begin
        j1 = commutator(commutator(Sx, Sy), Sz)
        j2 = commutator(commutator(Sy, Sz), Sx)
        j3 = commutator(commutator(Sz, Sx), Sy)
        @test iszero(simplify(j1 + j2 + j3))
    end
end
