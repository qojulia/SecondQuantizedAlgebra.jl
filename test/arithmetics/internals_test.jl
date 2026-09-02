using Test
using SecondQuantizedAlgebra
using Symbolics: @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, Transition, NO_INDEX,
    partial_sort!, site_compare, can_commute, commute_pair,
    reduce_pair, may_reduce, ground_state_expand,
    SiteCmp, Less, Greater, Equal, Undetermined,
    ReduceKind, NoReduction, ScalarReduction, OpReduction,
    CNUM_ONE, CNUM_ZERO, CNUM_NEG1, CNUM_IM, CNUM_NEG_IM, EMPTY_NE,
    mul_cnum, add_cnum, neg_cnum, to_cnum, NUM_ZERO, NUM_ONE,
    canonicalize_to_dict!, QTermDict,
    reduce_ops, commute_ops, expand_gs_ops, substitute_ops,
    stream!, canonicalize!, Index

# Tests of private hooks, passes, sort, and canonical-form invariants.

# Helper: is this op a ground-state projector of an NLevelSpace?
function is_gs_projector(op)
    is_transition(op) || return false
    return op.l1 == op.g && op.l2 == op.g
end

# Helper: every dict key in `expr` is free of ground-state projectors
function no_gs_projectors(expr::QAdd)
    for term in keys(expr.arguments)
        any(is_gs_projector, term.ops) && return false
    end
    return true
end

# Helper: a term is canonical when no neighboring pair could swap or reduce.
function is_canonical(t)
    ops = t.ops
    ne = t.ne
    for i in 1:(length(ops) - 1)
        cmp = site_compare(ops[i], ops[i + 1], ne)
        if cmp === Greater
            return false
        end
        if cmp === Equal
            can_commute(ops[i], ops[i + 1]) || return false
            first(reduce_pair(ops[i], ops[i + 1])) === NoReduction || return false
        end
    end
    return true
end

@testset "Internals" begin

    @testset "Operator hooks" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = adjoint(a)
        @test site_compare(a, ad, EMPTY_NE) === Equal
        @test site_compare(ad, a, EMPTY_NE) === Equal
        @test can_commute(a, ad) === false
        @test can_commute(ad, a) === true
        sw = commute_pair(a, ad)
        @test sw[1] === ad && sw[2] === a

        ha = NLevelSpace(:atom, 2)
        σ12 = Transition(ha, :σ, 1, 2)
        σ21 = Transition(ha, :σ, 2, 1)
        @test site_compare(σ12, σ21, EMPTY_NE) === Equal
        @test can_commute(σ12, σ21) === false
        red = reduce_pair(σ21, σ12)   # σ²¹·σ¹² = σ²²
        @test red[1] === OpReduction
        @test red[2].l1 == 2 && red[2].l2 == 2
        @test red[3] === CNUM_ONE
        @test first(ground_state_expand(Transition(ha, :σ, 1, 1)))
    end

    @testset "Partial sort" begin
        @testset "distinct sites" begin
            h = FockSpace(:c) ⊗ NLevelSpace(:a, 2)
            a = Destroy(h, :a, 1)
            σ = Transition(h, :σ, 1, 2)
            ops = Op[σ, a]
            partial_sort!(ops, EMPTY_NE)
            @test is_destroy(ops[1])
            @test is_transition(ops[2])
        end

        @testset "same-site preserved" begin
            h = NLevelSpace(:a, 3)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            ops = Op[σ12, σ21]
            pre = copy(ops)
            partial_sort!(ops, EMPTY_NE)
            @test ops == pre
        end

        @testset "undetermined preserved" begin
            ha = NLevelSpace(:a, 2)
            i = Index(ha, :i, 5, ha)
            j = Index(ha, :j, 5, ha)
            σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
            σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
            ops = Op[σ_j, σ_i]
            partial_sort!(ops, EMPTY_NE)
            @test ops[1] === σ_j
            @test ops[2] === σ_i
        end

        @testset "ne resolves undetermined to distinct" begin
            ha = NLevelSpace(:a, 2)
            i = Index(ha, :i, 5, ha)
            j = Index(ha, :j, 5, ha)
            σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
            σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
            ops = Op[σ_j, σ_i]
            partial_sort!(ops, [(i, j)])
            @test ops[1] === σ_i
            @test ops[2] === σ_j
        end
    end

    @testset "Pass primitives" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = adjoint(a)
        hn = NLevelSpace(:a, 3)

        @testset "canonicalize_to_dict!: basic insert" begin
            out = QTermDict()
            canonicalize_to_dict!(out, Op[a], CNUM_ONE, EMPTY_NE)
            @test length(out) == 1
        end

        @testset "canonicalize_to_dict!: like-term collection" begin
            out = QTermDict()
            canonicalize_to_dict!(out, Op[a], CNUM_ONE, EMPTY_NE)
            canonicalize_to_dict!(out, Op[a], CNUM_ONE, EMPTY_NE)
            @test length(out) == 1
            @test first(values(out)) == 2 + 0im
        end

        @testset "canonicalize_to_dict!: zero coeff dropped" begin
            out = QTermDict()
            canonicalize_to_dict!(out, Op[a], CNUM_ZERO, EMPTY_NE)
            @test isempty(out)
        end

        @testset "canonicalize_to_dict!: symbolic prefactor cancellation" begin
            @variables γ D
            # γ/D + (-γ)/D is structurally non-zero to Symbolics (the `Div` nodes
            # are left un-combined), so exact-negation collection must drop the term.
            out = QTermDict()
            canonicalize_to_dict!(out, Op[a], to_cnum(γ / D), EMPTY_NE)
            canonicalize_to_dict!(out, Op[a], to_cnum(-γ / D), EMPTY_NE)
            @test isempty(out)
        end

        @testset "canonicalize_to_dict!: distinct symbolic prefactors kept" begin
            @variables γ D β
            out = QTermDict()
            canonicalize_to_dict!(out, Op[a], to_cnum(γ / D), EMPTY_NE)
            canonicalize_to_dict!(out, Op[a], to_cnum(β / D), EMPTY_NE)
            @test length(out) == 1
        end

        @testset "reduce_ops: Transition composition" begin
            σ12 = Transition(hn, :σ, 1, 2)
            σ23 = Transition(hn, :σ, 2, 3)
            emitted = Tuple{Vector{Op}, CNum}[]
            reduce_ops(Op[σ12, σ23], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test length(emitted[1][1]) == 1
            @test emitted[1][1][1].l1 == 1 && emitted[1][1][1].l2 == 3
        end

        @testset "reduce_ops: zero from incompatible composition" begin
            σ12 = Transition(hn, :σ, 1, 2)
            σ31 = Transition(hn, :σ, 3, 1)
            emitted = Tuple{Vector{Op}, CNum}[]
            reduce_ops(Op[σ12, σ31], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test isempty(emitted)
        end

        @testset "reduce_ops: no-op input passes through" begin
            emitted = Tuple{Vector{Op}, CNum}[]
            reduce_ops(Op[ad, a], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == Op[ad, a]
        end

        @testset "commute_ops: Fock aa† → a†a + 1" begin
            emitted = Tuple{Vector{Op}, CNum}[]
            commute_ops(Op[a, ad], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 2
            sort!(emitted, by = e -> length(e[1]))
            @test isempty(emitted[1][1]) && emitted[1][2] == CNUM_ONE
            @test emitted[2][1] == Op[ad, a]
        end

        @testset "commute_ops: no-op on already-ordered pair" begin
            emitted = Tuple{Vector{Op}, CNum}[]
            commute_ops(Op[ad, a], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == Op[ad, a]
        end

        @testset "expand_gs_ops: σ¹¹ → 1 - σ²²" begin
            h2 = NLevelSpace(:a, 2)
            σ11 = Transition(h2, :σ, 1, 1)
            emitted = Tuple{Vector{Op}, CNum}[]
            expand_gs_ops(Op[σ11], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 2
            sort!(emitted, by = e -> length(e[1]))
            @test isempty(emitted[1][1]) && emitted[1][2] == CNUM_ONE
            @test length(emitted[2][1]) == 1 && emitted[2][2] == -CNUM_ONE
            @test emitted[2][1][1].l1 == 2 && emitted[2][1][1].l2 == 2
        end

        @testset "expand_gs_ops: passthrough when no σᵍᵍ" begin
            emitted = Tuple{Vector{Op}, CNum}[]
            expand_gs_ops(Op[a], CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == Op[a]
        end

        @testset "substitute_ops: operator → scalar" begin
            emitted = Tuple{Vector{Op}, CNum}[]
            substitute_ops(Op[a, ad], CNUM_ONE, Dict(a => 2)) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][2] == 2 + 0im
            @test emitted[1][1] == Op[ad]
        end

        @testset "stream!: idempotent on canonical input" begin
            out = QTermDict()
            stream!(out, Op[ad, a], CNUM_ONE, EMPTY_NE)
            @test length(out) == 1
        end

        @testset "canonicalize!: aa† → a†a + 1" begin
            out = QTermDict()
            canonicalize!(out, Op[a, ad], CNUM_ONE, EMPTY_NE)
            @test length(out) == 2
        end
    end

    @testset "Canonical-form invariants" begin
        @testset "after `*`" begin
            h = FockSpace(:f) ⊗ NLevelSpace(:a, 3)
            a = Destroy(h, :a, 1)
            σ12 = Transition(h, :σ, 1, 2, 2)
            σ21 = Transition(h, :σ, 2, 1, 2)
            q = a * adjoint(a) * σ12 * σ21
            for (t, _) in q
                @test is_canonical(t)
            end
        end

        @testset "after normal_order" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = a * adjoint(a) * a
            for (t, _) in normal_order(q)
                @test is_canonical(t)
            end
        end

        @testset "after commutator" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = commutator(a, adjoint(a))
            for (t, _) in q
                @test is_canonical(t)
            end
        end

        @testset "after substitute" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = adjoint(a) * a
            sub = SecondQuantizedAlgebra.SymbolicUtils.substitute(q, Dict{Any, Any}(a => 0.5))
            for (t, _) in sub
                @test is_canonical(t)
            end
        end

        @testset "after expand_completeness" begin
            h = NLevelSpace(:a, 3)
            σ11 = Transition(h, :σ, 1, 1)
            exp = expand_completeness(σ11)
            for (t, _) in exp
                @test is_canonical(t)
            end
        end
    end

    @testset "CNum arithmetic type stability" begin
        @test @inferred(mul_cnum(CNUM_ONE, CNUM_ONE)) isa CNum
        @test @inferred(mul_cnum(CNUM_ONE, CNUM_ZERO)) isa CNum
        @test @inferred(mul_cnum(CNUM_NEG1, CNUM_IM)) isa CNum
        @test @inferred(add_cnum(CNUM_ONE, CNUM_ONE)) isa CNum
        @test @inferred(add_cnum(CNUM_ONE, CNUM_NEG1)) isa CNum
        @test @inferred(neg_cnum(CNUM_ONE)) isa CNum
        @test @inferred(neg_cnum(CNUM_IM)) isa CNum

        @variables x y
        c1 = to_cnum(Complex(SecondQuantizedAlgebra.Num(x), NUM_ZERO))
        c2 = to_cnum(Complex(SecondQuantizedAlgebra.Num(y), NUM_ZERO))
        @test @inferred(mul_cnum(c1, c2)) isa CNum
        @test @inferred(add_cnum(c1, c2)) isa CNum
        @test @inferred(neg_cnum(c1)) isa CNum
    end

    @testset "to_cnum type stability" begin
        @test @inferred(to_cnum(1)) isa CNum
        @test @inferred(to_cnum(0)) isa CNum
        @test @inferred(to_cnum(-1)) isa CNum
        @test @inferred(to_cnum(1.5)) isa CNum
        @test @inferred(to_cnum(im)) isa CNum
        @test @inferred(to_cnum(-im)) isa CNum
        @test @inferred(to_cnum(1 + 2im)) isa CNum
    end

    @testset "Operator hook type stability" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = Create(h, :a)
        @test @inferred(site_compare(a, ad, EMPTY_NE)) isa SiteCmp
        @test @inferred(can_commute(a, ad)) isa Bool
        @test @inferred(can_commute(ad, a)) isa Bool
        @test @inferred(commute_pair(a, ad)) isa Tuple{Op, Op, CNum, Vector{Op}, CNum, Vector{Op}}

        ha = NLevelSpace(:atom, 3)
        σ12 = Transition(ha, :σ, 1, 2)
        σ21 = Transition(ha, :σ, 2, 1)
        @test @inferred(reduce_pair(σ12, σ21)) isa Tuple{ReduceKind, Op, CNum}
        @test @inferred(ground_state_expand(σ12)) isa Tuple{Bool, Int, Int, Int}

        hp = PauliSpace(:s)
        px = Pauli(hp, :σ, 1)
        py = Pauli(hp, :σ, 2)
        @test @inferred(reduce_pair(px, py)) isa Tuple{ReduceKind, Op, CNum}

        hs = SpinSpace(:s)
        sx = Spin(hs, :S, 1)
        sy = Spin(hs, :S, 2)
        @test @inferred(commute_pair(sy, sx)) isa Tuple{Op, Op, CNum, Vector{Op}, CNum, Vector{Op}}

        hc = CollectiveNLevelSpace(:collective, 2)
        S12 = CollectiveTransition(hc, :S, 1, 2)
        S21 = CollectiveTransition(hc, :S, 2, 1)
        @test @inferred(commute_pair(S12, S21)) isa Tuple{Op, Op, CNum, Vector{Op}, CNum, Vector{Op}}

        # may_reduce: true iff the same-type pair's reduce_pair can fire
        @test @inferred(may_reduce(σ12, σ21)) === true
        @test @inferred(may_reduce(px, py)) === true
        @test may_reduce(a, ad) === false
        @test may_reduce(sx, sy) === false
        @test may_reduce(S12, S21) === false
        @test may_reduce(a, σ12) === false
    end

    @testset "Canonical form (NLevelSpace ground-state projection)" begin
        @testset "Transition carries ground_state and n_levels" begin
            h = NLevelSpace(:atom, 3, 2)
            σ = Transition(h, :σ, 1, 3)
            @test σ.g == 2
            @test σ.nlev == 3

            h2 = NLevelSpace(:atom, 2)
            σ2 = Transition(h2, :σ, 1, 2)
            @test σ2.g == 1
            @test σ2.nlev == 2

            hf = FockSpace(:c)
            hp = hf ⊗ h
            σp = Transition(hp, :σ, 1, 3, 2)
            @test σp.g == 2
            @test σp.nlev == 3
        end

        @testset "Field preservation through ops" begin
            h = NLevelSpace(:atom, 4, 3)
            σ = Transition(h, :σ, 1, 2)

            σadj = adjoint(σ)
            @test σadj.g == 3
            @test σadj.nlev == 4

            i = Index(h, :i, 10, h)
            σi = IndexedOperator(σ, i)
            @test σi.g == 3
            @test σi.nlev == 4

            j = Index(h, :j, 10, h)
            σj = change_index(σi, i, j)
            @test σj.g == 3
            @test σj.nlev == 4
        end

        @testset "Equality distinguishes ground_state and n_levels" begin
            for (h1, h2) in [
                    (NLevelSpace(:atom, 3, 1), NLevelSpace(:atom, 3, 2)),
                    (NLevelSpace(:atom, 2, 1), NLevelSpace(:atom, 3, 1)),
                ]
                σ1 = Transition(h1, :σ, 1, 2)
                σ2 = Transition(h2, :σ, 1, 2)
                @test σ1 != σ2
                @test hash(σ1) != hash(σ2)
            end
        end

        @testset "Composition+expand: σ¹²·σ²¹ → 1 - σ²² (2-level, g=1)" begin
            h = NLevelSpace(:atom, 2, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ22 = Transition(h, :σ, 2, 2)
            result = expand_completeness(σ12 * σ21)
            @test result isa QAdd
            @test isequal(result, simplify(1 - σ22))
            @test no_gs_projectors(result)
        end

        @testset "Composition+expand: σ¹²·σ²¹ → 1 - σ²² - σ³³ (3-level, g=1)" begin
            h = NLevelSpace(:atom, 3, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ22 = Transition(h, :σ, 2, 2)
            σ33 = Transition(h, :σ, 3, 3)
            result = expand_completeness(σ12 * σ21)
            @test isequal(result, simplify(1 - σ22 - σ33))
            @test no_gs_projectors(result)
        end

        @testset "Composition+expand: ground state ≠ 1 (g=2)" begin
            h = NLevelSpace(:atom, 3, 2)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ11 = Transition(h, :σ, 1, 1)
            σ33 = Transition(h, :σ, 3, 3)

            result = expand_completeness(σ21 * σ12)
            @test isequal(result, simplify(1 - σ11 - σ33))
            @test no_gs_projectors(result)

            result2 = σ12 * σ21
            @test length(result2) == 1
            (term, c) = only(collect(result2))
            ops = term.ops
            @test length(ops) == 1
            @test ops[1] == σ11
        end

        @testset "Composition+expand: indexed σᵢ¹²·σᵢ²¹ preserves index" begin
            ha = NLevelSpace(:atom, 2, 1)
            hf = FockSpace(:c)
            h = hf ⊗ ha

            @variables N
            i = Index(h, :i, N, ha)
            σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

            result = expand_completeness(σ(1, 2) * σ(2, 1))
            @test result isa QAdd
            @test length(result) == 2
            @test no_gs_projectors(result)
            for term in keys(result.arguments)
                for op in term.ops
                    if is_transition(op)
                        @test op.index == i
                        @test op.g == 1
                        @test op.nlev == 2
                    end
                end
            end
        end

        @testset "Composition+expand: σᵍᵍ via longer chain σ¹²·σ²³·σ³¹" begin
            h = NLevelSpace(:atom, 3, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ23 = Transition(h, :σ, 2, 3)
            σ31 = Transition(h, :σ, 3, 1)
            σ22 = Transition(h, :σ, 2, 2)
            σ33 = Transition(h, :σ, 3, 3)

            result = expand_completeness(σ12 * σ23 * σ31)
            @test isequal(result, simplify(1 - σ22 - σ33))
            @test no_gs_projectors(result)
        end

        @testset "Different-site σᵍᵍ_i · σᵍᵍ_j (expand via expand_completeness)" begin
            h = NLevelSpace(:atom, 3, 2)
            i = Index(h, :i, 10, 1)
            j = Index(h, :j, 10, 1)
            σi = IndexedOperator(Transition(h, :σ, 2, 2), i)
            σj = IndexedOperator(Transition(h, :σ, 2, 2), j)
            expanded = expand_completeness(σi * σj)
            @test no_gs_projectors(expanded)
            @test isequal(expanded, expand_completeness(normal_order(σi * σj)))
        end

        @testset "σᵍᵍ * σⁱʲ (expand opt-in)" begin
            h = NLevelSpace(:atom, 3, 1)
            σgg = Transition(h, :σ, 1, 1)
            σ12 = Transition(h, :σ, 1, 2)
            @test isequal(
                expand_completeness(σgg * σ12),
                expand_completeness(normal_order(σgg * σ12)),
            )
        end

        @testset "User-constructed σᵍᵍ stays atomic (no composition fired)" begin
            h = NLevelSpace(:atom, 2, 1)
            σgg = Transition(h, :σ, 1, 1)
            @test is_transition(σgg)
            @test σgg.l1 == 1 && σgg.l2 == 1
            @test isequal(simplify(σgg), simplify(σgg))

            σ22 = Transition(h, :σ, 2, 2)
            @test isequal(
                expand_completeness(simplify(σgg)),
                expand_completeness(simplify(1 - σ22)),
            )
            @test isequal(
                expand_completeness(normal_order(σgg)),
                expand_completeness(simplify(1 - σ22)),
            )
        end

        @testset "expand_completeness removes σᵍᵍ from commutators" begin
            ha = NLevelSpace(:atom, 2, 1)
            hf = FockSpace(:cavity)
            h = hf ⊗ ha
            @qnumbers a::Destroy(h, 1)
            σ(α, β) = Transition(h, :σ, α, β, 2)

            H_jc = a' * σ(1, 2) + a * σ(2, 1)
            for op in (σ(1, 2), σ(2, 1), σ(2, 2), a' * σ(1, 2), a * σ(2, 1))
                result = expand_completeness(commutator(H_jc, op))
                result isa QAdd || continue
                @test no_gs_projectors(result)
            end
        end

        @testset "expand_completeness removes σᵍᵍ from indexed sums" begin
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

            for op in (
                    a' * σ(1, 2, j),
                    a * σ(2, 1, j),
                    σ(1, 2, j) * σ(2, 1, k),
                    σ(2, 2, j),
                )
                result = expand_completeness(commutator(H, op))
                result isa QAdd || continue
                @test no_gs_projectors(result)
            end
        end

        @testset "expand_completeness for higher-level systems" begin
            h = NLevelSpace(:atom, 4, 2)
            σ(α, β) = Transition(h, :σ, α, β)

            r1 = σ(1, 2) * σ(2, 1)
            @test no_gs_projectors(r1)

            r2 = expand_completeness(σ(2, 1) * σ(1, 2))
            @test no_gs_projectors(r2)
            @test length(r2) == 4

            r3 = σ(3, 2) * σ(2, 3)
            (term, _) = only(collect(r3))
            @test term.ops == [Transition(h, :σ, 3, 3)]
        end

        @testset "expand_completeness in ProductSpace with multiple NLevelSpaces" begin
            ha = NLevelSpace(:atomA, 2, 1)
            hb = NLevelSpace(:atomB, 3, 2)
            h = ha ⊗ hb

            σA(α, β) = Transition(h, :σA, α, β, 1)
            σB(α, β) = Transition(h, :σB, α, β, 2)

            rA = expand_completeness(σA(1, 2) * σA(2, 1))
            @test no_gs_projectors(rA)

            rB = expand_completeness(σB(2, 1) * σB(1, 2))
            @test no_gs_projectors(rB)

            rB2 = σB(1, 2) * σB(2, 1)
            term, _ = only(collect(rB2))
            @test length(term.ops) == 1
            @test term.ops[1].l1 == 1 && term.ops[1].l2 == 1
            @test term.ops[1].g == 2
        end
    end
end

@testset "Intern tables: out-of-range id is bounds-checked" begin
    # A stale/out-of-range Int32 name id (e.g. an Op deserialized across a session
    # boundary) must throw on read, not index out of bounds and corrupt memory.
    @test_throws BoundsError SecondQuantizedAlgebra.name_from_id(Int32(1_000_000))
    @test_throws BoundsError SecondQuantizedAlgebra.name_rank(Int32(1_000_000))
end

@testset "Intern tables: incremental name rank is lexicographic" begin
    # Names interned out of alphabetical order must still rank alphabetically
    # (the per-name rank splice must match a full re-sort).
    za = SecondQuantizedAlgebra.intern_name(:zzz_rank)
    an = SecondQuantizedAlgebra.intern_name(:aaa_rank)
    mi = SecondQuantizedAlgebra.intern_name(:mmm_rank)
    r(id) = SecondQuantizedAlgebra.name_rank(id)
    @test r(an) < r(mi) < r(za)
end
