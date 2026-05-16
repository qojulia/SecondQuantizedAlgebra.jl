using Test
using SecondQuantizedAlgebra
using Symbolics: @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, Transition, NO_INDEX,
    _partial_sort!, _site_compare, _can_commute, _commute_pair,
    _reduce_pair, _ground_state_expand,
    SiteCmp, Less, Greater, Equal, Undetermined,
    _CNUM_ONE, _CNUM_ZERO, _EMPTY_NE,
    _canonicalize_to_dict!, QTermDict,
    _reduce_ops, _commute_ops, _expand_gs_ops, _substitute_ops,
    _stream!, _canonicalize!, Index

# Tests of private hooks, passes, sort, and canonical-form invariants.

# Helper: is this op a ground-state projector of an NLevelSpace?
function _is_gs_projector(op)
    op isa Transition || return false
    return op.i == op.ground_state && op.j == op.ground_state
end

# Helper: every dict key in `expr` is free of ground-state projectors
function _no_gs_projectors(expr::QAdd)
    for term in keys(expr.arguments)
        any(_is_gs_projector, term.ops) && return false
    end
    return true
end

# Helper: a term is canonical when no neighboring pair could swap or reduce.
function _is_canonical(t)
    ops = t.ops
    ne = t.ne
    for i in 1:(length(ops) - 1)
        cmp = _site_compare(ops[i], ops[i + 1], ne)
        if cmp === Greater
            return false
        end
        if cmp === Equal
            _can_commute(ops[i], ops[i + 1]) || return false
            _reduce_pair(ops[i], ops[i + 1]) === nothing || return false
        end
    end
    return true
end

@testset "Internals" begin

    @testset "Operator hooks" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        ad = adjoint(a)
        @test _site_compare(a, ad, _EMPTY_NE) === Equal
        @test _site_compare(ad, a, _EMPTY_NE) === Equal
        @test _can_commute(a, ad) === false
        @test _can_commute(ad, a) === true
        sw = _commute_pair(a, ad)
        @test sw[1] === ad && sw[2] === a

        ha = NLevelSpace(:atom, 2)
        σ12 = Transition(ha, :σ, 1, 2)
        σ21 = Transition(ha, :σ, 2, 1)
        @test _site_compare(σ12, σ21, _EMPTY_NE) === Equal
        @test _can_commute(σ12, σ21) === false
        red = _reduce_pair(σ21, σ12)   # σ²¹·σ¹² = σ²²
        @test red isa Tuple
        @test red[1].i == 2 && red[1].j == 2
        @test _ground_state_expand(Transition(ha, :σ, 1, 1)) !== nothing
    end

    @testset "Partial sort" begin
        @testset "distinct sites" begin
            h = FockSpace(:c) ⊗ NLevelSpace(:a, 2)
            a = Destroy(h, :a, 1)
            σ = Transition(h, :σ, 1, 2)
            ops = QSym[σ, a]
            _partial_sort!(ops, _EMPTY_NE)
            @test ops[1] isa Destroy
            @test ops[2] isa Transition
        end

        @testset "same-site preserved" begin
            h = NLevelSpace(:a, 3)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            ops = QSym[σ12, σ21]
            pre = copy(ops)
            _partial_sort!(ops, _EMPTY_NE)
            @test ops == pre
        end

        @testset "undetermined preserved" begin
            ha = NLevelSpace(:a, 2)
            i = Index(ha, :i, 5, ha)
            j = Index(ha, :j, 5, ha)
            σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
            σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
            ops = QSym[σ_j, σ_i]
            _partial_sort!(ops, _EMPTY_NE)
            @test ops[1] === σ_j
            @test ops[2] === σ_i
        end

        @testset "ne resolves undetermined to distinct" begin
            ha = NLevelSpace(:a, 2)
            i = Index(ha, :i, 5, ha)
            j = Index(ha, :j, 5, ha)
            σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
            σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
            ops = QSym[σ_j, σ_i]
            _partial_sort!(ops, [(i, j)])
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
            _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
            @test length(out) == 1
        end

        @testset "canonicalize_to_dict!: like-term collection" begin
            out = QTermDict()
            _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
            _canonicalize_to_dict!(out, QSym[a], _CNUM_ONE, _EMPTY_NE)
            @test length(out) == 1
            @test first(values(out)) == 2 + 0im
        end

        @testset "canonicalize_to_dict!: zero coeff dropped" begin
            out = QTermDict()
            _canonicalize_to_dict!(out, QSym[a], _CNUM_ZERO, _EMPTY_NE)
            @test isempty(out)
        end

        @testset "_reduce_ops: Transition composition" begin
            σ12 = Transition(hn, :σ, 1, 2)
            σ23 = Transition(hn, :σ, 2, 3)
            emitted = Tuple{Vector{QSym}, CNum}[]
            _reduce_ops(QSym[σ12, σ23], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test length(emitted[1][1]) == 1
            @test emitted[1][1][1].i == 1 && emitted[1][1][1].j == 3
        end

        @testset "_reduce_ops: zero from incompatible composition" begin
            σ12 = Transition(hn, :σ, 1, 2)
            σ31 = Transition(hn, :σ, 3, 1)
            emitted = Tuple{Vector{QSym}, CNum}[]
            _reduce_ops(QSym[σ12, σ31], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test isempty(emitted)
        end

        @testset "_reduce_ops: no-op input passes through" begin
            emitted = Tuple{Vector{QSym}, CNum}[]
            _reduce_ops(QSym[ad, a], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == QSym[ad, a]
        end

        @testset "_commute_ops: Fock aa† → a†a + 1" begin
            emitted = Tuple{Vector{QSym}, CNum}[]
            _commute_ops(QSym[a, ad], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 2
            sort!(emitted, by = e -> length(e[1]))
            @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
            @test emitted[2][1] == QSym[ad, a]
        end

        @testset "_commute_ops: no-op on already-ordered pair" begin
            emitted = Tuple{Vector{QSym}, CNum}[]
            _commute_ops(QSym[ad, a], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == QSym[ad, a]
        end

        @testset "_expand_gs_ops: σ¹¹ → 1 - σ²²" begin
            h2 = NLevelSpace(:a, 2)
            σ11 = Transition(h2, :σ, 1, 1)
            emitted = Tuple{Vector{QSym}, CNum}[]
            _expand_gs_ops(QSym[σ11], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 2
            sort!(emitted, by = e -> length(e[1]))
            @test isempty(emitted[1][1]) && emitted[1][2] == _CNUM_ONE
            @test length(emitted[2][1]) == 1 && emitted[2][2] == -_CNUM_ONE
            @test emitted[2][1][1].i == 2 && emitted[2][1][1].j == 2
        end

        @testset "_expand_gs_ops: passthrough when no σᵍᵍ" begin
            emitted = Tuple{Vector{QSym}, CNum}[]
            _expand_gs_ops(QSym[a], _CNUM_ONE) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][1] == QSym[a]
        end

        @testset "_substitute_ops: operator → scalar" begin
            emitted = Tuple{Vector{QSym}, CNum}[]
            _substitute_ops(QSym[a, ad], _CNUM_ONE, Dict(a => 2)) do o, c
                push!(emitted, (copy(o), c))
            end
            @test length(emitted) == 1
            @test emitted[1][2] == 2 + 0im
            @test emitted[1][1] == QSym[ad]
        end

        @testset "_stream!: idempotent on canonical input" begin
            out = QTermDict()
            _stream!(out, QSym[ad, a], _CNUM_ONE, _EMPTY_NE)
            @test length(out) == 1
        end

        @testset "_canonicalize!: aa† → a†a + 1" begin
            out = QTermDict()
            _canonicalize!(out, QSym[a, ad], _CNUM_ONE, _EMPTY_NE)
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
                @test _is_canonical(t)
            end
        end

        @testset "after normal_order" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = a * adjoint(a) * a
            for (t, _) in normal_order(q)
                @test _is_canonical(t)
            end
        end

        @testset "after commutator" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = commutator(a, adjoint(a))
            for (t, _) in q
                @test _is_canonical(t)
            end
        end

        @testset "after substitute" begin
            h = FockSpace(:f)
            a = Destroy(h, :a)
            q = adjoint(a) * a
            sub = SecondQuantizedAlgebra.SymbolicUtils.substitute(q, Dict{Any, Any}(a => 0.5))
            for (t, _) in sub
                @test _is_canonical(t)
            end
        end

        @testset "after expand_completeness" begin
            h = NLevelSpace(:a, 3)
            σ11 = Transition(h, :σ, 1, 1)
            exp = expand_completeness(σ11)
            for (t, _) in exp
                @test _is_canonical(t)
            end
        end
    end

    @testset "Canonical form (NLevelSpace ground-state projection)" begin
        @testset "Transition carries ground_state and n_levels" begin
            h = NLevelSpace(:atom, 3, 2)
            σ = Transition(h, :σ, 1, 3)
            @test σ.ground_state == 2
            @test σ.n_levels == 3

            h2 = NLevelSpace(:atom, 2)
            σ2 = Transition(h2, :σ, 1, 2)
            @test σ2.ground_state == 1
            @test σ2.n_levels == 2

            hf = FockSpace(:c)
            hp = hf ⊗ h
            σp = Transition(hp, :σ, 1, 3, 2)
            @test σp.ground_state == 2
            @test σp.n_levels == 3
        end

        @testset "Field preservation through ops" begin
            h = NLevelSpace(:atom, 4, 3)
            σ = Transition(h, :σ, 1, 2)

            σadj = adjoint(σ)
            @test σadj.ground_state == 3
            @test σadj.n_levels == 4

            i = Index(h, :i, 10, h)
            σi = IndexedOperator(σ, i)
            @test σi.ground_state == 3
            @test σi.n_levels == 4

            j = Index(h, :j, 10, h)
            σj = change_index(σi, i, j)
            @test σj.ground_state == 3
            @test σj.n_levels == 4
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
            @test _no_gs_projectors(result)
        end

        @testset "Composition+expand: σ¹²·σ²¹ → 1 - σ²² - σ³³ (3-level, g=1)" begin
            h = NLevelSpace(:atom, 3, 1)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ22 = Transition(h, :σ, 2, 2)
            σ33 = Transition(h, :σ, 3, 3)
            result = expand_completeness(σ12 * σ21)
            @test isequal(result, simplify(1 - σ22 - σ33))
            @test _no_gs_projectors(result)
        end

        @testset "Composition+expand: ground state ≠ 1 (g=2)" begin
            h = NLevelSpace(:atom, 3, 2)
            σ12 = Transition(h, :σ, 1, 2)
            σ21 = Transition(h, :σ, 2, 1)
            σ11 = Transition(h, :σ, 1, 1)
            σ33 = Transition(h, :σ, 3, 3)

            result = expand_completeness(σ21 * σ12)
            @test isequal(result, simplify(1 - σ11 - σ33))
            @test _no_gs_projectors(result)

            result2 = σ12 * σ21
            @test length(result2) == 1
            (term, _c) = only(collect(result2))
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
            @test _no_gs_projectors(result)
            for term in keys(result.arguments)
                for op in term.ops
                    if op isa Transition
                        @test op.index == i
                        @test op.ground_state == 1
                        @test op.n_levels == 2
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
            @test _no_gs_projectors(result)
        end

        @testset "Different-site σᵍᵍ_i · σᵍᵍ_j (expand via expand_completeness)" begin
            h = NLevelSpace(:atom, 3, 2)
            i = Index(h, :i, 10, 1)
            j = Index(h, :j, 10, 1)
            σi = IndexedOperator(Transition(h, :σ, 2, 2), i)
            σj = IndexedOperator(Transition(h, :σ, 2, 2), j)
            expanded = expand_completeness(σi * σj)
            @test _no_gs_projectors(expanded)
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
            @test σgg isa Transition
            @test σgg.i == 1 && σgg.j == 1
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
                @test _no_gs_projectors(result)
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
                @test _no_gs_projectors(result)
            end
        end

        @testset "expand_completeness for higher-level systems" begin
            h = NLevelSpace(:atom, 4, 2)
            σ(α, β) = Transition(h, :σ, α, β)

            r1 = σ(1, 2) * σ(2, 1)
            @test _no_gs_projectors(r1)

            r2 = expand_completeness(σ(2, 1) * σ(1, 2))
            @test _no_gs_projectors(r2)
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
            @test _no_gs_projectors(rA)

            rB = expand_completeness(σB(2, 1) * σB(1, 2))
            @test _no_gs_projectors(rB)

            rB2 = σB(1, 2) * σB(2, 1)
            term, _ = only(collect(rB2))
            @test length(term.ops) == 1
            @test term.ops[1].i == 1 && term.ops[1].j == 1
            @test term.ops[1].ground_state == 2
        end
    end
end
