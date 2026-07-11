using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils, SymReal, symtype
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QField, CNum, _to_cnum, _single_qadd,
    AvgFunc, sym_average, SumFunc, sym_sum, is_indexed_sum, sorted_arguments,
    has_sum_metadata, get_sum_indices, get_sum_non_equal,
    QTerm, QTermDict, NonEqualPair, order_key, qadd_order_key

@testset "Average" begin

    @testset "Construction — basic" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # average of a leaf operator
        avg_a = average(a)
        @test avg_a isa SymbolicUtils.BasicSymbolic{SymReal}
        @test SymbolicUtils.iscall(avg_a)
        @test symtype(avg_a) === Number
        @test is_average(avg_a)
        @test SymbolicUtils.operation(avg_a) isa AvgFunc
        # SymbolicUtils wraps QField args as Const; extract .val to check identity
        @test SymbolicUtils.arguments(avg_a)[1].val === a

        # average of adjoint is a distinct average
        avg_ad = average(ad)
        @test avg_ad isa SymbolicUtils.BasicSymbolic{SymReal}
        @test symtype(avg_ad) === Number
        @test SymbolicUtils.arguments(avg_ad)[1].val === ad
        @test !isequal(avg_a, avg_ad)

        # is_average type-stable dispatches
        @test !is_average(a)
        @test !is_average(3)
        @test !is_average(nothing)
    end

    @testset "average(QAdd) — prefactor extraction" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # Unit prefactor: average(a'a) = ⟨a'a⟩
        ada = ad * a
        avg_ada = average(ada)
        @test is_average(avg_ada)
        # SymbolicUtils wraps the QAdd argument as Const; extract via .val
        inner_wrapped = SymbolicUtils.arguments(avg_ada)[1]
        inner = SymbolicUtils.isconst(inner_wrapped) ? inner_wrapped.val : inner_wrapped
        @test inner isa QAdd
        @test length(operators(only(sorted_arguments(inner)))) == 2

        # Numeric prefactor: average(3 * a'a) = 3⟨a'a⟩
        three_ada = 3 * ad * a
        avg_3 = average(three_ada)
        @test !is_average(avg_3)  # prefactor pulled out

        # Pure scalar QAdd (empty operator key)
        scalar_add = _single_qadd(_to_cnum(5), Op[])
        avg_scalar = average(scalar_add)
        @test SymbolicUtils.isconst(avg_scalar) && avg_scalar.val == _to_cnum(5)

        # Scalar passthrough
        @test average(3) === 3
        @test average(0) === 0
    end

    @testset "average(QAdd) — linearity" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # average(a + a') distributes
        s = a + ad
        avg_s = average(s)
        @test avg_s isa SymbolicUtils.BasicSymbolic
        @test !(avg_s isa QField)
    end

    @testset "average(QAdd) — indexed sum metadata" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        a = Destroy(h, :a, 1)
        ad = a'
        @variables N
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        σ12 = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        # Σ_i(a' * σ_i^{12})
        s = Σ(ad * σ12, i)
        @test !isempty(s.indices)

        avg_s = average(s)
        @test avg_s isa SymbolicUtils.BasicSymbolic

        # The summed term becomes a dedicated `sym_sum` node carrying the scope.
        @test is_indexed_sum(avg_s)
        @test SymbolicUtils.operation(avg_s) isa SumFunc
        @test get_sum_indices(avg_s) == [i]
        @test isempty(get_sum_non_equal(avg_s))

        # Non-indexed QAdd produces no sum node
        plain = a + ad
        avg_plain = average(plain)
        @test !is_indexed_sum(avg_plain)
    end

    @testset "average(QAdd): index-dependent coefficient stays in sum (issue #175)" begin
        ha = NLevelSpace(:atoms, 2)
        hb = FockSpace(:bus)
        h = hb ⊗ ha
        @variables L::Real
        k = Index(h, :k, L, ha)
        kk = Index(h, :kk, L, ha)
        i = Index(h, :i, L, ha)
        j = Index(h, :j, L, ha)
        u(x, y) = DoubleIndexedVariable(:u, x, y)
        g(z) = IndexedVariable(:g, z)
        σ(z) = IndexedOperator(Transition(h, :σ, 1, 1), z)

        # Symptom 1: a pure c-number indexed sum keeps its Σ scope.
        c1 = average(Σ(u(kk, k), k))
        @test is_indexed_sum(c1)
        @test get_sum_indices(c1) == [k]
        @test undo_average(c1) == Σ(u(kk, k), k)

        # Symptom 2: an index-dependent coefficient lives INSIDE the Σ node body,
        # not hoisted out so `k` dangles. Check the off-diagonal node directly
        # (factor order within the body is Symbolics-internal and not asserted).
        c2 = average(Σ(u(kk, k) * σ(kk), k))
        @test undo_average(c2) == Σ(u(kk, k) * σ(kk), k)
        summands = SymbolicUtils.operation(c2) === (+) ? SymbolicUtils.arguments(c2) : [c2]
        offdiag = only(filter(is_indexed_sum, summands))
        @test get_sum_indices(offdiag) == [k]
        @test !isempty(get_sum_non_equal(offdiag))
        @test occursin("u(kk, k)", repr(SecondQuantizedAlgebra._sum_body(offdiag)))

        # Cancellation correctness: two sums over the same body but different
        # scope must NOT be `isequal` (metadata-only nodes would wrongly cancel).
        A = average(Σ(g(i) * σ(i), i))
        B = average(Σ(g(i) * σ(i), i, [j]))
        @test !isequal(A, B)
        @test !isequal(Symbolics.simplify(A - B), 0)

        # Grouping: terms sharing (indices, ne) collapse into one node body.
        h2(z) = IndexedVariable(:h, z)
        G = average(Σ(g(i) * σ(i) + h2(i) * IndexedOperator(Transition(h, :σ, 2, 2), i), i))
        @test is_indexed_sum(G)
        @test get_sum_indices(G) == [i]

        # Nested sums flatten into one scope; round-trips (compared modulo the
        # `.indices` ordering, which addition does not canonicalize).
        N = average(Σ(g(i) * g(j) * σ(i) * σ(j), i, j))
        @test iszero(simplify(undo_average(N) - Σ(g(i) * g(j) * σ(i) * σ(j), i, j)))
    end

    @testset "average(QAdd): singleton single-op NE is vacuous" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        @variables N
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        σ_i = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        bare_avg = average(σ_i)

        ne_qadd = QAdd(
            QTermDict(QTerm(Op[σ_i], NonEqualPair[(i, j)]) => _to_cnum(1)),
            Index[],
        )
        ne_avg = average(ne_qadd)

        @test isequal(bare_avg, ne_avg)
        @test !is_indexed_sum(ne_avg)

        σ_j = IndexedOperator(Transition(h, :σ, 2, 2, 2), j)
        prod_qadd_ne = QAdd(
            QTermDict(QTerm(Op[σ_i, σ_j], NonEqualPair[(i, j)]) => _to_cnum(1)),
            Index[],
        )
        prod_qadd_no_ne = QAdd(
            QTermDict(QTerm(Op[σ_i, σ_j], NonEqualPair[]) => _to_cnum(1)),
            Index[],
        )
        @test !isequal(average(prod_qadd_ne), average(prod_qadd_no_ne))
    end

    @testset "QAdd: dead NE pairs pruned at construction" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        @variables N
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        k = Index(h, :k, N, ha)
        σ_i = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        dead = QAdd(
            QTermDict(QTerm(Op[σ_i], NonEqualPair[(j, k)]) => _to_cnum(1)),
            Index[],
        )
        @test all(isempty(t.ne) for t in keys(dead.arguments))

        live_op = QAdd(
            QTermDict(QTerm(Op[σ_i], NonEqualPair[(i, j)]) => _to_cnum(1)),
            Index[],
        )
        @test any(!isempty(t.ne) for t in keys(live_op.arguments))

        live_scope = QAdd(
            QTermDict(QTerm(Op[σ_i], NonEqualPair[(j, k)]) => _to_cnum(1)),
            Index[j],
        )
        @test any(!isempty(t.ne) for t in keys(live_scope.arguments))

        @test isequal(QAdd(deepcopy(dead.arguments), Index[]), dead)

        gjk = DoubleIndexedVariable(:g, j, k)
        live_coef = QAdd(
            QTermDict(
                QTerm(Op[σ_i], NonEqualPair[(j, k)]) => _to_cnum(gjk),
            ),
            Index[],
        )
        @test any(!isempty(t.ne) for t in keys(live_coef.arguments))

        clash = QAdd(
            QTermDict(
                QTerm(Op[σ_i], NonEqualPair[(j, k)]) => _to_cnum(1),
                QTerm(Op[σ_i], NonEqualPair[]) => _to_cnum(2),
            ),
            Index[],
        )
        @test length(clash.arguments) == 1
        (only_term, only_c) = first(clash.arguments)
        @test isempty(only_term.ne)
        @test isequal(only_c, _to_cnum(3))
    end

    @testset "undo_average: always returns QAdd" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # QSym round-trip → QAdd wrapping the operator
        r = undo_average(average(a))
        @test r isa QAdd
        @test isequal(r, _single_qadd(_to_cnum(1), Op[a]))

        r2 = undo_average(average(ad))
        @test r2 isa QAdd
        @test isequal(r2, _single_qadd(_to_cnum(1), Op[ad]))

        # Number → QAdd with scalar prefactor and empty ops
        @test undo_average(3) isa QAdd
        @test isequal(undo_average(3), _single_qadd(_to_cnum(3), Op[]))

        @test undo_average(0.5) isa QAdd
        @test isequal(undo_average(0.5), _single_qadd(_to_cnum(0.5), Op[]))

        # QField → QAdd
        @test undo_average(a) isa QAdd
        @test isequal(undo_average(a), _single_qadd(_to_cnum(1), Op[a]))

        # QAdd passthrough
        qadd = ad * a
        @test undo_average(qadd) === qadd

        # Num round-trip
        avg_num = Symbolics.Num(average(a))
        result_num = undo_average(avg_num)
        @test result_num isa QAdd

        # Equation round-trip: undo_average produces QAdd (not BasicSymbolic), returns Pair
        eq = Symbolics.Equation(average(a), average(ad))
        result_eq = undo_average(eq)
        @test result_eq isa Pair
        @test result_eq.first isa QAdd
        @test result_eq.second isa QAdd

        # Equation with sum-of-averages
        avg_sum_lhs = average(a) + average(ad)
        avg_sum_rhs = average(a) + average(ad)
        eq2 = Symbolics.Equation(avg_sum_lhs, avg_sum_rhs)
        result_eq2 = undo_average(eq2)
        @test result_eq2 isa Pair
        @test result_eq2.first isa QAdd
        @test result_eq2.second isa QAdd
    end

    @testset "undo_average — sum of averages" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # undo_average on sum of averages
        avg_sum = average(a) + average(ad)
        result = undo_average(avg_sum)
        @test result isa QAdd
    end

    @testset "undo_average — indexed round-trip" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        a = Destroy(h, :a, 1)
        ad = a'
        @variables N
        i = Index(h, :i, N, ha)
        σ12 = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        s = Σ(ad * σ12, i)
        avg_s = average(s)

        result = undo_average(avg_s)
        @test result isa QAdd
        @test result.indices == [i]
    end

    @testset "undo_average preserves distinct scoped terms" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        @variables N_scope
        i = Index(h, :i, N_scope, ha)
        j = Index(h, :j, N_scope, ha)
        bi = IndexedOperator(Destroy(h, :b, 1), i)

        expr = Σ(bi, i) + Σ(bi, i, [j])
        roundtrip = undo_average(average(expr))

        @test isequal(roundtrip, expr)
        @test length(roundtrip) == 2
    end

    @testset "Metadata helpers" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        a = Destroy(h, :a, 1)
        ad = a'
        @variables N
        i = Index(h, :i, N, ha)
        σ12 = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        s = Σ(ad * σ12, i)
        avg_s = average(s)

        @test has_sum_metadata(avg_s)
        @test get_sum_indices(avg_s) == [i]
        @test isempty(get_sum_non_equal(avg_s))

        # Non-indexed expression
        avg_a = average(a)
        @test !has_sum_metadata(avg_a)

        # Fallback for non-BasicSymbolic types
        @test !has_sum_metadata(3)
        @test !has_sum_metadata(a)
    end

    @testset "acts_on" begin
        hf = FockSpace(:f)
        hn = NLevelSpace(:n, 2, 1)
        h = hf ⊗ hn
        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)

        # QSym
        @test acts_on(a) == [1]
        @test acts_on(σ) == [2]

        # QAdd (product)
        @test acts_on(a' * σ) == [1, 2]

        # QAdd
        @test acts_on(a + a') == [1]

        # Average
        avg = average(a' * σ)
        @test acts_on(avg) == [1, 2]

        # Sum of averages
        avg_sum = average(a) + average(σ)
        @test acts_on(avg_sum) == [1, 2]
    end

    @testset "Printing" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        @test repr(average(a)) == "⟨a⟩"
        @test repr(average(ad)) == "⟨a'⟩"
        @test repr(average(ad * a)) == "⟨a' * a⟩"
        @test repr(average(3 * ad * a)) == "3⟨a' * a⟩"
        avg_sum_repr = repr(average(a + ad))
        @test avg_sum_repr == "⟨a⟩ + ⟨a'⟩" || avg_sum_repr == "⟨a'⟩ + ⟨a⟩"

        # Multi-space
        hf = FockSpace(:f)
        hn = NLevelSpace(:n, 2, 1)
        hp = hf ⊗ hn
        b = Destroy(hp, :b, 1)
        σ = Transition(hp, :σ, 1, 2, 2)
        @test repr(average(b' * σ)) == "⟨b' * σ₁₂⟩"
    end

    @testset "LaTeX" begin
        using Latexify: latexify
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        @test string(latexify(average(a))) ==
            "\\begin{equation}\n\\langle a \\rangle\n\\end{equation}\n"
        @test string(latexify(average(ad * a))) ==
            "\\begin{equation}\n\\langle a^{\\dagger}a \\rangle\n\\end{equation}\n"
        @test string(latexify(average(3 * ad * a))) ==
            "\\begin{equation}\n3 ~ \\langle a^{\\dagger}a \\rangle\n\\end{equation}\n"
    end

    @testset "Arithmetic composition" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        avg_a = average(a)
        avg_ad = average(ad)

        # Addition of averages
        s = avg_a + avg_ad
        @test s isa SymbolicUtils.BasicSymbolic

        # Multiplication by scalar
        p = 2 * avg_a
        @test p isa SymbolicUtils.BasicSymbolic

        # Average of average (passthrough — averages are BasicSymbolic scalars)
        @test isequal(average(avg_a), avg_a)
    end

    @testset "Average extracts indexed variable prefactor" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        @variables N
        i = Index(h, :i, N, ha)
        gi = IndexedVariable(:g, i)
        σi = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)

        # average of indexed variable alone is itself (it's a c-number)
        @test isequal(average(gi), gi)

        # average(2*σ) = 2*average(σ)
        @test isequal(average(2 * σi), 2 * average(σi))

        # average(g(i)*σ(i)) pulls g(i) out as prefactor
        avg = average(gi * Transition(h, :σ, 2, 2, 2))
        @test SymbolicUtils.iscall(avg)
    end

    @testset "Commutativity of average products" begin
        hf = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2, 1)
        h = hf ⊗ ha
        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)

        @test isequal(
            average(a * σ) * average(a) - average(a) * average(a * σ),
            Symbolics.Num(0),
        )
    end

    @testset "get_indices on averages" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:f)
        h = hf ⊗ ha
        @variables N_gi
        i = Index(h, :i, N_gi, ha)
        j = Index(h, :j, N_gi, ha)
        σi = IndexedOperator(Transition(h, :σ, 1, 2, 2), i)
        σj = IndexedOperator(Transition(h, :σ, 2, 1, 2), j)

        # get_indices on QAdd of indexed operators (dict iteration order undefined)
        @test Set(get_indices(σi + σj)) == Set([i, j])

        inds = get_indices(average(σi) + 3 + average(σj))
        @test i in inds
        @test j in inds
        @test length(inds) == 2
    end

    @testset "C-number handling" begin
        hf = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2, 1)
        h = hf ⊗ ha

        a = Destroy(h, :a, 1)

        @variables ωc::Real
        @test isequal(average(ωc), ωc)
        # average pulls scalar prefactors out
        avg_ωₐ = average(ωc * a)
        @test is_average(avg_ωₐ) || SymbolicUtils.iscall(avg_ωₐ)
    end

    @testset "Double average (QC#242)" begin
        @variables Δ_::Real g::Real κ::Real η::Real

        hf = FockSpace(:cavity)
        ha1 = NLevelSpace(:atom1, 2, 1)
        ha2 = NLevelSpace(:atom2, 2, 1)
        h = hf ⊗ ha1 ⊗ ha2

        a = Destroy(h, :a, 1)
        s1(i, j) = Transition(h, :s1, i, j, 2)
        s2(i, j) = Transition(h, :s2, i, j, 3)

        H =
            Δ_ * a' * a +
            g * (a' + a) * (s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2)) +
            η * (a' + a)

        imH = im * H
        op_ = a' * a
        rhs_ = commutator(imH, op_)
        rhs_avg = average(rhs_)
        rhs_avg_simplified = SymbolicUtils.simplify(rhs_avg)

        terms = undo_average(rhs_avg_simplified)
        @test terms isa QAdd
        @test terms isa QAdd
    end

    @testset "Average addition/multiplication (PR 28)" begin
        ha1 = NLevelSpace(:atom1, 2, 1)
        ha2 = NLevelSpace(:atom2, 2, 1)
        h = ha1 ⊗ ha2
        s1(i, j) = Transition(h, :s1, i, j, 1)
        s2(i, j) = Transition(h, :s2, i, j, 2)

        expr = average(s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2))
        @test expr isa SymbolicUtils.BasicSymbolic

        expr2 = simplify(average(s1(2, 1) + s1(1, 2)))
        @test expr2 isa SymbolicUtils.BasicSymbolic
    end

    @testset "Multi-space averaging" begin
        hf = FockSpace(:f)
        hn = NLevelSpace(:n, 2, 1)
        h = hf ⊗ hn
        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)

        avg1 = average(a)
        avg2 = average(σ)
        @test !isequal(avg1, avg2)

        # Product across spaces
        avg_prod = average(a' * σ)
        @test is_average(avg_prod)
    end

    @testset "average(QAdd) — no literal `complex(re, im)` in result" begin
        # Background: `SymbolicUtils.unwrap(::Complex{<:Num})` (Symbolics) builds
        # a `Term(complex, [re, img])` symbolic call whenever both re/img are
        # symbolic. That literal call is opaque to `simplify` / `expand`. It
        # also generates a runtime `complex(::Real, ::Complex)` call when MTK
        # codegens the equation, for which Base has no method.
        #
        # `average(::QAdd)` must avoid materialising such literals. It does so
        # by accumulating into a `Num` (or `BasicSymbolic{SymReal}`) result,
        # bringing in `im` via `Symbolics.IM` (the BasicSymbolic sym for `im`)
        # rather than via a `Complex{Num}` intermediate.

        # Walk an expression tree and collect all `Term`s whose operation is `complex`.
        function _find_complex_terms(x, out = Any[])
            if x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x)
                SymbolicUtils.operation(x) === complex && push!(out, x)
                for arg in SymbolicUtils.arguments(x)
                    _find_complex_terms(arg, out)
                end
            end
            return out
        end

        @testset "operator branch: `i * im * avg` does not leak literals" begin
            h = FockSpace(:c)
            a = Destroy(h, :a)
            @variables Δ::Real
            # `commutator(im * Δ * a' * a, a)` ends up as `-Δim * a` (i.e. a QAdd
            # whose single coefficient is `Complex{Num}(0, -Δ)`). Averaging it
            # previously produced `⟨a⟩ * complex(0, -Δ)`.
            qadd = (-im * Δ) * a
            avg = average(qadd)
            @test isempty(_find_complex_terms(avg))
            # The result still represents the same mathematical quantity. We
            # build the comparison target via `Symbolics.IM * average(a)`
            # (additive form) rather than via `(-im) * Δ * average(a)`, because
            # the latter promotes through `Complex{Num}` and itself materialises
            # a literal `complex(...)` Term that defeats the diff.
            target = -(Δ * Symbolics.IM * average(a))
            @test SymbolicUtils._iszero(
                SymbolicUtils.simplify(avg - target; expand = true)
            )
        end

        @testset "constant branch: `result += c` does not leak literals" begin
            h = FockSpace(:c)
            a = Destroy(h, :a)
            @variables η::Real
            # `commutator(im*η*(a + a'), a')` = `i*η`, a pure-scalar QAdd (no ops)
            # — exercises the `isempty(term.ops)` branch of `average(::QAdd)`.
            scalar = _single_qadd(_to_cnum(Complex(0, η)), Op[])
            avg = average(scalar)
            @test isempty(_find_complex_terms(avg))
            # Mixed re + im constant must also stay literal-free.
            @variables r::Real
            mixed = _single_qadd(_to_cnum(Complex(r, η)), Op[])
            avg_mixed = average(mixed)
            @test isempty(_find_complex_terms(avg_mixed))
        end

        @testset "operator + constant in the same QAdd" begin
            h = FockSpace(:c)
            a = Destroy(h, :a)
            @variables Δ::Real η::Real
            # Build a QAdd that has both a scalar (im*η) and an operator (im*Δ*a)
            # coefficient with imaginary parts. Previously, accumulating these in
            # one `result::Num` promoted to `Complex{Num}` on the first scalar
            # add and corrupted every subsequent addition.
            qa = (im * η) * one(a) + (-im * Δ) * a
            avg = average(qa)
            @test isempty(_find_complex_terms(avg))
        end

        @testset "real-only coefficients are unaffected" begin
            h = FockSpace(:c)
            a = Destroy(h, :a)
            @variables κ::Real
            avg = average((κ / 2) * a)
            @test isempty(_find_complex_terms(avg))
            @test isequal(avg, (κ / 2) * average(a))
        end
    end

    @testset "Type stability" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        i = Index(h, :i, 3, h)
        x = average(Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i))

        # average / undo_average: round-trip pins to BasicSymbolic / QAdd
        @test @inferred(average(a)) isa SymbolicUtils.BasicSymbolic
        @test @inferred(average(a + a')) isa SymbolicUtils.BasicSymbolic
        @test @inferred(average(2 * a)) isa SymbolicUtils.BasicSymbolic
        @test @inferred(undo_average(average(a))) isa QAdd
        @test @inferred(undo_average(SymbolicUtils.unwrap(average(a' * a)))) isa QAdd

        # metadata accessors infer to their concrete return types
        @test @inferred(get_sum_indices(x)) isa Vector{Index}
        @test @inferred(get_sum_non_equal(x)) isa Vector{Tuple{Index, Index}}

        # recursive walks reseal through `_idx_seal` / `_aon_seal`
        @test @inferred(SecondQuantizedAlgebra.get_indices(SymbolicUtils.unwrap(average(a' * a)))) isa Vector{Index}
        @test @inferred(SecondQuantizedAlgebra.acts_on(SymbolicUtils.unwrap(average(a' * a)))) isa Vector{Int}
        @test @inferred(SecondQuantizedAlgebra.acts_on(a' * a)) isa Vector{Int}

        # Predicates stay Bool from any input
        @test @inferred(is_average(average(a))) isa Bool
        @test @inferred(is_average(a)) isa Bool
        @test @inferred(has_sum_metadata(x)) isa Bool

        # Seals: canaries against a regression that re-introduces an `Any` path.
        @test Base.return_types(SecondQuantizedAlgebra._idx_seal, (Any,))[1] === Vector{Index}
        @test Base.return_types(SecondQuantizedAlgebra._aon_seal, (Any,))[1] === Vector{Int}
    end

    @testset "Display of summed averages (issue #204)" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @variables x::Number

        @test sprint(show, MIME("text/plain"), average(a) + average(a')) == "⟨a⟩ + ⟨a'⟩"
        @test sprint(show, MIME("text/plain"), x * average(a) + x * average(a')) == "⟨a⟩*x + ⟨a'⟩*x"
        @test sprint(show, MIME("text/plain"), conj(x) * average(a) + conj(x) * average(a')) ==
            "⟨a⟩*conj(x) + ⟨a'⟩*conj(x)"
        @test sprint(show, MIME("text/plain"), average(a' * a) + average(a * a')) == "1 + 2⟨a' * a⟩"

        @test isless(a, a') == isless(order_key(a), order_key(a'))
        @test !isless(a, a)
        @test isless(a, a') ⊻ isless(a', a)
        @test sort([a', a]) == [a, a']
        @test isless(a' * a, a * a') == isless(qadd_order_key(a' * a), qadd_order_key(a * a'))
    end

end
