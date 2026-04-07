using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils, SymReal, symtype
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QField, CNum, _to_cnum, _single_qadd,
    QTermDict, AvgSym, AvgFunc, sym_average, SumIndices, SumNonEqual, sorted_arguments

@testset "Average" begin

    @testset "Construction — basic" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # average of a leaf operator
        avg_a = average(a)
        @test avg_a isa SymbolicUtils.BasicSymbolic{SymReal}
        @test SymbolicUtils.iscall(avg_a)
        @test symtype(avg_a) === AvgSym
        @test is_average(avg_a)
        @test SymbolicUtils.operation(avg_a) isa AvgFunc
        # SymbolicUtils wraps QField args as Const; extract .val to check identity
        @test SymbolicUtils.arguments(avg_a)[1].val === a

        # average of adjoint is a distinct average
        avg_ad = average(ad)
        @test avg_ad isa SymbolicUtils.BasicSymbolic{SymReal}
        @test symtype(avg_ad) === AvgSym
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

        # Unit prefactor: average(a†a) = ⟨a†a⟩
        ada = ad * a
        avg_ada = average(ada)
        @test is_average(avg_ada)
        # SymbolicUtils wraps the QAdd argument as Const; extract via .val
        inner_wrapped = SymbolicUtils.arguments(avg_ada)[1]
        inner = SymbolicUtils.isconst(inner_wrapped) ? inner_wrapped.val : inner_wrapped
        @test inner isa QAdd
        @test length(operators(only(sorted_arguments(inner)))) == 2

        # Numeric prefactor: average(3 * a†a) = 3⟨a†a⟩
        three_ada = 3 * ad * a
        avg_3 = average(three_ada)
        @test !is_average(avg_3)  # prefactor pulled out

        # Pure scalar QAdd (empty operator key)
        scalar_add = QAdd(QTermDict(QSym[] => _to_cnum(5)), Index[], Tuple{Index, Index}[])
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

        # average(a + a†) distributes
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

        # Σ_i(a† * σ_i^{12})
        s = Σ(ad * σ12, i)
        @test !isempty(s.indices)

        avg_s = average(s)
        @test avg_s isa SymbolicUtils.BasicSymbolic

        # Metadata should be preserved
        @test SymbolicUtils.hasmetadata(avg_s, SumIndices)
        @test SymbolicUtils.getmetadata(avg_s, SumIndices) == [i]
        @test SymbolicUtils.hasmetadata(avg_s, SumNonEqual)

        # Non-indexed QAdd has no metadata
        plain = a + ad
        avg_plain = average(plain)
        if avg_plain isa SymbolicUtils.BasicSymbolic
            @test !SymbolicUtils.hasmetadata(avg_plain, SumIndices)
        end
    end

    @testset "undo_average — always returns QAdd" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # QSym round-trip → QAdd wrapping the operator
        r = undo_average(average(a))
        @test r isa QAdd
        @test isequal(r, _single_qadd(_to_cnum(1), QSym[a]))

        r2 = undo_average(average(ad))
        @test r2 isa QAdd
        @test isequal(r2, _single_qadd(_to_cnum(1), QSym[ad]))

        # Number → QAdd with scalar prefactor and empty ops
        @test undo_average(3) isa QAdd
        @test isequal(undo_average(3), _single_qadd(_to_cnum(3), QSym[]))

        @test undo_average(0.5) isa QAdd
        @test isequal(undo_average(0.5), _single_qadd(_to_cnum(0.5), QSym[]))

        # QField → QAdd
        @test undo_average(a) isa QAdd
        @test isequal(undo_average(a), _single_qadd(_to_cnum(1), QSym[a]))

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
        @test get_sum_non_equal(avg_s) == Tuple{Index, Index}[]

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
        @test repr(average(ad)) == "⟨a†⟩"
        @test repr(average(ad * a)) == "⟨a† * a⟩"
        @test repr(average(3 * ad * a)) == "3⟨a† * a⟩"
        avg_sum_repr = repr(average(a + ad))
        @test avg_sum_repr == "⟨a⟩ + ⟨a†⟩" || avg_sum_repr == "⟨a†⟩ + ⟨a⟩"

        # Multi-space
        hf = FockSpace(:f)
        hn = NLevelSpace(:n, 2, 1)
        hp = hf ⊗ hn
        b = Destroy(hp, :b, 1)
        σ = Transition(hp, :σ, 1, 2, 2)
        @test repr(average(b' * σ)) == "⟨b† * σ₁₂⟩"
    end

    @testset "LaTeX" begin
        using Latexify: latexify
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # LaTeX renders as \mathrm{avg}(...) — Symbolics' _toexpr pipeline
        # has no non-piracy hook for custom operations. Upstream fix needed.
        @test string(latexify(average(a))) ==
            "\\begin{equation}\n\\mathrm{avg}\\left( a \\right)\n\\end{equation}\n"
        @test string(latexify(average(ad * a))) ==
            "\\begin{equation}\n\\mathrm{avg}\\left( a^{\\dagger}a \\right)\n\\end{equation}\n"
        @test string(latexify(average(3 * ad * a))) ==
            "\\begin{equation}\n3 ~ \\mathrm{avg}\\left( a^{\\dagger}a \\right)\n\\end{equation}\n"
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

        # get_indices on QAdd of indexed operators
        @test isequal(get_indices(σi + σj), [i, j])

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
        avg_ωa = average(ωc * a)
        @test is_average(avg_ωa) || SymbolicUtils.iscall(avg_ωa)
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

end
