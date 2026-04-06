using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils, SymReal, symtype
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: QMul, QAdd, QSym, QField, CNum, _to_cnum, _to_qmul,
    AvgSym, AvgFunc, sym_average, SumIndices, SumNonEqual

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

    @testset "average(QMul) — prefactor extraction" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # Unit prefactor: average(a†a) = ⟨a†a⟩
        ada = ad * a
        avg_ada = average(ada)
        @test is_average(avg_ada)
        # SymbolicUtils wraps the QMul argument as Const; extract via .val
        inner_wrapped = SymbolicUtils.arguments(avg_ada)[1]
        inner = SymbolicUtils.isconst(inner_wrapped) ? inner_wrapped.val : inner_wrapped
        @test inner isa QMul
        @test length(inner.args_nc) == 2

        # Numeric prefactor: average(3 * a†a) = 3⟨a†a⟩
        three_ada = 3 * ad * a
        avg_3 = average(three_ada)
        @test !is_average(avg_3)  # prefactor pulled out

        # Pure scalar QMul (empty args_nc)
        scalar_mul = QMul(_to_cnum(5), QSym[])
        @test average(scalar_mul) == _to_cnum(5)

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

    @testset "undo_average — round-trip" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # QSym round-trip
        @test undo_average(average(a)) === a
        @test undo_average(average(ad)) === ad

        # Number passthrough
        @test undo_average(3) === 3
        @test undo_average(0.5) === 0.5

        # QField passthrough
        @test undo_average(a) === a

        # Catch-all passthrough
        @test undo_average(:foo) === :foo

        # Num round-trip: Num(avg(a)) → undo → recovers the QSym directly
        # (Num only wraps BasicSymbolic, so when inner result is QField it's returned unwrapped)
        avg_num = Symbolics.Num(average(a))
        result_num = undo_average(avg_num)
        @test result_num === a

        # Equation round-trip: simple avg terms un-average to QField,
        # which can't be stored in Equation (requires BasicSymbolic), so returns Pair
        eq = Symbolics.Equation(average(a), average(ad))
        result_eq = undo_average(eq)
        @test result_eq isa Pair
        @test result_eq.first === a
        @test result_eq.second === ad

        # Equation with sum-of-averages: undo_average produces QAdd (not BasicSymbolic)
        avg_sum_lhs = average(a) + average(ad)
        avg_sum_rhs = average(a) + average(ad)
        eq2 = Symbolics.Equation(avg_sum_lhs, avg_sum_rhs)
        result_eq2 = undo_average(eq2)
        @test result_eq2 isa Pair  # QAdd can't be stored in Equation
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

        # QMul
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
