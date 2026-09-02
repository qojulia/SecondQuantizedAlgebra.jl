using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils, SymReal, symtype
using Symbolics: Symbolics, @variables
using Test
import SecondQuantizedAlgebra: get_sum_body, get_sum_indices, get_sum_non_equal,
    has_sum_metadata, indexed_sum, QAdd

@testset "Expectation-value API" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = adjoint(a)

    @testset "construction and linearity" begin
        avg_a = average(a)
        avg_ad = average(ad)
        @test avg_a isa SymbolicUtils.BasicSymbolic{SymReal}
        @test symtype(avg_a) === Number
        @test is_average(avg_a)
        @test !is_average(a)
        @test !is_average(3)
        @test !is_average(nothing)
        @test !isequal(avg_a, avg_ad)

        @test isequal(average(3), 3)
        @test isequal(average(0), 0)
        @test iszero(
            undo_average(average(a + ad)) -
                undo_average(average(a) + average(ad)),
        )
        @test iszero(
            undo_average(average(3 * ad * a)) -
                undo_average(3 * average(ad * a)),
        )

        @variables κ::Real
        @test isequal(average(κ), κ)
        @test isequal(average(average(a)), average(a))
        @test isequal(average(3 * ad * a), 3 * average(ad * a))
        @test isequal(average(κ * a), κ * average(a))
        @test iszero(undo_average(average(κ * a)) - κ * a)
    end

    @testset "round trips" begin
        expressions = (a, ad, ad * a, 3 * ad * a, (2 + 3im) * a)
        for expression in expressions
            @test iszero(undo_average(average(expression)) - expression)
        end

        @test undo_average(average(a)) isa QAdd
        @test undo_average(3) isa QAdd

        avg_sum = average(a) + average(ad)
        @test iszero(undo_average(avg_sum) - (a + ad))
        @test iszero(undo_average(average(a)) - a)
        @test iszero(undo_average(3) - 3)

        equation = Symbolics.Equation(average(a), average(ad))
        restored = undo_average(equation)
        @test restored isa Pair
        @test iszero(restored.first - a)
        @test iszero(restored.second - ad)
    end

    @testset "indexed sums retain scope" begin
        ha = NLevelSpace(:atom, 2)
        hmix = FockSpace(:field) ⊗ ha
        field = Destroy(hmix, :a, 1)
        @variables N
        i = Index(hmix, :i, N, ha)
        j = Index(hmix, :j, N, ha)
        σ12(k) = IndexedOperator(Transition(hmix, :σ, 1, 2, 2), k)

        sum_expr = Σ(field' * σ12(i), i)
        avg_sum = average(sum_expr)
        @test is_indexed_sum(avg_sum)
        @test has_sum_metadata(avg_sum)
        @test get_sum_indices(avg_sum) == [i]
        @test isempty(get_sum_non_equal(avg_sum))
        @test iszero(undo_average(avg_sum) - sum_expr)

        body = get_sum_body(avg_sum)
        @test isequal(indexed_sum(body, [i]), avg_sum)
        @test isequal(indexed_sum(Symbolics.wrap(body), [i]), avg_sum)

        off_diagonal = average(Σ(σ12(i) * σ12(j)', i, [j]))
        @test is_indexed_sum(off_diagonal)
        @test get_sum_indices(off_diagonal) == [i]
        @test get_sum_non_equal(off_diagonal) == [(i, j)]

        plain = average(field + field')
        @test !is_indexed_sum(plain)
        @test !has_sum_metadata(plain)
        @test !has_sum_metadata(3)
    end

    @testset "scope is part of the public symbolic value" begin
        ha = NLevelSpace(:atom, 2)
        @variables N
        i = Index(ha, :i, N, ha)
        j = Index(ha, :j, N, ha)
        σ(k) = IndexedOperator(Transition(ha, :σ, 1, 2), k)

        with_scope = average(Σ(σ(i), i))
        with_constraint = average(Σ(σ(i), i, [j]))
        @test !isequal(with_scope, with_constraint)
        @test !isequal(Symbolics.simplify(with_scope - with_constraint), 0)

        rebuilt = indexed_sum(get_sum_body(with_scope), [i]; non_equal = [(i, j)])
        @test get_sum_indices(rebuilt) == [i]
        @test get_sum_non_equal(rebuilt) == [(i, j)]
        @test !isequal(rebuilt, with_scope)
    end

    @testset "indexed coefficients stay inside their sum" begin
        ha = NLevelSpace(:atom, 2)
        @variables N
        i = Index(ha, :i, N, ha)
        j = Index(ha, :j, N, ha)
        g(k) = IndexedVariable(:g, k)
        u(k, l) = DoubleIndexedVariable(:u, k, l)
        σ(k) = IndexedOperator(Transition(ha, :σ, 1, 2), k)

        coefficient_sum = Σ(u(j, i), i)
        averaged_coefficient_sum = average(coefficient_sum)
        @test is_indexed_sum(averaged_coefficient_sum)
        @test get_sum_indices(averaged_coefficient_sum) == [i]
        @test iszero(undo_average(averaged_coefficient_sum) - coefficient_sum)

        expression = Σ(g(i) * σ(i), i)
        averaged = average(expression)
        @test is_indexed_sum(averaged)
        @test get_sum_indices(averaged) == [i]
        @test iszero(undo_average(averaged) - expression)

        nested = Σ(g(i) * g(j) * σ(i) * σ(j), i, j)
        @test iszero(undo_average(average(nested)) - nested)
    end

    @testset "averaging nontrivial operator expressions" begin
        @variables Δ γ η
        h2 = FockSpace(:double_average)
        a2 = Destroy(h2, :a)
        H = Δ * a2' * a2 + γ * (a2' + a2) + η * (a2' * a2 * a2 + a2' * a2' * a2)
        rhs = commutator(im * H, a2' * a2)
        averaged = average(rhs)

        # A complex, multi-term QAdd must survive the average layer without
        # losing terms or changing their coefficients.
        @test iszero(simplify(undo_average(averaged) - rhs))
        @test iszero(
            simplify(
                undo_average(average(rhs + a2)) - (rhs + a2),
            ),
        )
    end

    @testset "averages across multiple Hilbert-space factors" begin
        ha = NLevelSpace(:average_atom, 2)
        hmix = h ⊗ ha
        field = Destroy(hmix, :field, 1)
        transition = Transition(hmix, :σ, 1, 2, 2)

        @test is_average(average(field' * transition))
        @test iszero(undo_average(average(field' * transition)) - field' * transition)
        @test iszero(
            undo_average(average(field + transition)) - (field + transition),
        )
    end

    @testset "acts_on and get_indices" begin
        ha = NLevelSpace(:atom, 2)
        hmix = h ⊗ ha
        field = Destroy(hmix, :a, 1)
        σ12 = Transition(hmix, :σ, 1, 2, 2)
        @test acts_on(field) == [1]
        @test acts_on(σ12) == [2]
        @test acts_on(field' * σ12) == [1, 2]
        @test acts_on(average(field' * σ12)) == [1, 2]

        @variables N
        i = Index(hmix, :i, N, ha)
        j = Index(hmix, :j, N, ha)
        σi = IndexedOperator(σ12, i)
        σj = IndexedOperator(σ12', j)
        @test Set(get_indices(σi + σj)) == Set([i, j])
        @test Set(get_indices(average(σi) + average(σj))) == Set([i, j])
    end

    @testset "Hermitian values expose their mathematical type" begin
        @test symtype(average(ad * a)) === Real
        @test symtype(average(ad * ad * a * a)) === Real
        @test symtype(average(a)) === Number
        @test symtype(average(ad * a * a)) === Number

        hn = NLevelSpace(:atom, 3)
        @test symtype(average(Transition(hn, :σ, 2, 2))) === Real
        @test symtype(average(Transition(hn, :σ, 1, 2))) === Number

        hp = PauliSpace(:pauli)
        hs = SpinSpace(:spin)
        hq = PhaseSpace(:phase)
        @test symtype(average(Pauli(hp, :σ, 1))) === Real
        @test symtype(average(Spin(hs, :S, 1))) === Real
        @test symtype(average(Position(hq, :x))) === Real
        @test symtype(average(Momentum(hq, :p))) === Real

        @test symtype(average(2 * ad * a + 5 * ad * ad * a * a)) === Real
        @test symtype(average(2 * ad * a + 3im * a)) === Number

        @test isequal(conj(average(ad * a)), average(ad * a))
        @test !isequal(conj(average(a)), average(a))
    end

    @testset "indexed sums inherit the body type" begin
        i = Index(h, :i, 3, h)
        @variables g::Real
        summed_number = average(Σ(IndexedOperator(a, i), i))
        summed_real = average(
            Σ(g * IndexedOperator(ad, i) * IndexedOperator(a, i), i),
        )
        @test is_indexed_sum(summed_number)
        @test is_indexed_sum(summed_real)
        @test symtype(summed_number) === Number
        @test symtype(summed_real) === Real

        rebuilt = Symbolics.substitute(summed_real, Dict(g => 2g))
        @test is_indexed_sum(rebuilt)
        @test symtype(rebuilt) === Real
    end

    @testset "time-dependent averages" begin
        @variables t
        time = SymbolicUtils.unwrap(t)
        hermitian = make_time_dependent(average(ad * a), time)
        nonhermitian = make_time_dependent(average(a), time)
        @test symtype(SymbolicUtils.unwrap(hermitian)) === Real
        @test symtype(SymbolicUtils.unwrap(nonhermitian)) === Number
        @test iszero(
            undo_average(average((2 + 3im) * a)) - (2 + 3im) * a
        )
    end

    @testset "public entry points remain inferable" begin
        @inferred average(a)
        @inferred average(ad * a)
        @inferred undo_average(average(a))
        i = Index(h, :i, 3, h)
        summed = average(Σ(IndexedOperator(a, i), i))
        @inferred get_sum_indices(summed)
        @inferred get_sum_non_equal(summed)
    end
end
