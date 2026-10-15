using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test

import SecondQuantizedAlgebra: QExpr, expim

@testset "Formal operator expression structural operations" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h) b::Destroy(h)
    @variables θ::Real ϕ::Real

    A = a + a'
    B = b + b'
    expr = θ * cos(A) + sin(A)

    substituted = @inferred substitute(expr, Dict(a => b))
    @test substituted == θ * cos(B) + sin(B)
    @test @inferred(substitute(expim(A), Dict(a => b))) == expim(B)

    @test get_operators(expr) == get_operators(A)
    @test isequal(get_variables(expr), get_variables(θ * A))
    @test acts_on(expr) == acts_on(A)
    @test isempty(get_indices(expr))

    vars_expr = θ * cos(ϕ * A) + sin(A)
    expected_vars = Set(Symbolics.unwrap.([θ, ϕ]))
    @test Set(get_variables(vars_expr)) == expected_vars
    @test Set(get_variables(vars_expr, [ϕ])) == Set([Symbolics.unwrap(ϕ)])

    buffer = Set{Any}()
    @test Symbolics.get_variables!(buffer, vars_expr) === buffer
    @test buffer == expected_vars

    filtered_buffer = Set{Any}()
    @test Symbolics.get_variables!(filtered_buffer, vars_expr, [ϕ]) === filtered_buffer
    @test filtered_buffer == Set([Symbolics.unwrap(ϕ)])

    i = Index(h, :i, 3, h)
    j = Index(h, :j, 3, h)
    ai = IndexedOperator(a, i)
    aj = IndexedOperator(a, j)
    indexed = cos(ai + ai')
    @test get_indices(indexed) == [i]
    @test change_index(indexed, i, j) == cos(aj + aj')

    swapped = cos(ai + aj')
    @test change_index(swapped, Dict(i => j, j => i)) == cos(aj + ai')

    @test normal_order(cos(a * a')) == cos(normal_order(a * a'))
    @test simplify(expr) == θ * cos(A) + sin(A)
    @test expand(expr) == θ * cos(A) + sin(A)

    formal_product = cos(A) * sin(B)
    @test normal_order(formal_product) == formal_product
    @test simplify(formal_product) == formal_product
    @test expand(formal_product) == formal_product

    @test @inferred(commutator(cos(A), a)) == cos(A) * a - a * cos(A)
    @test @inferred(commutator(a, cos(A))) == a * cos(A) - cos(A) * a
    @test @inferred(commutator(cos(A), sin(B))) ==
        cos(A) * sin(B) - sin(B) * cos(A)
    @test @inferred(anticommutator(cos(A), a)) == cos(A) * a + a * cos(A)
    @test @inferred(commutator(cos(A), 2)) isa QExpr
    @test iszero(commutator(cos(A), 2))
    @test iszero(commutator(2, cos(A)))
end

@testset "Exact structural rewrites through formal functions" begin
    h = NLevelSpace(:atom, 2)
    σ11 = Transition(h, :σ, 1, 1)
    @test expand_completeness(cos(σ11)) == cos(expand_completeness(σ11))

    j = Index(h, :j, 3, h)
    k = Index(h, :k, 3, h)
    σ(i, m, idx) = IndexedOperator(Transition(h, :σ, i, m), idx)
    polynomial = σ(2, 1, k) * σ(1, 2, j)
    @test assume_distinct_index(cos(polynomial), [(j, k)]) ==
        cos(assume_distinct_index(polynomial, [(j, k)]))
end
