using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: QExpr, expim

@testset "Formal operator expression structural operations" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h) b::Destroy(h)
    @variables θ::Real

    A = a + a'
    B = b + b'
    expr = θ * cos(A) + sin(A)

    substituted = @inferred substitute(expr, Dict(a => b))
    @test substituted == θ * cos(B) + sin(B)
    @test @inferred(substitute(expim(A), Dict(a => b))) == expim(B)

    @test get_operators(expr) == get_operators(A)
    @test get_variables(expr) == get_variables(θ * A)
    @test acts_on(expr) == acts_on(A)
    @test isempty(get_indices(expr))

    i = Index(h, :i, 3, h)
    j = Index(h, :j, 3, h)
    ai = IndexedOperator(a, i)
    aj = IndexedOperator(a, j)
    indexed = cos(ai + ai')
    @test get_indices(indexed) == [i]
    @test change_index(indexed, i, j) == cos(aj + aj')

    @test normal_order(cos(a * a')) == cos(normal_order(a * a'))
    @test simplify(expr) == θ * cos(A) + sin(A)
    @test expand(expr) == θ * cos(A) + sin(A)

    @test @inferred(commutator(cos(A), a)) == cos(A) * a - a * cos(A)
    @test @inferred(commutator(a, cos(A))) == a * cos(A) - cos(A) * a
    @test @inferred(commutator(cos(A), sin(B))) ==
        cos(A) * sin(B) - sin(B) * cos(A)
end
