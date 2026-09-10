using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: QAdd, QExpr, expim

@testset "Formal operator expressions" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    @variables θ::Real

    A = a + a'
    cA = @inferred cos(A)
    sA = @inferred sin(A)

    @test cA isa QExpr
    @test sA isa QExpr
    @test @inferred((a + a')^4) isa QAdd

    @test 2 * cA + 3 * cA == 5 * cA
    @test iszero(cA - cA)
    @test cA + sA == sA + cA
    @test one(QExpr) + 1 == 2 * one(QExpr)

    @test a * cA != cA * a
    @test (a * cA) * a == a * (cA * a)

    @test @inferred(expim(A)) isa QExpr
    @test @inferred(expim(θ * A)) isa QExpr
    @test_throws ArgumentError expim(a)

    @test adjoint(adjoint(cA)) == cA
    @test adjoint(adjoint(sA)) == sA
    @test adjoint(adjoint(expim(A))) == expim(A)

    @test iszero(sin(zero(QExpr)))
    @test isone(cos(zero(QExpr)))
    @test isone(expim(zero(QExpr)))
end
