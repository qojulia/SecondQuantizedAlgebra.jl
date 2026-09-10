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

    @test zero(cA) == zero(QExpr)
    @test one(cA) == one(QExpr)
    @test +cA === cA
    @test iszero(cA * 0)

    @test 2 * cA + 3 * cA == 5 * cA
    @test iszero(cA - cA)
    @test cA + sA == sA + cA
    @test cA + a == a + cA
    @test 2 + cA == cA + 2
    @test one(QExpr) + 1 == 2 * one(QExpr)
    @test cA / 2 == (1 // 2) * cA
    @test cA // 2 == (1 // 2) * cA

    @test a * cA != cA * a
    @test (a * cA) * a == a * (cA * a)
    @test cA * A isa QExpr
    @test A * (A * cA) == (A * A) * cA

    scalar_qadd = 2 * commutator(a, a')
    @test cA * scalar_qadd == 2 * cA

    @test cA^0 == one(QExpr)
    @test cA^1 == cA
    @test cA^3 == cA * cA * cA
    @test_throws ArgumentError cA^(-1)
    n = -2
    @test_throws ArgumentError cA^n
    @test hash(cA) == hash(cos(A))

    @test @inferred(expim(A)) isa QExpr
    @test @inferred(expim(θ * A)) isa QExpr
    @test_throws ArgumentError expim(a)

    @test adjoint(adjoint(cA)) == cA
    @test adjoint(adjoint(sA)) == sA
    @test adjoint(adjoint(expim(A))) == expim(A)

    formal_sum = cA + im * sA
    @test adjoint(formal_sum) == adjoint(cA) + adjoint(im * sA)
    formal_product = a * cA * a
    @test adjoint(formal_product) == a' * adjoint(cA) * a'

    @test iszero(sin(zero(QExpr)))
    @test isone(cos(zero(QExpr)))
    @test isone(expim(zero(QExpr)))
end
