using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: QExpr

@testset "Formal operator expression canonical form" begin
    h = FockSpace(:a) ⊗ FockSpace(:b)
    a = Destroy(h, :a, 1)
    b = Destroy(h, :b, 2)
    x = a + a'
    y = b + b'

    @testset "Polynomial residue cancels against polynomials" begin
        residue = (sin(x) + a) - sin(x)
        @test residue isa QExpr
        @test iszero(residue - a)
        @test iszero(a - residue)
        @test residue + a == 2 * ((sin(x) + a) - sin(x))
        @test 2 * residue == (sin(x) + 2a) - sin(x)
        @test hash(2 * residue) == hash((sin(x) + 2a) - sin(x))
        @test iszero(residue * 0)
    end

    @testset "Scalars have one home in a product" begin
        @variables θ::Real
        @test (θ * a) * cos(x) == θ * (a * cos(x))
        @test (θ * x) * cos(y) == θ * (x * cos(y))
        @test (2a + 4a') * cos(x) == 2 * ((a + 2a') * cos(x))
        @test iszero((θ * a) * cos(x) - θ * (a * cos(x)))
    end

    @testset "Factors on disjoint spaces commute" begin
        @test sin(x) * b == b * sin(x)
        @test sin(x) * cos(y) == cos(y) * sin(x)
        @test iszero(commutator(b' * b, cos(x)))
        @test iszero(commutator(sin(x), cos(y)))
        @test iszero(commutator(sin(x), y))
        @test a * cos(y) * a' == a * a' * cos(y)
        @test b * sin(x) * b' * sin(x) == b * b' * sin(x) * sin(x)
    end

    @testset "Factors on a shared space keep their order" begin
        @test a * cos(x) != cos(x) * a
        @test !iszero(commutator(sin(x), x))
        @test sin(x) * cos(x) != cos(x) * sin(x)
        @test (a + b) * cos(x) != cos(x) * (a + b)
    end

    @testset "Reordering is consistent with lowering" begin
        @variables θ::Real
        lhs = θ * b' * sin(x) * b + cos(y) * a
        rhs = θ * sin(x) * b' * b + a * cos(y)
        @test lhs == rhs
        @test taylor(lhs, 0:3) == taylor(rhs, 0:3)
        @test taylor(sin(x) * cos(y), 0:3) == taylor(sin(x), 0:3) * taylor(cos(y), 0:3)
    end
end
