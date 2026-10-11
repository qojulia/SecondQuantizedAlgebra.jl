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
        @test (θ * a + θ^2 * a') * cos(y) == θ * ((a + θ * a') * cos(y))
        @test (2θ * a + 4θ * a') * cos(y) == 2θ * ((a + 2a') * cos(y))
        @test ((θ + θ^2) * a + θ * a') * cos(y) == θ * (((1 + θ) * a + a') * cos(y))
        @test iszero((θ * a + θ^2 * a') * cos(y) - θ * ((a + θ * a') * cos(y)))
        # No common parameter: nothing is divided out.
        @test ((θ * a + a') * cos(y)).coeff == 1
    end

    @testset "Bound-index content stays inside its sum" begin
        hs = h ⊗ NLevelSpace(:s, 2)
        as = Destroy(hs, :a, 1)
        i = Index(hs, :i, 3, 3)
        σ = IndexedOperator(Transition(hs, :σ, 1, 2, 3), i)
        g = IndexedVariable(:g, i)
        expr = Σ(g * σ, i) * cos(as + as')
        @test expr.coeff == 1
        @test taylor(expr, 0:0) == Σ(g * σ, i)
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

    @testset "Polynomial merging preserves tensor commutation" begin
        C = cos(x)
        for local_factor in (a, a' * a, x, 2a + 4a'),
                remote_factor in (b, b' * b, y, 3b + b')
            lhs = (C * remote_factor) * local_factor
            rhs = C * (remote_factor * local_factor)
            @test lhs == rhs
            @test hash(lhs) == hash(rhs)
            @test iszero(lhs - rhs)
            @test (local_factor * remote_factor) * C == local_factor * (remote_factor * C)
            @test taylor(lhs, 0:2) == taylor(rhs, 0:2)
        end
        @test (C * b) * (a * a') == C * (b * (a * a'))
        @test C * (a + b) != (a + b) * C
        h4 = h ⊗ FockSpace(:c) ⊗ FockSpace(:d)
        c = Destroy(h4, :c, 3)
        d = Destroy(h4, :d, 4)
        local_pair = a + b
        remote_pair = c + d
        @test (C * remote_pair) * local_pair == C * (remote_pair * local_pair)
        @test iszero((C * remote_pair) * local_pair - C * (remote_pair * local_pair))
    end

    @testset "Floating content normalization preserves finite coefficients" begin
        for scale in (1.0e-320, 1.0e-308, 1.0e-200, 1.0e200),
                phase in (1, im)
            P = (phase * scale) * a + (phase * (2scale)) * a'
            @test taylor(P * cos(x), 0:0) == P
            @test taylor(cos(x) * P, 0:0) == P
        end
        P = 1.0e-320 * a + 1.0e300 * a'
        @test taylor(P * cos(x), 0:0) == P
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
