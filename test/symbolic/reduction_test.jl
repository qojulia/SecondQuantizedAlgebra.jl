using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test
import SecondQuantizedAlgebra: exponential_form, trigonometric_form

@testset "Public coefficient reduction" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables θ φ r ω t g

    @testset "trigonometric and hyperbolic identities" begin
        @test iszero(simplify((cos(θ)^2 + sin(θ)^2 - 1) * a))
        @test iszero(simplify((cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4 - 1) * a))
        @test iszero(simplify((cosh(r)^2 - sinh(r)^2 - 1) * a))
        @test iszero(simplify((cosh(r)^4 - 2cosh(r)^2 * sinh(r)^2 + sinh(r)^4 - 1) * a))
        @test iszero(
            simplify(
                (cos(θ)^2 * cosh(r)^2 + sin(θ)^2 * cosh(r)^2 - sinh(r)^2 - 1) * a,
            )
        )
    end

    @testset "independent arguments and products" begin
        expression = (cos(θ)^2 + sin(θ)^2) * (cos(φ)^2 + sin(φ)^2) * a
        @test iszero(simplify(expression - a))
        @test iszero(simplify((g + cos(θ)^2 + sin(θ)^2 - g - 1) * a))
        @test !iszero(simplify((cos(θ)^2 + sin(φ)^2 - 1) * a))
    end

    @testset "cancellation gate — public view" begin
        @test isequal(simplify(cos(θ)^3 * sin(θ) * a), cos(θ)^3 * sin(θ) * a)
        @test isequal(simplify((cos(θ)^4 + sin(θ)^4) * a), (cos(θ)^4 + sin(θ)^4) * a)
        @test !iszero(simplify((cos(θ)^4 + sin(θ)^4 - 1) * a))
        @test iszero(simplify((cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4 - 1) * a))
    end

    @testset "composite arguments and substitutions" begin
        identity = cos(ω * t)^2 + sin(ω * t)^2
        @test iszero(simplify((identity - 1) * a))
        @test iszero(simplify(substitute((identity - 1) * a, Dict(t => 2))))
        @test iszero(simplify((cos(θ) * sin(θ) - sin(2θ) / 2) * a))
        @test isequal(simplify(g * a), g * a)
    end

    @testset "zero coefficients disappear from expressions" begin
        expression = (cos(θ)^4 + 2cos(θ)^2 * sin(θ)^2 + sin(θ)^4 - 1) * a' * a
        @test iszero(simplify(expression))
        @test isequal(simplify(expression + 2 * a), simplify(2 * a))
    end

    @testset "public coefficient views" begin
        value = cos(ω * t) + sin(ω * t)
        exponential = exponential_form(value)
        @test isequal(trigonometric_form(exponential), value)
        @test isequal(trigonometric_form(exponential_form(value * a)), value * a)
        @inferred simplify(value * a)
        @inferred exponential_form(value * a)
        @inferred trigonometric_form(exponential)
    end

    @testset "simplification remains idempotent" begin
        expression = (cos(θ)^2 + sin(θ)^2) * (g + a' * a)
        @test isequal(simplify(simplify(expression)), simplify(expression))
        @test isequal(Symbolics.expand(simplify(expression)), simplify(expression))
    end
end
