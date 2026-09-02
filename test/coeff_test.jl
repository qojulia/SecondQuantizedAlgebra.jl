using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using Test
import SecondQuantizedAlgebra: expim, exponential_form, phase_terms, to_num,
    trigonometric_form

@testset "Symbolic coefficients and phases" begin
    h = FockSpace(:coefficients)
    a = Destroy(h, :a)

    coefficient(x) = prefactor(x * a)

    @testset "numeric and exact coefficients" begin
        @test coefficient(2) == 2
        @test coefficient(2 + 3im) == 2 + 3im
        @test isequal(coefficient(1 // 4), Complex(Num(1 // 4), Num(0)))
        @test isequal(
            coefficient(Complex(1 // 4, 1 // 2)),
            Complex(Num(1 // 4), Num(1 // 2))
        )
        @test string(real(coefficient(2))) == "2"
        @test occursin("1//4", string(real(coefficient(1 // 4))))

        large = Complex((2^53 + 1) // 1, 0 // 1)
        @test isequal(real(coefficient(large)), Num(2^53 + 1))
        @test isequal(imag(coefficient(large)), Num(0))
    end

    @testset "symbolic arithmetic stays faithful" begin
        @variables g κ r
        @test isequal(coefficient(g), Complex(Num(g), Num(0)))
        @test isequal(coefficient(im * g), Complex(Num(0), Num(g)))
        @test isequal(coefficient(g + im * κ), Complex(Num(g), Num(κ)))
        @test isequal(coefficient(g * κ), Complex(Num(g * κ), Num(0)))
        @test isequal(coefficient(g + g), Complex(Num(2g), Num(0)))
        @test isequal(coefficient(g * κ / g), Complex(Num(κ), Num(0)))
        @test iszero(simplify((g - g) * a))
        @test isequal(conj(coefficient(2 + 3im)), coefficient(2 - 3im))
        @test isequal(conj(conj(coefficient(g + im * κ))), coefficient(g + im * κ))

        expression = (g + g * r) / (1 + r) * a
        @test iszero(simplify(expression - g * a))
    end

    @testset "public arithmetic collects and cancels coefficients" begin
        @variables γ δ β
        opposite = (γ / δ) * a + (-γ / δ) * a
        collected = (γ / δ) * a + (β / δ) * a

        @test iszero(simplify(opposite))
        @test iszero(simplify(collected - ((γ + β) / δ) * a))
        @test !iszero(simplify(collected))
    end

    @testset "phase group" begin
        @variables ω θ t x y
        p = expim(ω * t)

        @test p * expim(-ω * t) == 1
        @test expim(x) * expim(y) == expim(x + y)
        @test expim(x) * expim(y - x) == expim(y)
        @test expim(x)^3 == expim(3x)
        @test expim(x)^(-3) == expim(-3x)
        @test inv(p) == conj(p) == expim(-ω * t)
        @test abs2(p) == 1
        @test abs(p) == 1
        @test isequal(real(expim(θ)), cos(θ))
        @test isequal(imag(expim(θ)), sin(θ))
        @test isequal(real(conj(expim(θ))), cos(θ))
        @test isequal(imag(conj(expim(θ))), -sin(θ))
        @test (@inferred expim(0.5)) ≈ exp(0.5im)
        @test_throws ArgumentError expim(1 + 0im)
        @test_throws ArgumentError expim(1 + 2im)
    end

    @testset "phase substitution and calculus" begin
        @variables ω t θ
        p = expim(ω * t)
        @test substitute(p, Dict()) === p
        @test substitute(p, Dict(ω => 2ω)) == expim(2ω * t)
        @test substitute(p, Dict(t => 0)) == 1
        @test Symbolics.simplify(Symbolics.derivative(p, t) - im * ω * p) == 0

        phase = expim(θ)
        @test Symbolics.simplify(Symbolics.derivative(phase, θ) - im * expim(θ)) == 0
        @test isequal(to_num(phase), to_num(expim(θ)))
    end

    @testset "explicit representations round-trip" begin
        @variables θ ω t g κ
        cosine = exponential_form(cos(θ))
        sine = exponential_form(sin(θ))
        composite = exponential_form(cos(ω * t) + sin(ω * t))

        @test cosine == (expim(θ) + expim(-θ)) / 2
        @test sine == (expim(θ) - expim(-θ)) / (2im)
        @test trigonometric_form(cosine) == cos(θ)
        @test trigonometric_form(sine) == sin(θ)
        @test trigonometric_form(expim(θ)) == cos(θ) + im * sin(θ)
        @test trigonometric_form(composite) == cos(ω * t) + sin(ω * t)
        @test iszero(
            simplify(
                exponential_form(cos(ω * t) * a) -
                    (expim(ω * t) + expim(-ω * t)) / 2 * a
            )
        )

        reconstruct(c) = sum(
            term.amplitude * expim(term.phase) for term in phase_terms(c);
            init = 0,
        )
        for value in (cosine, sine, exponential_form(cos(ω * t) + g * sin(κ * t)))
            @test iszero(simplify(reconstruct(value) - value))
        end

        @test_throws ArgumentError phase_terms(exponential_form(1 / (1 + expim(θ))))
    end

    @testset "phases inside operator expressions" begin
        @variables ω t g
        p = expim(ω * t)
        expression = g * p * a + conj(p) * a'
        @test iszero(simplify(expression - (g * expim(ω * t) * a + expim(-ω * t) * a')))
        @test iszero(
            simplify(
                trigonometric_form(exponential_form(expression)) -
                    trigonometric_form(expression)
            )
        )
    end

    @testset "public coefficient boundaries remain inferable" begin
        @variables ω t
        @test prefactor(expim(ω * t) * a) isa Complex{Num}
        @inferred expim(ω * t)
        @inferred exponential_form(cos(ω * t))
        @inferred trigonometric_form(expim(ω * t))
        @inferred phase_terms(expim(ω * t))
    end
end
