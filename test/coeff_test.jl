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
        @test occursin("1//4", string(real(coefficient((1 // 4) * a))))
        @test string(real(coefficient(0.25 * a))) == "0.25"
        @test iszero(simplify((1 // 4) * a + (1 // 4) * a - (1 // 2) * a))
        @test iszero(simplify((1 // 3) * a - (1 // 3) * a))
    end

    @testset "large and rational exactness boundaries" begin
        @variables r
        @test occursin("1//2", string(real(coefficient(r / 2 * a))))
        @test occursin("1//4", string(real(coefficient((1 // 4) * r * a))))
        large2 = big(2)^70
        @test coefficient(large2 * a) == large2
        @test isequal(real(coefficient((large2 + 1) * a)), Num(large2 + 1))
        @test string(real(coefficient(2))) == "2"
        @test string(real(coefficient(0.7))) == "0.7"
        # Rational that round-trips through Float64 stays exact via prefactor string
        @test isequal(coefficient(Complex(1 // 3, 0)), Complex(Num(1 // 3), Num(0)))
        @test iszero(simplify(((1 // 3) * a) * 3 - a))
    end

    @testset "conjugation and hash stability are observable" begin
        c = coefficient(2)
        @test isequal(conj(c), c)
        @test hash(conj(c)) == hash(c)
        doubled = conj(conj(coefficient(2 + 3im)))
        @test isequal(doubled, coefficient(2 + 3im))
        @test hash(doubled) == hash(coefficient(2 + 3im))
        @test isequal(conj(coefficient(2 + 3im)), coefficient(2 - 3im))
        q = 2 * a' * a
        @test isequal(q', q)
        @test hash(q') == hash(q)
        @variables g k
        simplified = simplify(((g + g * k) / (1 + k)) * a)
        @test isequal(simplified, g * a)
        @test hash(simplified) == hash(g * a)
    end

    @testset "wide sums fold and dedup through public algebra" begin
        @variables p[1:6] ω g κ
        distinct = [Num(k) * p[k] for k in 1:6]
        s = sum(distinct[k] * a for k in 1:6)
        @test length(s) == 1
        @test occursin("p[1]", string(prefactor(s)))
        @test occursin("p[6]", string(prefactor(s)))
        coalesced = g * a + 2g * a + 3g * a
        @test iszero(simplify(coalesced - 6g * a))
        @test length(coalesced) == 1
        different = g * a + κ * a
        @test length(different) == 1
        @test occursin("g", string(prefactor(different)))
        @test occursin("κ", string(prefactor(different)))
        @test iszero(simplify(different - (g + κ) * a))
        product = (2g * a) * (3κ * a')
        @test iszero(simplify(product - 6 * g * κ * a * a'))
        cancelling = g * a + κ * a - g * a - κ * a
        @test iszero(cancelling)
        @variables q
        two_ops = p[1] * a + p[2] * a'
        @test length(two_ops) == 2
    end

    @testset "radicals and exact powers stay canonical" begin
        @variables ω g κ
        @test iszero(simplify(sqrt(g) * sqrt(g) - g))
        @test iszero(simplify(g^(1 // 2) * g^(1 // 2) - g))
        @test iszero(simplify(sqrt(4) * a - 2 * a))
        @test iszero(simplify(g^2 * a - g * g * a))
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

    @testset "public phase and scalar boundary cases" begin
        @variables θ φ ω t z::Number
        @test isequal(exponential_form(exp(im * θ)), expim(θ))
        @test isequal(exponential_form(cis(θ)), expim(θ))
        @test isequal(exponential_form(conj(expim(θ))), expim(-θ))
        @test exponential_form(abs(expim(θ))) == 1
        @test exponential_form(abs2(expim(θ))) == 1
        @test isequal(exponential_form(expim(θ)^2), expim(2θ))
        @test isequal(exponential_form(expim(θ)^(-2)), expim(-2θ))

        product = expim(θ) * expim(φ) * exp(ω * t)
        quotient = expim(θ) / (expim(φ) * exp(ω * t))
        @test iszero(simplify(exponential_form(product) - product))
        @test iszero(simplify(exponential_form(quotient) - quotient))
        @test iszero(simplify(trigonometric_form(product) - trigonometric_form(exponential_form(product))))
        @test length(phase_terms(product)) == 1
        @test length(phase_terms(expim(θ) / 2)) == 1
        @test length(phase_terms(conj(expim(θ)))) == 1
        @test length(phase_terms(expim(θ)^(-1))) == 1

        trig = trigonometric_form(Complex(Num(z), Num(z)))
        @test iszero(simplify(trig * a - (z + im * z) * a))
    end
end
