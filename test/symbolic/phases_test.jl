using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
using Test
import SecondQuantizedAlgebra: expim, exponential_form, phase_terms, to_num,
    trigonometric_form

@testset "Symbolic phases" begin
    h = FockSpace(:coefficients)
    a = Destroy(h, :a)

    stored_coefficient(x) = only(x).second

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
        @test isequal(expim(SymbolicUtils.unwrap(θ)), expim(θ))
        @test isequal(expim(θ + ω)^2, expim(2θ + 2ω))
        @test isequal(expim(θ + ω)^(-2), expim(-2θ - 2ω))
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
        @variables θ φ ω t z::Number g
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
        complex_coefficient = Complex(Num(z), Num(z)) * a
        @test isequal(prefactor(complex_coefficient), Complex(Num(z), Num(z)))
        @test isequal(
            exponential_form(Complex(Num(z), Num(z))),
            Complex(Num(z), Num(z)),
        )
        @test isequal(real(-expim(θ)), -cos(θ))
        @test isequal(imag(-expim(θ)), -sin(θ))
        unit_phase = (3 // 5 + (4 // 5) * im) * expim(θ)
        @test iszero(
            simplify(real(unit_phase) - (3 // 5) * cos(θ) + (4 // 5) * sin(θ)),
        )
        @test iszero(
            simplify(imag(unit_phase) - (4 // 5) * cos(θ) - (3 // 5) * sin(θ)),
        )
        @test isempty(phase_terms(exponential_form(Num(0))))

        radical = exponential_form(sqrt(g))
        @test isequal(radical, sqrt(g))
        @test isequal(exponential_form(cbrt(g)), cbrt(g))
    end

    @testset "public iterator exposes coefficient arithmetic" begin
        @variables θ φ ω t g κ z::Number
        poly = stored_coefficient((g + κ) * a)
        raw = stored_coefficient(Complex(Num(z), Num(z)) * a)
        phase_product = stored_coefficient(expim(θ) * exp(ω * t) * a)
        phase_quotient = stored_coefficient(
            expim(θ) / (expim(φ) * exp(ω * t)) * a,
        )

        @test poly + 2 == stored_coefficient((g + κ + 2) * a)
        @test 2 + poly == stored_coefficient((2 + g + κ) * a)
        @test poly - 2 == stored_coefficient((g + κ - 2) * a)
        @test 2 - poly == stored_coefficient((2 - g - κ) * a)
        @test poly * 2 == stored_coefficient((2g + 2κ) * a)
        @test 2 * poly == stored_coefficient((2g + 2κ) * a)
        @test poly / 2 == stored_coefficient(((g + κ) / 2) * a)
        @test inv(poly) == stored_coefficient((1 / (g + κ)) * a)
        @test poly^3 == stored_coefficient((g + κ)^3 * a)
        @test poly^(-2) == stored_coefficient((g + κ)^(-2) * a)
        @test iszero(poly + -poly)

        native = stored_coefficient(2 * a)
        @test native == 2
        @test 2 == native
        @test native != 3
        @test native^3 == stored_coefficient(8 * a)
        @test native^(-2) == stored_coefficient((1 / 4) * a)

        large = stored_coefficient(complex(big(2)^70, big(3)^70) * a)
        @test isequal(real(large), Num(big(2)^70))
        @test isequal(imag(large), Num(big(3)^70))

        @test isequal(real(raw), z)
        @test isequal(imag(raw), z)
        @test isequal(real(conj(raw)), conj(z))
        @test isequal(imag(conj(raw)), -conj(z))
        @test iszero(raw + -raw)
        @test iszero(phase_product + -phase_product)
        @test length(phase_terms(phase_product)) == 1
        @test isequal(phase_terms(phase_product)[1].phase, θ)
        @test length(phase_terms(phase_quotient)) == 1
        @test isequal(phase_terms(phase_quotient)[1].phase, θ - φ)

        @test isequal(real(stored_coefficient(cos(ω * t) * a)), cos(ω * t))
        @test isequal(imag(stored_coefficient(cos(ω * t) * a)), 0)
        @test Symbolics.simplify(
            to_num(inv(stored_coefficient(cos(ω * t) * a))) - 1 / cos(ω * t),
        ) == 0
        @test isequal(
            to_num(Symbolics.derivative(raw, z)),
            Complex(Num(1), Num(1)),
        )
    end

    @testset "public raw coefficient paths" begin
        @variables θ φ g z::Number
        ordinary = stored_coefficient(exp(1 + g) * a)
        raw_phase = stored_coefficient(exp(1 + g) * expim(θ) * a)
        raw_phase_other = stored_coefficient(exp(2 + g) * expim(φ) * a)
        raw_phase_negative = stored_coefficient(exp(2 + g) * expim(-θ) * a)

        mixed_phase = stored_coefficient(g * expim(θ) * expim(φ) * a)
        mixed_terms = phase_terms(mixed_phase)
        @test length(mixed_terms) == 1
        @test Symbolics.simplify(mixed_terms[1].phase - θ - φ) == 0
        @test Symbolics.simplify(to_num(mixed_terms[1].amplitude) - g) == 0

        duplicated_phase = stored_coefficient(g * expim(θ) * a) *
            stored_coefficient(g * expim(φ) * a)
        duplicated_terms = phase_terms(duplicated_phase)
        @test length(duplicated_terms) == 1
        @test Symbolics.simplify(duplicated_terms[1].phase - θ - φ) == 0
        @test Symbolics.simplify(to_num(duplicated_terms[1].amplitude) - g^2) == 0

        constant_phase = expim(θ + 1) * expim(-θ)
        @test constant_phase == exp(im)
        mixed_constant_phase = stored_coefficient(
            g * expim(θ + 1) * expim(-θ) * a,
        )
        @test Symbolics.simplify(to_num(mixed_constant_phase) - g * exp(im)) == 0

        raw_sum = raw_phase + raw_phase_other
        sum_terms = phase_terms(raw_sum)
        @test length(sum_terms) == 2
        @test any(isequal(term.phase, θ) for term in sum_terms)
        @test any(isequal(term.phase, φ) for term in sum_terms)
        @test iszero(raw_sum + -raw_sum)
        @test hash(raw_phase) == hash(
            stored_coefficient(exp(1 + g) * expim(θ) * a),
        )

        square_terms = phase_terms(raw_phase^2)
        @test length(square_terms) == 1
        @test Symbolics.simplify(square_terms[1].phase - 2θ) == 0
        @test Symbolics.simplify(
            real(square_terms[1].amplitude) - exp(1 + g)^2,
        ) == 0

        single_quotient = ordinary / stored_coefficient(expim(θ) * a)
        expected_quotient = stored_coefficient(exp(1 + g) * expim(-θ) * a)
        @test Symbolics.simplify(
            to_num(single_quotient) - to_num(expected_quotient),
        ) == 0

        denominator_product = raw_phase * raw_phase_other
        denominator_terms = phase_terms(ordinary / denominator_product)
        @test length(denominator_terms) == 1
        @test Symbolics.simplify(denominator_terms[1].phase + θ + φ) == 0

        product_terms = phase_terms((raw_phase * raw_phase_other) / ordinary)
        @test length(product_terms) == 1
        @test Symbolics.simplify(product_terms[1].phase - θ - φ) == 0

        cancellation_terms = phase_terms(
            (raw_phase * raw_phase_negative) / ordinary,
        )
        @test length(cancellation_terms) == 1
        @test iszero(cancellation_terms[1].phase)
        @test Symbolics.simplify(
            real(cancellation_terms[1].amplitude) - exp(2 + g),
        ) == 0

        power_terms = phase_terms(raw_phase^2 / ordinary)
        @test length(power_terms) == 1
        @test Symbolics.simplify(power_terms[1].phase - 2θ) == 0

        raw_complex = stored_coefficient((2 + 3im) * exp(1 + g) * a)
        conjugate_complex = conj(raw_complex)
        @test Symbolics.simplify(real(conjugate_complex) - 2exp(1 + g)) == 0
        @test Symbolics.simplify(imag(conjugate_complex) + 3exp(1 + g)) == 0

        raw_nonreal = stored_coefficient(exp(1 + z) * a)
        conjugate_nonreal = conj(raw_nonreal)
        @test isequal(real(conjugate_nonreal), real(exp(1 + z)))
        @test isequal(imag(conjugate_nonreal), -imag(exp(1 + z)))

        raw_complex_slots = stored_coefficient(
            Complex(Num(exp(1 + g)), Num(exp(1 + g))) * a,
        )
        @test Symbolics.simplify(real(raw_complex_slots) - exp(1 + g)) == 0
        @test Symbolics.simplify(imag(raw_complex_slots) - exp(1 + g)) == 0

        complex_constant = stored_coefficient(Complex(1 // 3, 1 // 2) * a)
        complex_quotient = complex_constant / ordinary
        @test Symbolics.simplify(
            real(complex_quotient) - (1 // 3) / exp(1 + g),
        ) == 0
        @test Symbolics.simplify(
            imag(complex_quotient) - (1 // 2) / exp(1 + g),
        ) == 0
        @test Symbolics.simplify(
            to_num(Symbolics.derivative(complex_quotient, g)) -
                ((-1 // 3) - im * (1 // 2)) / exp(1 + g),
        ) == 0
        conjugate_quotient = conj(complex_quotient)
        @test Symbolics.simplify(
            real(conjugate_quotient) - (1 // 3) / exp(1 + g),
        ) == 0
        @test Symbolics.simplify(
            imag(conjugate_quotient) + (1 // 2) / exp(1 + g),
        ) == 0
        conjugate_constant = conj(complex_constant)
        @test isequal(real(conjugate_constant), Num(1 // 3))
        @test isequal(imag(conjugate_constant), -Num(1 // 2))

        cube_root = stored_coefficient(cbrt(g) * a)
        square_root = stored_coefficient(sqrt(g) * a)
        inverse_square_root = stored_coefficient(g^(-1 // 2) * a)
        @test isequal(real(cube_root), g^(1 // 3))
        @test isequal(real(square_root), sqrt(g))
        @test isequal(real(inverse_square_root), 1 / sqrt(g))
        negative_integer_power = stored_coefficient(g^(-2) * a)
        @test Symbolics.simplify(real(negative_integer_power) - 1 / g^2) == 0

        polynomial = stored_coefficient(g * a)
        @test isequal(
            2 / polynomial,
            stored_coefficient((2 / g) * a),
        )
        @test convert(SecondQuantizedAlgebra.Coeff, polynomial) === polynomial
        @test isequal(
            convert(SecondQuantizedAlgebra.Coeff, Complex(Num(g), Num(0))),
            polynomial,
        )
        negative_exponent = parse(Int, "-2")
        @test isequal(
            polynomial^negative_exponent,
            stored_coefficient((1 / g^2) * a),
        )
    end

    @testset "public symbolic conversion boundaries" begin
        @variables θ ω t z::Number
        raw_exp = SymbolicUtils.term(
            exp,
            SymbolicUtils.term(
                *,
                Symbolics.IM,
                SymbolicUtils.unwrap(θ);
                type = Complex{Real}, shape = UnitRange{Int}[],
            );
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_cis = SymbolicUtils.term(
            cis,
            SymbolicUtils.unwrap(θ);
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_phase = SymbolicUtils.term(
            expim,
            SymbolicUtils.unwrap(θ);
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_power = SymbolicUtils.term(
            ^,
            raw_phase,
            2;
            type = Complex{Real}, shape = UnitRange{Int}[],
        )

        @test isequal(exponential_form(Num(raw_exp)), expim(θ))
        @test isequal(exponential_form(Num(raw_cis)), expim(θ))
        @test isequal(exponential_form(Num(raw_power)), expim(2θ))
        raw_complex = SymbolicUtils.term(
            complex,
            SymbolicUtils.unwrap(θ),
            SymbolicUtils.unwrap(ω);
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        @test isequal(
            to_num(exponential_form(Num(raw_complex))),
            Complex(θ, ω),
        )
        @test isequal(
            exponential_form(Complex(Num(raw_phase), Num(0))),
            expim(θ),
        )
        raw_conj_phase = SymbolicUtils.term(
            conj,
            raw_phase;
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_real_phase = SymbolicUtils.term(
            real,
            raw_phase;
            type = Real, shape = UnitRange{Int}[],
        )
        raw_imag_phase = SymbolicUtils.term(
            imag,
            raw_phase;
            type = Real, shape = UnitRange{Int}[],
        )
        raw_abs_phase = SymbolicUtils.term(
            abs,
            raw_phase;
            type = Real, shape = UnitRange{Int}[],
        )
        raw_abs2_phase = SymbolicUtils.term(
            abs2,
            raw_phase;
            type = Real, shape = UnitRange{Int}[],
        )
        raw_unary_sum = SymbolicUtils.term(
            +,
            raw_conj_phase,
            raw_real_phase,
            raw_imag_phase,
            raw_abs_phase,
            raw_abs2_phase,
            SymbolicUtils.unwrap(z);
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        raw_unary_outer = SymbolicUtils.term(
            exp,
            raw_unary_sum;
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        unary = to_num(exponential_form(Num(raw_unary_outer)))
        @test occursin("cos(θ)", string(unary))
        @test occursin("sin(θ)", string(unary))
        @test occursin("exp(-im*θ)", string(unary))

        raw_nonpolynomial = SymbolicUtils.unwrap(exp(1 + z))
        raw_conj_nonpolynomial = SymbolicUtils.term(
            conj,
            raw_nonpolynomial;
            type = Complex{Real}, shape = UnitRange{Int}[],
        )
        conjugate_terms = phase_terms(
            exponential_form(Num(raw_conj_nonpolynomial)),
        )
        @test length(conjugate_terms) == 1
        @test iszero(conjugate_terms[1].phase)
        @test isequal(
            real(conjugate_terms[1].amplitude),
            real(exp(1 + z)),
        )
        @test isequal(
            imag(conjugate_terms[1].amplitude),
            -imag(exp(1 + z)),
        )
        @test isequal(exponential_form(Num(θ)), θ)
        @test isequal(exponential_form(SymbolicUtils.unwrap(θ)), θ)
        @test isequal(trigonometric_form(Num(θ)), θ)
        @test isequal(trigonometric_form(SymbolicUtils.unwrap(θ)), θ)
        @test_throws ArgumentError phase_terms(
            exponential_form(
                Num(
                    SymbolicUtils.term(
                        ^,
                        raw_phase,
                        sqrt(2);
                        type = Complex{Real}, shape = UnitRange{Int}[],
                    )
                )
            )
        )
        @test isequal(exponential_form(exp(ω * t)), exp(ω * t))
    end

    @testset "raw phase substitutions reject non-real angles" begin
        @variables θ
        raw_phase = SymbolicUtils.term(
            expim,
            SymbolicUtils.unwrap(θ);
            type = Complex{Real},
            shape = UnitRange{Int}[],
        )
        raw_product = SymbolicUtils.term(
            *,
            raw_phase,
            SymbolicUtils.unwrap(1 + θ);
            type = Complex{Real},
            shape = UnitRange{Int}[],
        )
        @test_throws ArgumentError substitute(
            Num(raw_product) * a,
            Dict(θ => 1im),
        )
    end

    @testset "public phase projections and expression trees" begin
        @variables θ φ ω t g κ z::Number
        phase = expim(θ)

        @test isequal(exponential_form(exp(im * θ)), phase)
        @test isequal(exponential_form(cis(θ)), phase)
        @test isequal(exponential_form(conj(phase)), expim(-θ))
        @test isequal(real(phase), cos(θ))
        @test isequal(imag(phase), sin(θ))
        @test isequal(real(conj(phase)), cos(θ))
        @test isequal(imag(conj(phase)), -sin(θ))
        @test iszero(
            simplify((conj(cos(z))^2 + conj(sin(z))^2) * a - a),
        )
        @test isequal(real(im * phase), -sin(θ))
        @test isequal(imag(im * phase), cos(θ))
        @test abs(phase) == 1
        @test abs2(phase) == 1
        @test_throws MethodError abs(exponential_form(g))
        @test_throws MethodError abs2(exponential_form(g))

        raw = exponential_form(Complex(Num(z), Num(z)))
        @test isequal(to_num(raw), Complex(Num(z), Num(z)))
        @test isequal(real(raw), z)
        @test isequal(imag(raw), z)
        @test isequal(real(conj(raw)), conj(z))
        @test isequal(imag(conj(raw)), -conj(z))
        product = phase * expim(φ) * exp(ω * t)
        quotient = phase / (expim(φ) * exp(ω * t))
        @test isequal(phase_terms(product)[1].phase, θ + φ)
        @test isequal(phase_terms(quotient)[1].phase, θ - φ)
        @test length(phase_terms(exponential_form(exp(ω * t)))) == 1
        @test length(phase_terms(exponential_form(exp(ω * t)^(1 // 2)))) == 1
        @test_throws ArgumentError phase_terms(1 / (1 + phase))

        @test_throws ArgumentError substitute(phase, Dict(θ => z))
        @test_throws ArgumentError substitute(phase, Dict(θ => 1 + 2im))
        @test_throws ArgumentError substitute(phase, Dict(θ => 1 + 0im))

        differentiated = Symbolics.derivative(phase, θ)
        @test Symbolics.simplify(differentiated - im * phase) == 0
    end
end
