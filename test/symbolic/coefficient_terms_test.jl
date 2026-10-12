using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: CoefficientTerm, coefficient_terms, has_coefficient_terms,
    n_levels, ground_state, transition_levels, operator_axis, to_cnum

rebuild(terms) = sum((t.scalar * prod((x^e for (x, e) in t.powers); init = 1) for t in terms); init = 0)

@testset "coefficient terms" begin
    @variables g κ

    @testset "terms rebuild the coefficient" begin
        for x in (2g^2 * κ + 1 // 3, (1 + 2im) * g - κ^3 / 7, sqrt(2) * g, 1 / g, sin(g) * κ, big(2)^70 * g + 3)
            @test isequal(to_cnum(rebuild(coefficient_terms(x))), to_cnum(x))
        end
    end

    @testset "scalars are exact or floating" begin
        exact = only(coefficient_terms(3 // 4 * g))
        @test exact.scalar isa Complex{Rational{BigInt}}
        @test exact.scalar == 3 // 4
        @test isequal(only(exact.powers), g => 1)
        inexact = only(coefficient_terms(0.25 * g))
        @test inexact.scalar isa ComplexF64
        @test inexact.scalar == 0.25
        big_term = only(t for t in coefficient_terms(big(2)^70 * g + 3) if !isempty(t.powers))
        @test big_term.scalar == big(2)^70
    end

    @testset "numbers and zero" begin
        @test isempty(coefficient_terms(0))
        @test isempty(coefficient_terms(g - g))
        number = only(coefficient_terms(2 - im))
        @test isempty(number.powers)
        @test number.scalar == 2 - im
        @test isequal(only(only(coefficient_terms(1 / g)).powers), g => -1)
    end

    @testset "trees without terms are refused" begin
        for x in (g^g, 1 / (g + κ), (g + κ)^(1 // 3))
            @test !has_coefficient_terms(x)
            @test_throws ArgumentError coefficient_terms(x)
        end
        @test has_coefficient_terms(exp(g))
        @test has_coefficient_terms(2)
    end
end

@testset "operator accessors" begin
    h = NLevelSpace(:atom, 3, 2)
    σ = Transition(h, :σ, 1, 3)
    @test n_levels(σ) == 3
    @test ground_state(σ) == 2
    @test n_levels(Transition(NLevelSpace(:qubit, (:g, :e)), :s, :g, :e)) == 2
    @test_throws ArgumentError n_levels(Destroy(FockSpace(:c), :a))
    @test_throws ArgumentError ground_state(Pauli(PauliSpace(:p), :σ, 1))
    @test transition_levels(σ) == (1, 3)
    @test transition_levels(σ') == (3, 1)
    @test operator_axis(Pauli(PauliSpace(:p), :σ, 2)) == 2
    @test operator_axis(Spin(SpinSpace(:s), :S, 3)) == 3
    @test_throws ArgumentError operator_axis(σ)
end
