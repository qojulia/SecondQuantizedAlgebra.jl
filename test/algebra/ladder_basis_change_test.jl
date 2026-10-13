using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: quadratures, ladder_operator, FrequencyChange

@testset "Ladder bases of one oscillator" begin
    a = Destroy(FockSpace(:a), :a)
    @variables ω0::Real ωr::Real
    oscillator(x, p, ω) = p * p / 2 + ω^2 * x * x / 2
    # H0 = p²/2 + ω0² x²/2 in the ladder basis at ωr, with the mode `a` standing for b.
    reference_form(ω0, ωr) =
        (ω0^2 + ωr^2) / (2ωr) * (a' * a + 1 // 2) + (ω0^2 - ωr^2) / (4ωr) * (a * a + a' * a')

    @testset "quadratures are canonical" begin
        for ω in (1, 2, 1 // 8, 3, ω0)
            x, p = quadratures(a, ω)
            @test iszero(simplify(commutator(x, p) - im))
            @test iszero(simplify(x - adjoint(x)))
            @test iszero(simplify(p - adjoint(p)))
        end
        @test isequal(quadratures(a'), quadratures(a))
    end

    @testset "exact coefficients for rational frequencies" begin
        x, p = quadratures(a, 1 // 8)
        @test isequal(x, 2 * a + 2 * a')
        @test isequal(p, (im // 4) * (a' - a))
        x, p = quadratures(a, 2)
        @test isequal(x, (a + a') // 2)
        @test isequal(p, im * (a' - a))
    end

    @testset "own basis diagonalizes the oscillator" begin
        for ω in (2, 1 // 8, ω0)
            H = simplify(oscillator(quadratures(a, ω)..., ω))
            @test iszero(simplify(H - ω * (a' * a + 1 // 2)))
        end
    end

    @testset "reference basis" begin
        H = oscillator(quadratures(a, ωr)..., ω0)
        @test iszero(simplify(H - reference_form(ω0, ωr)))
        H = simplify(oscillator(quadratures(a, 1 // 2)..., 2))
        @test isequal(H, reference_form(2, 1 // 2))
        @test isequal(H, 17 // 8 + 17 // 4 * a' * a + 15 // 8 * (a * a + a' * a'))
    end

    @testset "ladder_operator inverts quadratures" begin
        for ω in (1, 2, 1 // 8, ω0)
            @test iszero(simplify(ladder_operator(quadratures(a, ω)..., ω) - a))
        end
        @test isequal(ladder_operator(quadratures(a, 2)..., 1 // 2), 5 // 4 * a - 3 // 4 * a')

        b = ladder_operator(quadratures(a, ω0)..., ωr)
        @test iszero(simplify(commutator(b, adjoint(b)) - 1))
        root = sqrt(2ω0) * sqrt(2ωr)
        @test iszero(simplify(b - ((ω0 + ωr) / root * a - (ω0 - ωr) / root * a')))
    end

    @testset "ladder_operator on Position and Momentum" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        c = ladder_operator(x, p, 2)
        @test isequal(c, x + (im // 2) * p)
        @test iszero(simplify(commutator(c, adjoint(c)) - 1))
    end

    @testset "FrequencyChange" begin
        U = FrequencyChange(a, 2, 1 // 2)
        @test U isa UnitaryTransform
        @test isequal(conjugate(a, U), 5 // 4 * a + 3 // 4 * a')
        H0 = simplify(oscillator(quadratures(a, 2)..., 2))
        @test isequal(simplify(conjugate(H0, U)), reference_form(2, 1 // 2))
        @test iszero(simplify(conjugate(conjugate(a, U), inv(U)) - a))

        U = FrequencyChange(a', ω0, ωr)
        H0 = simplify(oscillator(quadratures(a, ω0)..., ω0))
        @test iszero(simplify(conjugate(H0, U) - reference_form(ω0, ωr)))
        b = conjugate(a, U)
        @test iszero(simplify(commutator(b, adjoint(b)) - 1))
        @test iszero(simplify(conjugate(b, inv(U)) - a))
        @test iszero(simplify(conjugate(a, FrequencyChange(a, ω0, ω0)) - a))
        # The ω0 ladder operator written in the ωr basis, as built from the quadratures.
        a_in_b = ladder_operator(quadratures(a, ωr)..., ω0)
        @test iszero(simplify(conjugate(a, U) - a_in_b))
    end

    @testset "FrequencyChange is Squeeze(a, log(ω0/ωr)/2)" begin
        U = FrequencyChange(a, 2, 1 // 2)
        S = Squeeze(a, log(4.0) / 2)
        lhs = conjugate(a, U)
        rhs = conjugate(a, S)
        coefficients(q) = Dict(
            term.ops => ComplexF64(SecondQuantizedAlgebra.native_scalar(c)) for (term, c) in q
        )
        @test !isempty(coefficients(lhs))
        @test keys(coefficients(lhs)) == keys(coefficients(rhs))
        for (ops, c) in coefficients(lhs)
            @test c ≈ coefficients(rhs)[ops]
        end
    end

    @testset "refusals" begin
        @test_throws ArgumentError quadratures(a, 0)
        @test_throws ArgumentError quadratures(a, -1)
        @test_throws ArgumentError ladder_operator(quadratures(a)..., -2)
        @test_throws ArgumentError FrequencyChange(a, 1, 0)
        @test_throws ArgumentError quadratures(Position(PhaseSpace(:q), :x), 1)
        @test_throws ArgumentError FrequencyChange(Pauli(PauliSpace(:s), :σ, 1), 1, 2)
    end
end
