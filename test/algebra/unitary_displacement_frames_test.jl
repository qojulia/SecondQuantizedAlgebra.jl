using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim

@testset "Hamiltonian-derived displacement frames" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    @variables ω Ω η g t dx dp ωd K
    @variables envelope(t)

    @testset "static Fock displacement frame" begin
        reference = ω * a' * a + η * (a + a') + g
        U = DisplacementFrame(a, reference)

        @test U isa UnitaryTransform
        @test iszero(simplify(conjugate(a, U) - (a - η / ω)))
        @test iszero(
            simplify(transform(reference, U) - (ω * a' * a - η^2 / ω + g)),
        )
        @test iszero(simplify(conjugate(a, inv(U)) - (a + η / ω)))

        undriven = DisplacementFrame(a, ω * a' * a)
        @test iszero(simplify(conjugate(a, undriven) - a))
    end

    @testset "bounded harmonic Fock displacement frame" begin
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')
        U = DisplacementFrame(a, reference, t)

        expected =
            im * Ω / (2 * (ωd - ω)) * expim(-ωd * t) -
            im * Ω / (2 * (ωd + ω)) * expim(ωd * t)
        @test iszero(simplify(conjugate(a, U) - a - expected))

        transformed = transform(reference, U)
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops == Op[a', a]
        end

        constant = DisplacementFrame(a, ω * a' * a + η * (a + a'), t)
        @test iszero(simplify(conjugate(a, constant) - (a - η / ω)))

        numeric_constant = DisplacementFrame(a, ω * a' * a + a + a', t)
        @test iszero(simplify(conjugate(a, numeric_constant) - (a - 1 / ω)))

        multitone_drive = η * cos(ωd * t) + g * sin(2ωd * t)
        multitone_reference = ω * a' * a + multitone_drive * (a + a')
        multitone = DisplacementFrame(a, multitone_reference, t)
        for (term, _) in transform(multitone_reference, multitone)
            isempty(term.ops) && continue
            @test term.ops == Op[a', a]
        end

        Kerr = (K / 2) * a'^2 * a^2
        @test iszero(
            simplify(
                transform(reference + Kerr, U) -
                    transform(reference, U) - conjugate(Kerr, U),
            ),
        )
    end

    @testset "Fock frame refusal is exact" begin
        other = Destroy(FockSpace(:other), :b)
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')

        @test_throws ArgumentError DisplacementFrame(a, reference + K * a'^2 * a^2, t)
        @test_throws ArgumentError DisplacementFrame(a, reference + g * other, t)
        @test_throws ArgumentError DisplacementFrame(a, ω * a' * a + η * a', t)
        @test_throws ArgumentError DisplacementFrame(a, im * a' * a + η * (a + a'))
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + η * (a + a') + im,
        )
        @test_throws ArgumentError DisplacementFrame(
            a, (ω + g * cos(ωd * t)) * a' * a + η * (a + a'), t,
        )
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + envelope * (a + a'), t,
        )

        indexed_fock = FockSpace(:indexed_reference)
        i = Index(indexed_fock, :i, 3, indexed_fock)
        indexed_a = IndexedOperator(Destroy(indexed_fock, :c), i)
        @test_throws ArgumentError DisplacementFrame(
            a, Σ(indexed_a' * indexed_a, i),
        )

        nonlinear_phase = expim(t^2)
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + nonlinear_phase * a' + conj(nonlinear_phase) * a, t,
        )

        resonant = expim(-ω * t)
        @test_throws ArgumentError DisplacementFrame(
            a, ω * a' * a + resonant * a' + conj(resonant) * a, t,
        )
        @test_throws ArgumentError DisplacementFrame(a, η * (a + a'), t)
        @test_throws ArgumentError DisplacementFrame(a, η * (a + a'))
    end

    @testset "static quadrature displacement frame" begin
        reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 + η * x + dx * p
        U = DisplacementFrame(x, p, reference)
        determinant = ω * Ω - g^2

        @test U isa UnitaryTransform
        @test iszero(
            simplify(conjugate(x, U) - x - (g * dx - Ω * η) / determinant),
        )
        @test iszero(
            simplify(conjugate(p, U) - p - (g * η - ω * dx) / determinant),
        )

        transformed = simplify(transform(reference, U))
        allowed = Set((Op[x, x], Op[x, p], Op[p, p]))
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops in allowed
        end
    end

    @testset "bounded harmonic quadrature displacement frame" begin
        reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        U = DisplacementFrame(x, p, reference, t)
        @test U isa UnitaryTransform

        transformed = transform(reference, U)
        allowed = Set((Op[x, x], Op[p, p]))
        for (term, _) in transformed
            isempty(term.ops) && continue
            @test term.ops in allowed
        end

        numeric_constant = DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + x + p, t,
        )
        @test numeric_constant isa UnitaryTransform

        multitone_reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 +
            (η * cos(ωd * t) + dx * sin(2ωd * t)) * x +
            dp * cos(3ωd * t) * p
        multitone = DisplacementFrame(x, p, multitone_reference, t)
        transformed_multitone = transform(multitone_reference, multitone)
        allowed_multitone = Set((Op[x, x], Op[x, p], Op[p, p]))
        for (term, _) in transformed_multitone
            isempty(term.ops) && continue
            @test term.ops in allowed_multitone
        end
    end

    @testset "quadrature frame refusal is exact" begin
        reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        other_phase = PhaseSpace(:selected) ⊗ PhaseSpace(:other)
        other_x = Position(other_phase, :other_x, 2)
        nonlinear_phase = expim(t^2)

        @test_throws ArgumentError DisplacementFrame(p, x, reference, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + other_x, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + K * x^3, t)
        @test_throws ArgumentError DisplacementFrame(x, p, reference + im * x, t)
        @test_throws ArgumentError DisplacementFrame(
            x, p,
            ((ω + g * cos(ωd * t)) / 2) * x^2 + (Ω / 2) * p^2 + η * x,
            t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + envelope * x, t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p,
            (ω / 2) * (x^2 + p^2) +
                (nonlinear_phase + conj(nonlinear_phase)) * x,
            t,
        )
        @test_throws ArgumentError DisplacementFrame(
            x, p, (ω / 2) * (x^2 + p^2) + η * cos(ω * t) * x, t,
        )
        @test_throws ArgumentError DisplacementFrame(x, p, (ω / 2) * x^2 + η * p)

        nonunique = DisplacementFrame(x, p, (ω / 2) * x^2)
        @test iszero(simplify(conjugate(x + p, nonunique) - (x + p)))
    end
end
