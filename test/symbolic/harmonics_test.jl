using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: QAdd, Coeff, SuperOperator, spre, spost, sprepost, expim,
    exponential_form, harmonic, harmonics, secular

@testset "Harmonic projection" begin
    @variables ω0::Real ωd::Real ωp::Real γ::Real F::Real θ::Real t::Real
    h = FockSpace(:harmonics)
    a = Destroy(h, :a)
    x = a + a'
    H = ω0 * a' * a - F * x * (expim(ωd * t) + expim(-ωd * t)) / 2
    U = Rotation(a, ωd * t, t)
    reconstruct(parts) = sum(parts[n] * expim(n * ωd * t) for n in keys(parts))

    @testset "rotating-frame Hamiltonian" begin
        Ht = transform(H, U)
        parts = @inferred harmonics(Ht, ωd, t)
        @test sort!(collect(keys(parts))) == [-2, 0, 2]
        @test isequal(parts[0], (ω0 - ωd) * a' * a - F / 2 * x)
        @test isequal(parts[2], -F / 2 * a')
        @test isequal(parts[-2], -F / 2 * a)
        @test iszero(simplify(reconstruct(parts) - exponential_form(Ht)))
        @test isequal(@inferred(secular(Ht, ωd, t)), parts[0])
        @test isequal(@inferred(harmonic(Ht, ωd, t, 2)), parts[2])
        @test iszero(harmonic(Ht, ωd, t, 1))
        @test isequal(secular(H, ωd, t), ω0 * a' * a)
    end

    @testset "scalars and operators" begin
        @test @inferred(harmonic(3, ωd, t, 0)) isa Coeff
        @test isequal(harmonic(3, ωd, t, 0), 3)
        @test iszero(harmonic(3, ωd, t, 1))
        @test isequal(secular(a, ωd, t), 1 * a)
        @test isequal(harmonic(F * expim(-2ωd * t), ωd, t, -2), F)
        @test isequal(harmonics(cos(2ωd * t) * a, ωd, t), Dict(2 => a / 2, -2 => a / 2))
        @test isempty(harmonics(0 * a, ωd, t))
    end

    @testset "offsets and composite base frequencies" begin
        @test isequal(harmonics(expim(ωd * t + θ) * a, ωd, t), Dict(1 => expim(θ) * a))
        @test isequal(harmonics(expim((ωd - ω0) * t) * a, ωd - ω0, t), Dict(1 => 1 * a))
        @test isequal(harmonics(expim(2 * t) * a, 1, t), Dict(2 => 1 * a))
        @test isequal(harmonics(expim(2.0 * t) * a, 1.0, t), Dict(2 => 1 * a))
        @test isequal(harmonics(expim(t) * a, 1 // 2, t), Dict(2 => 1 * a))
    end

    @testset "refusals" begin
        # A second frequency, a fractional or nonlinear phase and a non-phase envelope
        # cannot be assigned to a harmonic of ωd without an assumption.
        @test_throws ArgumentError harmonics(expim(ωp * t) * a, ωd, t)
        @test_throws ArgumentError harmonic(expim(ωp * t) * a + a', ωd, t, 0)
        @test_throws ArgumentError harmonics(expim(ωd * t / 2) * a, ωd, t)
        @test_throws ArgumentError harmonics(expim(ωd * t^2) * a, ωd, t)
        @test_throws ArgumentError harmonics(t * a, ωd, t)
        @test_throws ArgumentError harmonics(a, 0, t)
        @test_throws ArgumentError harmonics(a, ωd * t, t)
        @test_throws ArgumentError harmonics(a, ωd, 2t)
        @test_throws ArgumentError harmonics(expim(t) * a, 2, t)
    end

    @testset "indexed sums keep their scope" begin
        hc = FockSpace(:harmonics_sites)
        i = Index(hc, :i, 3, hc)
        b = IndexedOperator(Destroy(hc, :b), i)
        q = Σ(expim(ωd * t) * b + b' * b, i)
        parts = harmonics(q, ωd, t)
        @test isequal(parts[1], Σ(b, i))
        @test isequal(parts[0], Σ(b' * b, i))
    end

    @testset "superoperator" begin
        D = sprepost(x, x) - spre(x * x) / 2 - spost(x * x) / 2
        S = -im * spre(H) + im * spost(H) + γ * D
        St = transform(S, U)
        parts = @inferred harmonics(St, ωd, t)
        @test sort!(collect(keys(parts))) == [-2, 0, 2]
        @test iszero(simplify(reconstruct(parts) - exponential_form(St)))
        @test isequal(@inferred(harmonic(St, ωd, t, 2)), parts[2])
        Hs = (ω0 - ωd) * a' * a - F / 2 * x
        rwa = -im * spre(Hs) + im * spost(Hs) + γ * (
            sprepost(a, a') - spre(a' * a) / 2 - spost(a' * a) / 2 +
                sprepost(a', a) - spre(a * a') / 2 - spost(a * a') / 2
        )
        @test iszero(simplify(secular(St, ωd, t) - rwa))
    end

    @testset "Lindbladian" begin
        L = hamiltonian_action(H) + γ * dissipator(x)
        Lt = transform(L, U)
        parts = @inferred harmonics(Lt, ωd, t)
        @test sort!(collect(keys(parts))) == [-2, 0, 2]
        lowered = SuperOperator(Lt)
        for n in -3:3
            Ln = harmonic(Lt, ωd, t, n)
            @test Ln isa Lindbladian{QAdd}
            @test isequal(SuperOperator(Ln), harmonic(lowered, ωd, t, n))
        end
        Hs = (ω0 - ωd) * a' * a - F / 2 * x
        @test secular(Lt, ωd, t) ==
            hamiltonian_action(Hs) + γ * dissipator(a) + γ * dissipator(a')

        # Time-dependent jump with a phase-carrying partner: the partner enters conjugated.
        J = a + F * expim(ωd * t) * a'
        Lj = dissipator(J) + γ * expim(-ωd * t) * dissipator(a, a')
        lowered_j = SuperOperator(Lj)
        for n in -2:2
            @test isequal(SuperOperator(harmonic(Lj, ωd, t, n)), harmonic(lowered_j, ωd, t, n))
        end
    end
end
