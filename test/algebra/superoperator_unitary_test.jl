using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, hamiltonian, channels,
    exponential_form

# Constant part of the Heisenberg drift of `a`, i.e. the residual coherent drive.
function drift_offset(L, a)
    d = transpose(L)(a)
    for (term, c) in d
        isempty(term.ops) && return c
    end
    return nothing
end

# Squeezing carries float coefficients, so compare maps up to roundoff.
function approx_equal(A, B; atol = 1.0e-12)
    for (_, c) in A - B
        abs(SecondQuantizedAlgebra.to_complex(c)) <= atol || return false
    end
    return true
end

@testset "Superoperator frames" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    @variables t::Real ω::Real ωd::Real F::Real Δ::Real χ::Real

    @testset "frame invariant" begin
        H = Δ * a' * a + F * (a + a') + χ * a' * a' * a * a
        for U in (Displace(a, 2 + 1im), Squeeze(a, 0.3), RotatingFrame(Δ * a' * a, t))
            @test transform(hamiltonian_action(H), U) == hamiltonian_action(transform(H, U))
        end
    end

    @testset "conjugation is a homomorphism" begin
        A = sprepost(a, a') + spre(a' * a)
        B = spost(a * a) - 2 * sprepost(a', a)
        for U in (Displace(a, 1 // 2), RotatingFrame(ω * a' * a, t))
            @test conjugate(A ∘ B, U) == conjugate(A, U) ∘ conjugate(B, U)
        end
        S = Squeeze(a, 0.2)
        @test approx_equal(conjugate(A ∘ B, S), conjugate(A, S) ∘ conjugate(B, S))
        L = hamiltonian_action(Δ * a' * a) + dissipator(a)
        U = Displace(a, 3)
        @test SuperOperator(conjugate(L, U)) == conjugate(SuperOperator(L), U)
    end

    @testset "rotating frame of a damped driven cavity" begin
        H = ω * a' * a + F * (a' * exp(-im * ωd * t) + a * exp(im * ωd * t))
        L = hamiltonian_action(H) + dissipator(a)
        L2 = exponential_form(transform(L, RotatingFrame(ωd * a' * a, t)))
        @test L2 == hamiltonian_action((ω - ωd) * a' * a + F * (a + a')) + dissipator(a)
    end

    @testset "dissipative displacement frame" begin
        # Pure loss: the drift is λa + c with λ = -iΔ - κ/2 and c = -iF.
        κ = 1 // 2
        L = hamiltonian_action(2 * a' * a + 3 * (a + a')) + κ * dissipator(a)
        U = DisplacementFrame(a, L)
        @test drift_offset(transform(L, U), a) === nothing

        # Thermal bath: only κ₋ - κ₊ damps the mean field.
        κm, κp = 1 // 2, 1 // 10
        Lth = hamiltonian_action(2 * a' * a + 3 * (a + a')) + κm * dissipator(a) + κp * dissipator(a')
        λ = -2im - (κm - κp) // 2
        @test transpose(Lth)(a) == λ * a - 3im
        @test drift_offset(transform(Lth, DisplacementFrame(a, Lth)), a) === nothing

        # Squeezed bath 𝒟[μa + νa†] damps at (|μ|² - |ν|²)/2.
        μ, ν = 2, 1
        Lsq = hamiltonian_action(a' * a + (a + a')) + dissipator(μ * a + ν * a')
        @test transpose(Lsq)(a) == (-1im - (μ^2 - ν^2) // 2) * a - 1im
        @test drift_offset(transform(Lsq, DisplacementFrame(a, Lsq)), a) === nothing

        # Time-dependent drive.
        Ld = hamiltonian_action(Δ * a' * a + F * (a * exp(im * ωd * t) + a' * exp(-im * ωd * t))) +
            (1 // 2) * dissipator(a)
        Ud = DisplacementFrame(a, Ld, t)
        @test Ud isa UnitaryTransform
        Lt = exponential_form(simplify(transform(Ld, Ud)))
        @test drift_offset(Lt, a) === nothing
    end

    @testset "non-affine drift throws" begin
        L = hamiltonian_action(χ * a' * a' * a * a) + dissipator(a)
        @test_throws ArgumentError DisplacementFrame(a, L)
        @test_throws ArgumentError DisplacementFrame(a, dissipator(a * a))
        # Dephasing is affine: transpose(𝒟[a†a])(a) = -a/2.
        @test transpose(dissipator(a' * a))(a) == -(1 // 2) * a
    end

    @testset "squeezed-frame Kerr oscillator" begin
        L = hamiltonian_action(Δ * a' * a + χ * a' * a' * a * a) + dissipator(a)
        U = Squeeze(a, 0.25)
        @test approx_equal(SuperOperator(transform(L, U)), transform(SuperOperator(L), U))
        @test approx_equal(transpose(transform(L, U))(1), zero(transpose(L)(1)))
    end
end
