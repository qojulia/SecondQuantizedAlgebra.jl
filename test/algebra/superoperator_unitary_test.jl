using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, hamiltonian, channels,
    exponential_form, expim, simplify, to_complex

# Constant part of the Heisenberg drift of `a`, i.e. the residual coherent drive.
function drift_offset(L, a)
    d = transpose(L)(a)
    for (term, c) in d
        isempty(term.ops) && return c
    end
    return nothing
end

scalar_coefficient(q) = (
    for (term, c) in q
        isempty(term.ops) && return c
    end; 0
)
operator_coefficient(q, o) = (
    for (term, c) in q
        term.ops == [o] && return c
    end; 0
)
# Coefficients of the sandwiches with a single operator, i.e. the residual drive of a map.
linear_coefficients(S) =
    [c for (term, c) in S.arguments if length(term.left) + length(term.right) == 1]
drive_free(S) = all(iszero ∘ simplify, linear_coefficients(S))

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

    @testset "displacement frame of a general linear map" begin
        # Caldeira-Leggett (quantum Brownian motion) generator: not GKSL, and its damping acts
        # on p only, so the drift of `b` mixes in `b'`.
        @variables q0::Real γ::Real D::Real
        hq = FockSpace(:qbm)
        @qnumbers b::Destroy(hq)
        ω0 = 2q0^2
        x = (b + b') / (2q0)
        p = im * q0 * (b' - b)
        ad(A) = spre(A) - spost(A)
        anti(A) = spre(A) + spost(A)
        qbm(H) = SuperOperator(hamiltonian_action(H)) - D / 2 * ad(x) ∘ ad(x) -
            im * γ / 2 * ad(x) ∘ anti(p)
        drive = (expim(ωd * t) + expim(-ωd * t)) / 2
        L = qbm(p^2 / 2 + ω0^2 * x^2 / 2 - F * x * drive)
        @test transpose(L)(b) == (-γ / 2 - im * ω0) * b + γ / 2 * b' + im * F * drive / (2q0)

        # The displacement is the classical driven damped oscillator, x_c = Re[X e^{-iωd t}].
        U = DisplacementFrame(b, L, t)
        X = F / (ω0^2 - ωd^2 - im * γ * ωd)
        e = expim(-ωd * t)
        xc = (X * e + conj(X) * conj(e)) / 2
        pc = -im * ωd * (X * e - conj(X) * conj(e)) / 2
        α = scalar_coefficient(conjugate(b, U) - b)
        @test iszero(simplify((α + conj(α)) / (2q0) - xc))
        @test iszero(simplify(im * q0 * (conj(α) - α) - pc))
        @test drive_free(transform(L, U))

        # The gauge is the exact moving-frame gauge of `Displace(b, α, t)`.
        αdot = q0 * pc - im * ωd^2 * (X * e + conj(X) * conj(e)) / (4q0)
        @test iszero(simplify(operator_coefficient(U.gauge, b') + im * αdot))
        values = Dict(q0 => 3 // 4, ωd => 7 // 5, γ => 1 // 5, D => 2 // 5, F => 9 // 10, t => 3 // 5)
        reference = Displace(b, q0 * xc + im * pc / (2q0), t)
        @test abs(to_complex(scalar_coefficient(substitute(U.gauge - reference.gauge, values)))) <
            1.0e-12

        # Static force: x_c = F/ω0², p_c = 0.
        Ls = qbm(p^2 / 2 + ω0^2 * x^2 / 2 - F * x)
        Us = DisplacementFrame(b, Ls)
        @test iszero(simplify(scalar_coefficient(conjugate(b, Us) - b) - q0 * F / ω0^2))
        @test drive_free(transform(Ls, Us))

        # A squeezed Lindbladian, whose drift mixes `b` and `b'`, now also has a frame.
        Lsq = hamiltonian_action(Δ * b' * b + χ / 2 * (b^2 + b'^2) + F * (b + b')) +
            (1 // 2) * dissipator(b)
        @test drive_free(SuperOperator(transform(Lsq, DisplacementFrame(b, Lsq))))
        Lsqt = hamiltonian_action(
            Δ * b' * b + χ / 2 * (b^2 + b'^2) + F * (b * expim(ωd * t) + b' * expim(-ωd * t)),
        ) + (1 // 2) * dissipator(b)
        @test drive_free(SuperOperator(transform(Lsqt, DisplacementFrame(b, Lsqt, t))))
    end

    @testset "quadrature displacement frame of a linear map" begin
        @variables ω::Real γ::Real D::Real
        hq = PhaseSpace(:q)
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        ad(A) = spre(A) - spost(A)
        anti(A) = spre(A) + spost(A)
        qbm(H) = SuperOperator(hamiltonian_action(H)) - D / 2 * ad(x) ∘ ad(x) -
            im * γ / 2 * ad(x) ∘ anti(p)
        L = qbm(p^2 / 2 + ω^2 * x^2 / 2 - F * x * (expim(ωd * t) + expim(-ωd * t)) / 2)
        U = DisplacementFrame(x, p, L, t)
        X = F / (ω^2 - ωd^2 - im * γ * ωd)
        e = expim(-ωd * t)
        xc = (X * e + conj(X) * conj(e)) / 2
        pc = -im * ωd * (X * e - conj(X) * conj(e)) / 2
        @test iszero(simplify(scalar_coefficient(conjugate(x, U) - x) - xc))
        @test iszero(simplify(scalar_coefficient(conjugate(p, U) - p) - pc))
        @test drive_free(transform(L, U))

        Ls = qbm(p^2 / 2 + ω^2 * x^2 / 2 - F * x)
        Us = DisplacementFrame(x, p, Ls)
        @test conjugate(x, Us) == x + F / ω^2
        @test conjugate(p, Us) == 1 * p
        @test drive_free(transform(Ls, Us))
        @test_throws ArgumentError DisplacementFrame(x, p, Ls + SuperOperator(hamiltonian_action(x^4)))
    end

    @testset "linear-map displacement frame refusals" begin
        for H in (a' * a + a', a' * a + a)
            S = SuperOperator(hamiltonian_action(H))
            @test_throws ArgumentError DisplacementFrame(a, S)
            @test_throws ArgumentError DisplacementFrame(a, S, t)
        end
        # Normal-mode frequency sqrt(5² - 3²) = 4: resonant at ν = ±4, singular when Δ = χ.
        resonant = 5 * a' * a + (3 // 2) * (a^2 + a'^2) + F * (a + a') * (expim(4t) + expim(-4t))
        @test_throws ArgumentError DisplacementFrame(a, SuperOperator(hamiltonian_action(resonant)), t)
        singular = 3 * a' * a + (3 // 2) * (a^2 + a'^2) + F * (a + a')
        @test_throws ArgumentError DisplacementFrame(a, SuperOperator(hamiltonian_action(singular)))
        # A non-Hermitian squeezing generator does not preserve Hermiticity.
        S = -im * (spre(a'^2) - spost(a'^2)) + SuperOperator(hamiltonian_action(a' * a + F * (a + a')))
        @test_throws ArgumentError DisplacementFrame(a, S)
        # A time-dependent drift matrix has no harmonic bounded solution.
        modulated = hamiltonian_action(a' * a + χ * t * (a^2 + a'^2) + F * (a + a'))
        @test_throws ArgumentError DisplacementFrame(a, modulated, t)
    end
end
