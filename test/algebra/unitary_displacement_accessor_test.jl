using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: Coeff
using Test
using Symbolics: @variables, Num

@testset "displacement accessor" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)
    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)
    @variables α::Number θ ω F Ω t
    A = Coeff(Num(α))

    @testset "static displacement" begin
        U = Displace(a, α)
        @test displacement(U, a) == A
        @test displacement(U, a') == conj(A)
        @test iszero(displacement(Displace(a, 0), a))
    end

    @testset "moving displacement" begin
        U = Displace(a, α * t, t)
        @test displacement(U, a) == A * t
    end

    @testset "pure rotation" begin
        @test iszero(displacement(Rotation(a, θ), a))
        @test iszero(displacement(Rotation(a, ω * t, t), a))
    end

    @testset "untouched site" begin
        b = Destroy(FockSpace(:other), :b)
        @test iszero(displacement(Displace(b, α), a))
    end

    @testset "composition order" begin
        # The left factor acts first: a -> a + α -> exp(-iθ)a + α.
        U = Displace(a, α) * Rotation(a, θ)
        @test displacement(U, a) == A
        # Reversed: a -> exp(-iθ)a -> exp(-iθ)(a + α).
        V = Rotation(a, θ) * Displace(a, α)
        @test iszero(simplify(conjugate(a, V) - expim(-θ) * (a + α)))
        @test displacement(V, a) == expim(-θ) * A
    end

    @testset "displacement frame" begin
        U = DisplacementFrame(a, ω * a' * a + F * a' + F * a)
        @test displacement(U, a) == Coeff(-F / ω)
        Hlin = ω * a' * a + F * expim(-Ω * t) * a' + F * expim(Ω * t) * a
        moving = DisplacementFrame(a, Hlin, t)
        amplitude = displacement(moving, a)
        @test iszero(simplify(amplitude * a + F / (ω - Ω) * expim(-Ω * t) * a))
        @test iszero(simplify(conjugate(a, moving) - a - amplitude))
    end

    @testset "phase space" begin
        U = Displace(x, p, 1, 2)
        @test displacement(U, x) == 1
        @test displacement(U, p) == 2
    end

    @testset "refusal" begin
        spin = SpinSpace(:spin)
        Sx = Spin(spin, :S, 1)
        @test_throws ArgumentError displacement(Rotation(Sx, 3, θ), Sx)
    end
end
