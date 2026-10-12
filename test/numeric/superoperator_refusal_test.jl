# Deliberate numeric refusals for maps.

using SecondQuantizedAlgebra
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperExpr, liouvillian, to_numeric_master

@testset "Superoperator numeric refusals" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    hφ = PhaseSpace(:t)
    φ = Position(hφ, :φ)
    @variables t κ(t) Δ
    b = QOB.FockBasis(3)
    L = hamiltonian_action(a' * a) + dissipator(a)

    @testset "formal maps" begin
        Lf = hamiltonian_action(cos(φ))
        @test_throws ArgumentError to_numeric(Lf, b)
        @test_throws ArgumentError to_numeric(SuperExpr(Lf), b)
        @test_throws ArgumentError to_numeric_master(Lf, hφ, 3; backend = QuantumOpticsBackend())
    end

    @testset "QuantumOptics time dependence" begin
        Lt = κ * dissipator(a)
        tp = Dict(κ => (t -> 1.0))
        @test_throws ArgumentError to_numeric(Lt, b; time_parameter = tp)
        @test_throws ArgumentError to_numeric_master(Lt, h, 3; backend = QuantumOpticsBackend(), time_parameter = tp)
    end

    @testset "averages of maps" begin
        ψ = QOB.fockstate(b, 1)
        @test_throws ArgumentError SecondQuantizedAlgebra.expect(L, ψ)
        @test_throws ArgumentError numeric_average(L, ψ)
        @test_throws ArgumentError average(L)
        @test_throws ArgumentError normal_to_symmetric(L)
        @test_throws ArgumentError symmetric_to_normal(SecondQuantizedAlgebra.SuperOperator(L))
    end

    @testset "non-physical solver inputs" begin
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        Lneg = liouvillian(0, [a1, b1]; rates = [1.0 2.0; 2.0 1.0])
        @test_throws ArgumentError to_numeric_master(Lneg, hp, (2, 2); backend = QuantumOpticsBackend())
        @test_throws ArgumentError to_numeric_master(-1 * dissipator(a), h, 3; backend = QuantumOpticsBackend())
        @test_throws ArgumentError to_numeric_master(hamiltonian_action(im * a' * a), h, 3; backend = QuantumOpticsBackend())
        @test_throws ArgumentError to_numeric_master(Δ * dissipator(a), h, 3; backend = QuantumOpticsBackend())
    end

    @testset "time-dependent Hamiltonian Hermiticity" begin
        @variables Ω::Number η(t)
        QTBB = QuantumToolboxBackend()
        tp = Dict(Ω => (t -> cos(t) + im * sin(t)), η => (t -> 1 + t))
        # A conjugate-paired drive and a real-valued coefficient are Hermitian.
        Lok = hamiltonian_action(Ω * a' + conj(Ω) * a + η * a' * a) + dissipator(a)
        @test to_numeric_master(Lok, h, 3; backend = QTBB, time_parameter = tp).collapse_operators isa Vector
        # The same drive without its conjugate partner is not.
        Lbad = hamiltonian_action(Ω * a' + Ω * a) + dissipator(a)
        @test_throws ArgumentError to_numeric_master(Lbad, h, 3; backend = QTBB, time_parameter = tp)
        Lim = hamiltonian_action(im * η * a' * a) + dissipator(a)
        @test_throws ArgumentError to_numeric_master(Lim, h, 3; backend = QTBB, time_parameter = tp)
        # Solver-parameter coefficients are checked symbolically, with real parameters.
        tpp = Dict(η => ((p, t) -> p[1] * t))
        @test to_numeric_master(hamiltonian_action(η * a' * a), h, 3; backend = QTBB, time_parameter = tpp).hamiltonian isa QTB.QuantumObjectEvolution
        @test_throws ArgumentError to_numeric_master(
            hamiltonian_action(im * η * a' * a), h, 3; backend = QTBB, time_parameter = tpp,
        )
    end

    @testset "time-dependent correlated block" begin
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        Lc = liouvillian(0, [a1, b1]; rates = [1 κ; κ 1])
        tp = Dict(κ => (t -> 0.1))
        @test_throws ArgumentError to_numeric_master(
            Lc, hp, (2, 2); backend = QuantumToolboxBackend(), time_parameter = tp,
        )
    end
end
