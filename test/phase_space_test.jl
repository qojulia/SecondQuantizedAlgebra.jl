using SecondQuantizedAlgebra
using QuantumOpticsBase
using SymbolicUtils
using Test

@testset "phase_space" begin
    @testset "Single Pauli Space" begin
        h_ps1 = PhaseSpace(:motion)
        @qnumbers x::Position(h_ps1, 1)
        p = Momentum(h_ps1, :p, 1)

        # Basic Phase Space Algebra
        @test isequal(p*x, x*p - 1im)
        @test isequal(simplify(p*x*x), simplify(x*x*p - 2im*x))
        @test isequal(substitute(x*x*p, Dict([x, p] .=> [2, 3])), 2*2*3)
        @test isequal(adjoint(x), x)
        @test isequal(adjoint(p), p)
        @test SecondQuantizedAlgebra.ismergeable(p, x)
    end

    @testset "Multi-Pauli Space Operations" begin
        hps1 = PhaseSpace(:motion1)
        hps2 = PhaseSpace(:motion2)
        h = hps1 ⊗ hps2

        x1 = Position(h, :x_1, 1)
        p1 = Momentum(h, :p_1, 1)
        x2 = Position(h, :x_2, 2)
        p2 = Momentum(h, :p_2, 2)

        @test isequal(p1*x2, x2*p1)
        @test isequal(simplify(p1*x2*x1), simplify((x1*p1-1im)*x2))
        @test isequal(adjoint(x1), x1)
        @test isequal(adjoint(p1), p1)
        @test SecondQuantizedAlgebra.ismergeable(p1, x1)
        @test !SecondQuantizedAlgebra.ismergeable(p1, x2)
    end

    @testset "Pauli Hamiltonian Construction and Equations" begin
        hps1 = PhaseSpace(:motion1)
        hps2 = PhaseSpace(:motion2)
        h = hps1 ⊗ hps2

        x1 = Position(h, :x_1, 1)
        p1 = Momentum(h, :p_1, 1)
        x2 = Position(h, :x_2, 2)
        p2 = Momentum(h, :p_2, 2)

        @cnumbers ω m miv c
        H = 0.5*miv*p1^2 + 0.5m*ω^2*x1^2 + c*x1^4

        @test commutator(H, x2) == 0
        @test commutator(H, p2) == 0
        @test isequal(simplify(1im*commutator(H, x1)), simplify(p1*miv))
        @test isequal(simplify(1im*commutator(H, p1)), simplify(-x1*m*ω^2 - 4*c*x1^3))

        # Test that Hamiltonian construction works
        @test H isa SecondQuantizedAlgebra.QAdd
    end
end
