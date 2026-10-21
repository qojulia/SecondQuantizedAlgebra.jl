import SecondQuantizedAlgebra as SQA
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
import SciMLOperators
import SparseArrays
using SecondQuantizedAlgebra: NLevelSpace, Transition
using Symbolics: @variables
using Test: @test, @testset, @test_throws

atom = NLevelSpace(:atom, 2)
lowering = Transition(atom, :σ, 1, 2)
A = 0.2 + 0.4 * (lowering + lowering')
B = 0.7 * Transition(atom, :σ, 2, 2)
Am = ComplexF64[0.2 0.4; 0.4 0.2]
Bm = ComplexF64[0 0; 0 0.7]
backends = (SQA.QuantumOpticsBackend(), SQA.QuantumToolboxBackend())
operator_matrix(op::Union{QOB.Operator, QTB.QuantumObject}) = Matrix{ComplexF64}(op.data)

@testset "Operator functions on a finite basis" begin
    @variables θ

    @test all(backends) do be
        result = SQA.to_numeric(cos(A + B), atom, 2; backend = be)
        operator_matrix(result) ≈ cos(Am + Bm)
    end

    @test all(backends) do be
        result = SQA.to_numeric(cos(A) * sin(B) + sin(cos(A + B)), atom, 2; backend = be)
        operator_matrix(result) ≈ cos(Am) * sin(Bm) + sin(cos(Am + Bm))
    end

    @test all(backends) do be
        result = SQA.to_numeric(
            (2 + 3im) * cos(θ * A + B), atom, 2;
            backend = be, parameter = Dict(θ => 0.3)
        )
        operator_matrix(result) ≈ (2 + 3im) * cos(0.3Am + Bm)
    end

    @test_throws ArgumentError SQA.to_numeric(
        cos(θ * A), atom, 2;
        backend = SQA.QuantumOpticsBackend()
    )
end
