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

@testset "Hermitian formal exponential" begin
    C = im * (lowering - lowering')
    Cm = ComplexF64[0 im; -im 0]
    for be in backends
        for (expr, expected) in (
                (SQA.expim(C), exp(im * Cm)),
                (SQA.expim(cos(A)), exp(im * cos(Am))),
            )
            result = SQA.to_numeric(expr, atom, 2; backend = be)
            matrix = operator_matrix(result)
            @test matrix ≈ expected
            @test matrix' * matrix ≈ ComplexF64[1 0; 0 1]
        end
        projector = Transition(atom, :σ, 2, 2)
        replacement = SQA.to_numeric(lowering, atom, 2; backend = be)
        result = SQA.to_numeric(
            SQA.expim(projector), atom, 2; backend = be,
            operators = Dict(projector => replacement), adjoint_ops = false,
        )
        @test operator_matrix(result) ≈ exp(im * operator_matrix(replacement))
    end
end
