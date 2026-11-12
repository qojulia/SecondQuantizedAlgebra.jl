import SecondQuantizedAlgebra as SQA
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
using Test: @test, @testset

@testset "Formal positional operator dictionaries" begin
    h = SQA.FockSpace(:formal_positional)
    a = SQA.Destroy(h, :a)
    expr = cos(a + a')
    b = QOB.FockBasis(5)
    qob_a = QOB.destroy(b)
    qtb_a = QTB.destroy(6)
    qob_state = QOB.fockstate(b, 0)
    qtb_state = QTB.fock(6, 0)

    for (target, operator, state) in (
            (b, qob_a, qob_state),
            (6, qtb_a, qtb_state), ([6], qtb_a, qtb_state), ((6,), qtb_a, qtb_state),
        )
        for replacements in (
                Dict{SQA.QSym, Union{}}(), Dict(a => 2 * operator),
                Dict(a => 2 * operator, a' => 2 * operator'),
            )
            actual = SQA.to_numeric(expr, target, replacements)
            expected = SQA.to_numeric(expr, target; operators = replacements, adjoint_ops = false)
            @test Matrix(actual.data) ≈ Matrix(expected.data)
            @test Matrix(SQA.to_numeric(expr, state, replacements).data) ≈ Matrix(actual.data)
        end
        replacements = Dict(a => 2 * operator)
        @test Matrix(SQA.to_numeric(expr, target, replacements).data) ≈
            cos(Matrix((2 * operator + operator').data))
        @test Matrix(SQA.to_numeric(expr, state).data) ≈ Matrix(SQA.to_numeric(expr, target).data)
    end
end
