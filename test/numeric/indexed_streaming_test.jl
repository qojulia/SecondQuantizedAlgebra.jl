import SecondQuantizedAlgebra as SQA
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
using Test: @test, @testset, @test_throws
using Symbolics: @variables

@testset "Streaming indexed sums" begin
    h = SQA.NLevelSpace(:stream_atom, 2)
    i = SQA.Index(h, :i, 2, h)
    j = SQA.Index(h, :j, 2, h)
    P = SQA.Transition(h, :σ, 2, 2)
    Pi = SQA.IndexedOperator(P, i)
    Pj = SQA.IndexedOperator(P, j)
    S = SQA.Σ(Pi, i)
    D = SQA.Σ(SQA.Σ(Pi * Pj, j), i)
    excluded = SQA.Σ(SQA.Σ(Pi * Pj, i, [j]), j)
    gi = SQA.IndexedVariable(:site_g, i)
    weighted = SQA.Σ(gi * Pi, i)
    parameters = Dict(SQA.substitute(gi, Dict(SQA.index_sym(i) => k)) => k for k in 1:2)
    @variables g
    b = QOB.tensor(QOB.NLevelBasis(2), QOB.NLevelBasis(2))
    for target in (b, (2, 2))
        expected = SQA.to_numeric(S, target, Dict(1 => [1, 2]))
        actual = SQA.to_numeric(S, target)
        @test Matrix(actual.data) ≈ Matrix(expected.data)
        @test Matrix(actual.data) ≈ ComplexF64[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2]
        @test Matrix(SQA.to_numeric(g * S, target; parameter = Dict(g => 3)).data) ≈ 3 * Matrix(expected.data)
        @test Matrix(SQA.to_numeric(excluded, target).data) ≈ Matrix(expected.data)^2 - Matrix(expected.data)
        @test_throws ArgumentError SQA.to_numeric(g * S, target)
        weighted_matrix = Matrix(SQA.to_numeric(weighted, target, Dict(1 => [1, 2]), SQA.NO_SUBS, parameters).data)
        @test sort(real.([weighted_matrix[k, k] for k in 1:4])) == [0, 1, 2, 3]
        @test Matrix(SQA.to_numeric(weighted, target; parameter = parameters).data) ≈ weighted_matrix
        @test Matrix(SQA.to_numeric(cos(weighted), target; parameter = parameters).data) ≈ cos(weighted_matrix)
        timed = SQA.to_numeric(g * D, target; time_parameter = Dict(g => t -> 1 + t))
        @test Matrix((target isa Tuple ? timed(nothing, 0.5) : QOB.sparse(QOB.static_operator(timed(0.5)))).data) ≈ 1.5 * Matrix(expected.data)^2
        @test Matrix(SQA.to_numeric(D, target).data) ≈ Matrix(expected.data)^2
        @test Matrix(SQA.to_numeric(cos(S), target).data) ≈ cos(Matrix(expected.data))
    end
end

struct StreamingProbe <: SQA.NumericBackend end
mutable struct StreamingCounter
    assignments::Int
end
SQA.numeric_num_subsystems(::StreamingProbe, n::Int) = n
SQA.numeric_zero(::StreamingProbe, ::Int) = StreamingCounter(0)
SQA.numeric_materialize(::StreamingProbe, counter::StreamingCounter, ::SQA.EagerOperator) = counter
function SQA.numeric_add_term!(counter::StreamingCounter, coefficient::ComplexF64, ops::Vector{SQA.Op}, ctx::SQA.NumericContext{StreamingProbe})
    @test coefficient == 1
    @test all(op -> SQA.index_slot(op.index) > 0, ops)
    counter.assignments += 1
    return counter
end

@testset "Indexed assignment accumulation" begin
    h = SQA.NLevelSpace(:probe_atom, 2)
    i = SQA.Index(h, :i, 3, h)
    P = SQA.IndexedOperator(SQA.Transition(h, :σ, 2, 2), i)
    ctx = SQA.NumericContext(StreamingProbe(), 3)
    result = SQA.static_operator(SQA.Σ(P, i), ctx, SQA.EagerOperator())
    @test result.assignments == 3
end
