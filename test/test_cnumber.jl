using SymbolicUtils, Symbolics, SecondQuantizedAlgebra, Test

@testset "Parameter" begin
    @rnumbers ω
    @cnumbers G G1
    @rnumbers r1 r2

    @test first(typeof(G).parameters) == Complex{Real}
    @test first(typeof(ω).parameters) == Real

    @test SymbolicUtils.getmetadata(ω, Symbolics.VariableSource) == (:RealParameter, :ω)
    @test SymbolicUtils.getmetadata(G, Symbolics.VariableSource) == (:Parameter, :G)

    @testset "rnumbers_cnumbers" begin
        @test isequal(rnumbers(:r1, :r2), (r1, r2))
        @test isequal(rnumbers("r1"), (r1,))
        @test isequal(rnumber("r1"), (r1))
        @test isequal(rnumber(:r1), (r1))
        @test isequal(cnumber("c"), cnumber(:c))
    end

    @testset "one_zero" begin
        @test one(ω) == 1.0
        @test one(G) == 1.0 + 0.0im
        @test zero(ω) == 0.0
        @test zero(G) == 0.0 + 0.0im
    end

    @testset "adjoint" begin
        @test isequal(adjoint(ω), ω)
        @test isequal(adjoint(3ω), 3ω)
        @test isequal(adjoint(G), conj(G))
        @test_broken isequal(conj(3*G), 3*conj(G))
        @test isequal(simplify(adjoint(G*G1)), conj(G*G1))
    end
end
