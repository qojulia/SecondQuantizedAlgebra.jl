using SymbolicUtils, Symbolics, SecondQuantizedAlgebra, Test

@testset "Parameter" begin
    @rnumbers ω
    @cnumbers G

    @test first(typeof(G).parameters) == Complex{Real}
    @test first(typeof(ω).parameters) == Real

    @test SymbolicUtils.getmetadata(ω, Symbolics.VariableSource) == (:RealParameter, :ω)
    @test SymbolicUtils.getmetadata(G, Symbolics.VariableSource) == (:Parameter, :G)
end
