using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: needs_pf_parens, show_display

@testset "Generic symbolic number printer compatibility" begin
    @variables x y

    @test sprint(io -> show_display(io, x)) == "x"
    @test !needs_pf_parens(1)
    @test !needs_pf_parens(x)
    @test needs_pf_parens(x + y)
    @test !needs_pf_parens(x * y)
end
