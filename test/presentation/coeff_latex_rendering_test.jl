using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: expim, to_cnum

@testset "Standalone Coeff LaTeX" begin
    @variables g κ θ

    cases = [
        (to_cnum(2im), L"2\mathit{i}"),
        (to_cnum(g), L"g"),
        (to_cnum(im * g), L"g ~ \mathit{i}"),
        (to_cnum(g + im * κ), L"g + \kappa ~ \mathit{i}"),
        (expim(θ), L"e^{i \theta}"),
    ]

    for (c, expected) in cases
        @test c isa SecondQuantizedAlgebra.Coeff
        @test latexify(c) == expected
        @test repr(MIME"text/latex"(), c) == expected
    end
end
