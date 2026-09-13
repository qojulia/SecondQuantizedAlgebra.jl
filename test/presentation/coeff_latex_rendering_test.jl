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

@testset "Coeff arrays render as LaTeX arrays" begin
    @variables E γ ω

    M = fill(to_cnum(0), 3, 3)
    M[2, 2] = to_cnum(0.25 * E^4 * γ / ω^4)
    matrix_tex = raw"""
    \begin{equation}
    \left[
    \begin{array}{ccc}
    0 & 0 & 0 \\
    0 & \frac{0.25 ~ E^{4} ~ \gamma}{\omega^{4}} & 0 \\
    0 & 0 & 0 \\
    \end{array}
    \right]
    \end{equation}
    """
    @test String(latexify(M)) == matrix_tex
    @test repr(MIME"text/latex"(), M) == matrix_tex

    # One column per prefactor branch: pure imaginary, mixed, and native zero.
    v = [to_cnum(2im), to_cnum(E + im * γ), to_cnum(0)]
    vector_tex = raw"""
    \begin{equation}
    \left[
    \begin{array}{c}
    2\mathit{i} \\
    E + \gamma ~ \mathit{i} \\
    0 \\
    \end{array}
    \right]
    \end{equation}
    """
    @test String(latexify(v)) == vector_tex
    @test repr(MIME"text/latex"(), v) == vector_tex
end
