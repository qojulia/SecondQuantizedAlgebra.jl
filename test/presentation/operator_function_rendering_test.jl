using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test

@testset "Formal operator expression rendering" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    @variables θ::Real

    N = a' * a
    cN = cos(N)

    @test repr(cN) == "cos(a' * a)"
    @test repr(2 * cN) == "2 * cos(a' * a)"
    @test repr(a * cN) == "a * cos(a' * a)"
    @test repr(cN * a) == "cos(a' * a) * a"
    @test repr((a + a') * cN) == "(a + a') * cos(a' * a)"
    @test repr(expim(-θ * N)) == "exp(im*(-θ * a' * a))"

    @test latexify(cN) == L"\cos\left( a^{\dagger}a \right)"
    @test latexify(2 * cN) == L"2 \cos\left( a^{\dagger}a \right)"
    @test latexify((a + a') * cN) ==
        L"\left( a + a^{\dagger} \right) \cos\left( a^{\dagger}a \right)"
    @test latexify(expim(-θ * N)) == L"e^{i\left( -\theta a^{\dagger}a \right)}"
    @test repr(MIME"text/latex"(), cN) == latexify(cN)
end
