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
    sN = sin(N)

    @test repr(cN) == "cos(a' * a)"
    @test repr(sN) == "sin(a' * a)"
    @test repr(2 * cN) == "2 * cos(a' * a)"
    @test repr(a * cN) == "a * cos(a' * a)"
    @test repr(cN * a) == "cos(a' * a) * a"
    @test repr((a + a') * cN) == "(a + a') * cos(a' * a)"
    @test repr(expim(-θ * N)) == "exp(im*(-θ * a' * a))"

    formal_sum = sN + cN
    formal_difference = sN - cN
    @test repr(formal_sum) == "sin(a' * a) + cos(a' * a)"
    @test repr(formal_difference) == "sin(a' * a) - cos(a' * a)"
    @test repr(zero(SecondQuantizedAlgebra.QExpr)) == "0"
    @test repr(one(SecondQuantizedAlgebra.QExpr)) == "1"
    @test repr(-cN) == "-cos(a' * a)"

    rational = repr((1 // 2) * cN)
    @test occursin("1//2", rational)
    @test !occursin("0.5", rational)

    raw_prefactor = sqrt(θ + 1) * cN
    @test occursin("sqrt", repr(raw_prefactor))
    @test occursin("cos", repr(raw_prefactor))

    grouped_formal = formal_sum * a
    @test repr(grouped_formal) == "(sin(a' * a) + cos(a' * a)) * a"

    @test latexify(cN) == L"\cos\left( a^{\dagger}a \right)"
    @test latexify(sN) == L"\sin\left( a^{\dagger}a \right)"
    @test latexify(2 * cN) == L"2 \cos\left( a^{\dagger}a \right)"
    @test latexify((a + a') * cN) ==
        L"\left( a + a^{\dagger} \right) \cos\left( a^{\dagger}a \right)"
    @test latexify(formal_sum) ==
        L"\sin\left( a^{\dagger}a \right) + \cos\left( a^{\dagger}a \right)"
    @test latexify(formal_difference) ==
        L"\sin\left( a^{\dagger}a \right) - \cos\left( a^{\dagger}a \right)"
    @test latexify(zero(SecondQuantizedAlgebra.QExpr)) == L"0"
    @test latexify(-cN) == L"-\cos\left( a^{\dagger}a \right)"
    @test occursin("\\frac", String(latexify((1 // 2) * cN)))
    @test occursin("\\sqrt", String(latexify(raw_prefactor)))
    @test occursin("\\left", String(latexify(grouped_formal)))
    # Latexify spaces unary minus as an operator. The exact contract here is that the
    # complete negative operator argument remains inside the exponential grouping.
    @test latexify(expim(-θ * N)) == L"e^{i\left(  - \theta a^{\dagger}a \right)}"
    @test repr(MIME"text/latex"(), cN) == latexify(cN)
end
