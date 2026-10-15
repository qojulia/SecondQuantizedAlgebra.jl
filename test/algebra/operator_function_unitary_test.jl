using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Formal operator functions under exact unitary transforms" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    @variables θ::Real ω::Real t::Real

    A = a + a'
    U = Rotation(a, θ)
    image = conjugate(A, U)

    @test @inferred(conjugate(cos(A), U)) == cos(image)
    @test @inferred(conjugate(sin(A), U)) == sin(image)
    @test @inferred(conjugate(expim(A), U)) == expim(image)
    @test conjugate(conjugate(cos(A), U), inv(U)) == cos(A)
    @test transform(cos(A), U) == conjugate(cos(A), U)

    composite = cos(A) + sin(A)
    @test conjugate(composite, U) == cos(image) + sin(image)

    moving = Rotation(a, ω * t, t)
    @test transform(cos(A), moving) ==
        conjugate(cos(A), moving) + gauge_term(moving)

    h2 = FockSpace(:left) ⊗ FockSpace(:right)
    @qnumbers a2::Destroy(h2, 1) b2::Destroy(h2, 2)
    A2 = a2 + a2'
    mixing = Rotation(a2, b2, θ)
    unit_coeff = cos(θ)^2 + sin(θ)^2
    weighted = unit_coeff * cos(A2)
    @test conjugate(weighted, mixing) == cos(conjugate(A2, mixing))
end
