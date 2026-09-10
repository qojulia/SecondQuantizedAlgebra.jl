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

    moving = Rotation(a, ω * t, t)
    @test transform(cos(A), moving) ==
        conjugate(cos(A), moving) + gauge_term(moving)
end
