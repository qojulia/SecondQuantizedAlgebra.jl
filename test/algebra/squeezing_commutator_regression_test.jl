using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Squeezing commutator cancels exact rational coefficients" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables Δ λ

    H = Δ * a' * a + (λ / 2) * (a' * a' + a * a)
    rhs = commutator(im * H, a * a)
    expected = -2im * Δ * a * a - im * λ * (2 * a' * a + 1)

    # Exact public-API equality: this must fail if a zero-coefficient fourth-order term survives.
    @test iszero(rhs - expected)
end
