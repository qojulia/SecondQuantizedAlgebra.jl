using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Transmon/Josephson Taylor workflow" begin
    h = PhaseSpace(:transmon) ⊗ FockSpace(:resonator)
    @qnumbers φ::Position(h, 1) n::Momentum(h, 1) a::Destroy(h, 2)

    hb = FockSpace(:transmon_boson) ⊗ FockSpace(:resonator)
    @qnumbers b::Destroy(hb, 1)

    @variables E_C::Real E_J::Real ω::Real g::Real φ_zpf::Real n_zpf::Real

    H = 4 * E_C * n^2 - E_J * cos(φ) + ω * a' * a + g * n * (a' + a)
    Hb = substitute(
        H,
        Dict(
            φ => φ_zpf * (b + b'),
            n => im * n_zpf * (b' - b),
        ),
    )
    H4 = @inferred taylor(Hb, 0:4)

    X2_normal = 1 + b'^2 + 2 * b' * b + b^2
    X4_normal =
        3 + b'^4 + 4 * b'^3 * b + 6 * b'^2 * b^2 + 4 * b' * b^3 + b^4 +
        6 * b'^2 + 12 * b' * b + 6 * b^2
    n2_normal = 1 + 2 * b' * b - b'^2 - b^2

    expected =
        4 * E_C * n_zpf^2 * n2_normal -
        E_J * (
            1 - (1 // 2) * φ_zpf^2 * X2_normal +
                (1 // 24) * φ_zpf^4 * X4_normal
        ) +
        ω * a' * a +
        im * g * n_zpf * (b' * a' + b' * a - b * a' - b * a)

    @test H4 isa SecondQuantizedAlgebra.QAdd
    @test iszero(simplify(H4 - expected))
end
