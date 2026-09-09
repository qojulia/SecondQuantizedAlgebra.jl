using SecondQuantizedAlgebra
using Test
using Symbolics: @variables

@testset "Public unitary refusal boundaries" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    @variables ω t
    nonharmonic = sqrt(1 + t^2)

    @test_throws ArgumentError DisplacementFrame(
        a, ω * a' * a + nonharmonic * (a + a'), t,
    )
    @test_throws ArgumentError DisplacementFrame(
        x, p, (ω / 2) * (x^2 + p^2) + nonharmonic * x, t,
    )

    @test_throws ArgumentError UnitaryTransform((2 * x^2 + 3 * p^2) / 2, 1)

    collective = CollectiveNLevelSpace(:collective, 2)
    S11 = CollectiveTransition(collective, :S, 1, 1)
    @test_throws ArgumentError UnitaryTransform(S11, 1)

    moving = Rotation(a, ω * t, t)
    unchanged = substitute(moving, Dict(:unused => 1))
    @test iszero(simplify(conjugate(a, unchanged) - conjugate(a, moving)))
    @test iszero(simplify(gauge_term(unchanged) - gauge_term(moving)))
end
