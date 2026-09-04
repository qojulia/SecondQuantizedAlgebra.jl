using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Public substitution" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables x::Real y::Real

    @testset "scalar substitutions" begin
        @test substitute(42, Dict(x => 0)) == 42
        @test iszero(substitute(x * a, Dict(x => 0)))
        @test iszero(substitute(x * a, Dict(x => 1)) - a)
        @test iszero(simplify(substitute(x * a, Dict(x => y)) - y * a))
        @test iszero(
            simplify(
                substitute(x * (a + a'), Dict(x => y)) - y * (a + a'),
            )
        )
    end

    @testset "operator substitutions" begin
        @test iszero(substitute(x * a, Dict(a => 0)))
        @test iszero(simplify(substitute(x * a, Dict(a => y)) - x * y))
        @test iszero(simplify(substitute(x * a, Dict(a => a')) - x * a'))
        @test iszero(substitute(a, Dict(a' => a)) - a')
        @test iszero(
            substitute(a', Dict(a => 2 * a); replace_adjoint = false) - a',
        )
        @test_throws ArgumentError substitute(a, Dict(a => :unsupported); replace_adjoint = false)
    end

    @testset "replacement expressions are expanded once" begin
        h3 = FockSpace(:one) ⊗ FockSpace(:two) ⊗ FockSpace(:three)
        a1 = Destroy(h3, :a, 1)
        a2 = Destroy(h3, :a, 2)
        a3 = Destroy(h3, :a, 3)
        @variables g₂::Real g₃::Real
        replacement = g₂ * a1 + g₃ * a3

        @test iszero(simplify(substitute(a1, Dict(a1 => replacement)) - replacement))
        @test iszero(
            simplify(
                substitute(a1 * a2, Dict(a1 => g₂ * a2 + g₃ * a3)) -
                    (g₂ * a2 + g₃ * a3) * a2,
            )
        )
        @test iszero(
            simplify(
                substitute(a1', Dict(a1 => g₂ * a2)) - g₂ * a2',
            )
        )
        @test iszero(
            substitute(a1', Dict(a1 => g₂ * a2); replace_adjoint = false) - a1',
        )
    end

    @testset "substitution also handles composite spaces" begin
        h2 = h ⊗ FockSpace(:second)
        b = Destroy(h2, :b, 2)
        @test iszero(substitute(x * a' * b, Dict(x => 2)) - 2 * a' * b)

        @variables gc::Number
        @test iszero(substitute((gc * a)', Dict(gc => 0)))
        @test isequal(substitute((gc * a)', Dict(gc => 2)), 2 * a')
    end

    @testset "public substitution is inferable" begin
        @inferred substitute(x * a, Dict(x => 2))
        @inferred substitute(a, Dict(a => a'))
    end
end
