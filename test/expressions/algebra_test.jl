using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

@testset "Algebraic expressions" begin
    h = FockSpace(:fock)
    a = Destroy(h, :a)
    ad = adjoint(a)

    @testset "linear arithmetic" begin
        @test iszero(simplify((a + ad) - a - ad))
        @test iszero(simplify((a + ad) + (a + ad) - 2 * a - 2 * ad))
        @test iszero(simplify((a + ad) - (ad + a)))
        @test iszero(simplify(a - a))
        @test iszero(simplify(3 * a - a - a - a))
    end

    @testset "distribution and collection" begin
        left = (a + ad) * (a - ad)
        right = a * a - a * ad + ad * a - ad * ad
        @test iszero(simplify(left - right))

        @variables g ω
        @test iszero(simplify(g * a + ω * a - (g + ω) * a))
        @test iszero(simplify((2 * a + 3 * ad) + (4 * a - ad) - 6 * a - 2 * ad))
    end

    @testset "eager canonicalization is observable through public operations" begin
        # The commutation relation is already applied by `*`; the explicit
        # finalizers therefore preserve the same mathematical expression.
        expected = ad * a + 1
        @test iszero(a * ad - expected)
        @test iszero(normal_order(a * ad) - expected)
        @test iszero(simplify(a * ad) - expected)
        @test isequal(simplify(simplify(a * ad)), simplify(a * ad))
    end

    @testset "powers, division, and adjoints" begin
        @test iszero(a^0 - 1)
        @test iszero(a^3 - a * a * a)
        @test iszero(a / 2 - (1 // 2) * a)
        @test iszero(
            adjoint(adjoint((a + 2im * ad) * (a - ad))) -
                (a + 2im * ad) * (a - ad)
        )
        @test iszero(adjoint(a * ad) - adjoint(ad) * adjoint(a))
    end

    @testset "public QAdd arithmetic boundaries" begin
        @variables z::Number
        @test iszero(simplify((a + ad) + 3 - (3 + a + ad)))
        @test iszero(simplify((a + ad) * 3 - 3 * (a + ad)))
        @test iszero(simplify((a + ad) * a - (a * a + ad * a)))
        @test iszero(
            simplify(
                (a + ad) * (a - ad) -
                    (a * a - a * ad + ad * a - ad * ad)
            )
        )
        @test iszero(simplify((a + ad) / 2 - (1 // 2) * (a + ad)))
        @test iszero(simplify(+a - a))
        @test iszero(simplify((z + a) - (a + z)))

        expression = (2 + 3im) * ad * a
        @test iszero(adjoint(expression) - (2 - 3im) * ad * a)
    end

    @testset "independent sites commute" begin
        h2 = FockSpace(:left) ⊗ FockSpace(:right)
        left = Destroy(h2, :left, 1)
        right = Destroy(h2, :right, 2)
        @test iszero(left * right - right * left)
        @test iszero(left' * right - right * left')
        @test iszero(commutator(left' * left, right' * right))
    end

    @testset "public reductions" begin
        terms = [a + ad, 2 * a, 3 * ad, a * ad, ad * a]
        manual = terms[1] + terms[2] + terms[3] + terms[4] + terms[5]
        @test iszero(sum(terms) - manual)
        @test iszero(reduce(+, terms) - manual)

        @test iszero(zero(a + ad))
        @test iszero(sum([a - a]))
        @test iszero(sum(x -> x * a, [a, ad]) - (a * a + ad * a))
    end

    @testset "public entry points remain inferable" begin
        @inferred(a + ad)
        @inferred(a * ad)
        @inferred(commutator(a, ad))
        @inferred(simplify(a + ad))
        @inferred(normal_order(a * ad))
        @inferred(adjoint(a + 2 * ad))
    end
end
