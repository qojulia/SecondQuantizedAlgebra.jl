using SecondQuantizedAlgebra
using MutableArithmetics: add_mul, operate!!
using Symbolics: @variables
using Test
import MutableArithmetics
import SecondQuantizedAlgebra: constraint_pairs

@testset "Additive reductions" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = adjoint(a)

    @testset "sum and reduce preserve algebra" begin
        @variables x y
        terms = [x * a, y * ad, a' * a, 2 * a * ad]
        expected = foldl(+, terms)

        @test sum(terms) == expected
        @test reduce(+, terms) == expected
        @test iszero(sum([((x + y)^2) * a, -(x + y)^2 * a]))
        @test iszero(sum([a - a]))
        @test sum(y -> 2 * y, terms) == 2 * expected
        @test sum(identity, terms) == expected

        @test iszero(sum(x -> a, [a + ad]) - a)
        @test iszero(sum(x -> 2, [a + ad]) - 2)
        @test iszero(sum(x -> 0, [a + ad]))
    end

    @testset "indexed expressions retain public scope" begin
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        terms = [Σ(ai, i), Σ(aj' * ai, i)]

        reduced = sum(terms)
        @test get_indices(reduced) == get_indices(terms[1] + terms[2])
        @test Set(constraint_pairs(reduced)) ==
            Set(constraint_pairs(terms[1] + terms[2]))
        @test reduce(+, terms) == reduced
    end

    @testset "inputs remain values" begin
        terms = [a + ad, 2 * a, 3 * ad]
        before = terms[1]
        sum(terms)
        reduce(+, terms)
        @test isequal(terms[1], before)
        @test iszero(sum([a - a]))
    end

    @testset "inference at the public boundary" begin
        terms = [a + ad, 2 * a, 3 * ad]
        @test (@inferred sum(terms)) isa typeof(terms[1])
        @test (@inferred reduce(+, terms)) isa typeof(terms[1])
        @test (@inferred sum(y -> 2 * y, terms)) isa typeof(terms[1])
    end

    @testset "MutableArithmetics public integration" begin
        base = a + ad
        @test operate!!(+, base, a) == 2 * a + ad
        @test operate!!(+, base, 2) == 2 + a + ad
        @test operate!!(add_mul, base, 2, a) == 3 * a + ad

        # The rewrite macro is a public MutableArithmetics integration boundary.
        # In particular, this exercises bare operators and scalar terms, not only QAdd inputs.
        @test MutableArithmetics.@rewrite(a + ad + 2a) == 3a + ad
        @test MutableArithmetics.@rewrite(a + 2) == a + 2
        @test MutableArithmetics.@rewrite((a + ad) + 2) == a + ad + 2

        manual = MutableArithmetics.@rewrite sum([a + ad, 2a, 3ad])
        expected = sum([a + ad, 2a, 3ad])
        @test manual == expected
        @test expected == manual
        @test isequal(manual, expected)
        @test isequal(expected, manual)

        @test MutableArithmetics.promote_operation(sum, Vector{typeof(a + ad)}) === typeof(a + ad)
        @test MutableArithmetics.operate(sum, [a + ad, 2a, 3ad]) == expected
        @test MutableArithmetics.operate!!(sum, [a + ad, 2a, 3ad]) == expected
        @test (@inferred MutableArithmetics.operate(sum, [a + ad, 2a, 3ad])) isa typeof(expected)
        @test (@inferred MutableArithmetics.operate!!(sum, [a + ad, 2a, 3ad])) isa typeof(expected)
        @test MutableArithmetics.operate(sum, [a + ad, 2a, 3ad]; init = 2) == expected + 2
        @test MutableArithmetics.operate(sum, typeof(expected)[]; init = 2) == 2

        seeded = 1.0e16 * a
        values = [1.0 * a, -1.0e16 * a, 1.0 * a]
        @test iszero(MutableArithmetics.operate(sum, values; init = seeded) - a)
    end
end
