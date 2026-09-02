using SecondQuantizedAlgebra
using MutableArithmetics: add_mul, operate!!
using Symbolics: @variables
using Test
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
    end
end
