using SecondQuantizedAlgebra
using Test

@testset "utils" begin
    @testset "operators" begin
        h1 = FockSpace(:s1)
        h2 = NLevelSpace(:s2, 2)
        h3 = FockSpace(:s3)
        h = h1 ⊗ h2 ⊗ h3
        a = Destroy(h, :a, 1)
        σ(i, j) = Transition(h, :σ, i, j)
        b = Destroy(h, :b, 3)
        @testset "fundamental_operators" begin
            @test fundamental_operators(h1; names=[:a]) == [Destroy(h1, :a, 1)]
            @test fundamental_operators(h2; names=[:σ]) ==
                [Transition(h2, :σ, 1, 2), Transition(h2, :σ, 2, 2)]
            @test fundamental_operators(h; names=[:a, :σ, :b]) == [
                Destroy(h, :a, 1),
                Transition(h, :σ, 1, 2),
                Transition(h, :σ, 2, 2),
                Destroy(h, :b, 3),
            ]
        end
        @testset "embed" begin
            @test SecondQuantizedAlgebra.embed(h, Destroy(h1, :a, 1), 1) == a
            @test SecondQuantizedAlgebra.embed(h, Create(h1, :a, 1), 1) == a'
            @test SecondQuantizedAlgebra.embed(h, Transition(h2, :σ, 1, 2), 2) == σ(1, 2)
        end
        @testset "unique_ops" begin
            ops = [a, a', σ(1, 2), σ(2, 2), b, σ(2, 1)]
            @test isequal(unique_ops(ops), [a, σ(1, 2), σ(2, 2), b])
            @test isequal(unique_ops!(ops), [a, σ(1, 2), σ(2, 2), b])
        end
    end
end
