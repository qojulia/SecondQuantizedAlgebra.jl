using SecondQuantizedAlgebra
using Test

@testset "utils" begin
    @testset "operators" begin
        h1 = FockSpace(:s1)
        h2 = NLevelSpace(:s2, 2)
        h3 = FockSpace(:s3)
        h4 = PauliSpace(:s4)
        h5 = SpinSpace(:s5)
        h6 = PhaseSpace(:s6)
        h = h1 ⊗ h2 ⊗ h3 ⊗ h4 ⊗ h5 ⊗ h6
        a = Destroy(h, :a, 1)
        σ(i, j) = Transition(h, :σ, i, j)
        b = Destroy(h, :b, 3)
        pauli(i) = Pauli(h, :σP, i, 4)
        spin(i) = Spin(h, :S, i, 5)
        x = Position(h, :x, 6)
        p = Momentum(h, :p, 6)
        @testset "fundamental_operators" begin
            @test fundamental_operators(h1; names=[:a]) == [Destroy(h1, :a, 1)]
            @test fundamental_operators(h2; names=[:σ]) ==
                [Transition(h2, :σ, 1, 2), Transition(h2, :σ, 2, 2)]
            @test fundamental_operators(h4; names=[:σP]) == [Pauli(h4, :σP, i) for i in 1:3]
            @test fundamental_operators(h5; names=[:S]) == [Spin(h5, :S, i) for i in 1:3]
            @test fundamental_operators(h6; names=[(:x, :p)]) ==
                [Position(h6, :x, 1), Momentum(h6, :p, 1)]
            @test fundamental_operators(h; names=[:a, :σ, :b, :σP, :S, (:x, :p)]) == [
                a,
                σ(1, 2),
                σ(2, 2),
                b,
                pauli(1),
                pauli(2),
                pauli(3),
                spin(1),
                spin(2),
                spin(3),
                x,
                p,
            ]
        end
        @testset "embed" begin
            @test SecondQuantizedAlgebra.embed(h, Destroy(h1, :a, 1), 1) == a
            @test SecondQuantizedAlgebra.embed(h, Create(h1, :a, 1), 1) == a'
            @test SecondQuantizedAlgebra.embed(h, Transition(h2, :σ, 1, 2), 2) == σ(1, 2)
            @test SecondQuantizedAlgebra.embed(h, Pauli(h4, :σP, 1), 4) == pauli(1)
            @test SecondQuantizedAlgebra.embed(h, Spin(h5, :S, 1), 5) == spin(1)
            @test SecondQuantizedAlgebra.embed(h, Position(h6, :x, 1), 6) == x
            @test SecondQuantizedAlgebra.embed(h, Momentum(h6, :p, 1), 6) == p
        end
        @testset "unique_ops" begin
            ops = [
                a,
                a',
                σ(1, 2),
                σ(2, 2),
                b,
                σ(2, 1),
                pauli(1),
                pauli(2),
                pauli(3),
                spin(1),
                spin(2),
                spin(3),
                x,
                p,
            ]
            @test isequal(
                unique_ops(ops),
                [
                    a,
                    σ(1, 2),
                    σ(2, 2),
                    b,
                    pauli(1),
                    pauli(2),
                    pauli(3),
                    spin(1),
                    spin(2),
                    spin(3),
                    x,
                    p,
                ],
            )
            @test isequal(
                unique_ops!(ops),
                [
                    a,
                    σ(1, 2),
                    σ(2, 2),
                    b,
                    pauli(1),
                    pauli(2),
                    pauli(3),
                    spin(1),
                    spin(2),
                    spin(3),
                    x,
                    p,
                ],
            )
        end
    end
end
