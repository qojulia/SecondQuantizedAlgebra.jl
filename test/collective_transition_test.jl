using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "collective transition" begin
    h = NLevelSpace(:atom, 3)
    S(i, j) = CollectiveTransition(h, :S, i, j)

    @testset "basic properties" begin
        @test S(1, 2)' == S(2, 1)
        @test S(2, 2)' == S(2, 2)
    end

    @testset "ordered products stay unchanged" begin
        @test isequal(S(2, 1) * S(1, 3), SecondQuantizedAlgebra.QMul(1, [S(2, 1), S(1, 3)]))
        @test isequal(S(2, 3) * S(2, 1), SecondQuantizedAlgebra.QMul(1, [S(2, 3), S(2, 1)]))
        @test isequal(S(2, 2) * S(2, 2), SecondQuantizedAlgebra.QMul(1, [S(2, 2), S(2, 2)]))
    end

    @testset "rewrite rule" begin
        @test isequal(S(1, 2) * S(2, 3), S(1, 3) + S(2, 3) * S(1, 2))
        @test isequal(S(1, 2) * S(1, 3), S(1, 3) * S(1, 2))
        @test isequal(S(1, 2) * S(3, 1), -S(3, 2) + S(3, 1) * S(1, 2))
        @test isequal(S(1, 2) * S(1, 1), simplify(- S(1, 2) + S(1, 1) * S(1, 2)))
    end

    @testset "product spaces still commute by acts_on ordering" begin
        h2 = NLevelSpace(:atom2, 3)
        hp = h ⊗ h2
        S1(i, j) = CollectiveTransition(hp, :S1, i, j, 1)
        S2(i, j) = CollectiveTransition(hp, :S2, i, j, 2)

        @test isequal(S2(2, 1) * S1(1, 3), S1(1, 3) * S2(2, 1))
    end

    @testset "two collective transitions with one destroy operator" begin
        hc = FockSpace(:cavity)
        hp = hc ⊗ h
        a = Destroy(hp, :a, 1)
        Sp(i, j) = CollectiveTransition(hp, :S, i, j, 2)

        @testset "commuting transitions" begin
            @test isequal(Sp(1, 2) * Sp(2, 3) * a, a*Sp(1, 3) + Sp(2, 3) * Sp(1, 2) * a)
        end

        @testset "noncommuting transitions" begin
            @test isequal(Sp(1, 2) * Sp(1, 3) * a, Sp(1, 3) * Sp(1, 2) * a)
            @test isequal(simplify(Sp(1, 2) * Sp(1, 3) * a - (a * Sp(1, 3) * Sp(1, 2))), 0)
        end
    end
end
