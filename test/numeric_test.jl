using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

@testset "numeric conversion" begin
    @testset "Single space — basic" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
    end

    @testset "QMul" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a' * a, b) == create(b) * destroy(b)
        @test to_numeric(2 * a, b) == 2 * destroy(b)
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        expected = destroy(b) + create(b)
        @test result == expected
    end

    @testset "Scalar" begin
        b = FockBasis(7)
        @test to_numeric(3, b) == 3 * one(b)
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)

        a1_num = to_numeric(a1, b)
        @test a1_num isa LazyTensor
    end

    @testset "numeric_average" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end
end
