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

    @testset "NLevel numeric" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        b = NLevelBasis(3)
        @test to_numeric(σ12, b) == transition(b, 1, 2)
        @test to_numeric(σ12', b) == transition(b, 2, 1)
    end

    @testset "Pauli numeric" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)
        b = SpinBasis(1 // 2)
        @test to_numeric(σx, b) == sigmax(b)
        @test to_numeric(σy, b) == sigmay(b)
        @test to_numeric(σz, b) == sigmaz(b)
    end

    @testset "Spin numeric" begin
        h = SpinSpace(:s, 5 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        b = SpinBasis(5 // 2)
        @test to_numeric(Sx, b) == 0.5 * sigmax(b)
        @test to_numeric(Sy, b) == 0.5 * sigmay(b)
        @test to_numeric(Sz, b) == 0.5 * sigmaz(b)
    end

    @testset "Composite NLevel + Fock" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        @qnumbers a::Destroy(h, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        bf = FockBasis(3)
        bn = NLevelBasis(3)
        bc = bf ⊗ bn
        @test to_numeric(σ12, bc) isa LazyTensor
    end

    # TODO: to_numeric is not type-stable due to QuantumOpticsBase dispatch.
    # Fix upstream in QuantumOpticsBase.
    # @testset "Type stability" begin
    #     h = FockSpace(:fock)
    #     @qnumbers a::Destroy(h)
    #     b = FockBasis(7)
    #     @inferred to_numeric(a, b)
    #     @inferred to_numeric(a', b)
    # end

    @testset "Allocations — to_numeric" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        # Warmup
        to_numeric(a, b)
        to_numeric(a', b)
        to_numeric(a' * a, b)

        # Single operators should be bounded
        @test @allocations(to_numeric(a, b)) < 50
        @test @allocations(to_numeric(a', b)) < 50
        @test @allocations(to_numeric(a' * a, b)) < 100
    end
end
