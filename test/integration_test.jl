using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: @variables
using QuantumOpticsBase
using Test
import SecondQuantizedAlgebra: simplify

@testset "integration" begin
    @testset "Basic operator creation and algebra" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'

        @test hash(a) != hash(ad)
        @test isequal(a, ad')

        # Lazy multiplication
        @test a * ad isa QMul
        @test ad * a isa QMul
    end

    @testset "Commutation via normal_order" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a * a† normal-ordered → a†a + 1
        m = a * a'
        result = simplify(normal_order(m))
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "Numeric conversion round-trip" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
        @test to_numeric(a' * a, b) == create(b) * destroy(b)

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end

    @testset "Multi-space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

        # Operators on different spaces — normal_order doesn't commute them
        m = a * b'
        result = normal_order(m)
        @test length(result.arguments) == 1

        # Same-space commutation still works
        m2 = a * a'
        result2 = simplify(normal_order(m2))
        @test length(result2.arguments) == 2
    end

    @testset "Symbolic parameters" begin
        @variables g ω
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        expr = g * a' * a + ω * a' * a
        result = simplify(expr)
        @test length(result.arguments) == 1
    end
end
