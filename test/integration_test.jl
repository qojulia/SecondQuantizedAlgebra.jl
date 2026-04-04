using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
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

    @testset "Distinct modes on same space don't commute" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        a = Destroy(h1, :a)
        b = Destroy(h2, :b)
        # a and b have same space_index=1 but different names
        m = a * b'
        result = normal_order(m)
        @test length(result.arguments) == 1  # no commutation applied
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

        # numeric_average on QAdd
        expr = a + a' * a
        @test numeric_average(expr, ψ) ≈ α + abs2(α)

        # to_numeric with scalar QMul (empty args_nc)
        @test to_numeric(QMul(3, QSym[]), b) == 3 * one(b)
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

        # Printing with symbolic prefactors should not crash
        @test repr(g * a) isa String
        @test repr(g * a + ω * a') isa String
    end

    @testset "Division with non-integer" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test (a / 2.0) isa QMul
        @test (a / 2.0).arg_c ≈ 0.5

        @test (a / 2) isa QMul
        @test (a / 2).arg_c == 1 // 2

        @test ((a + a') / 2.0) isa QAdd
    end

    @testset "Power edge cases" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^0 = identity (scalar QMul)
        m0 = a^0
        @test m0 isa QMul
        @test isempty(m0.args_nc)
        @test m0.arg_c == 1
    end

    @testset "Adjoint with complex prefactor" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        m = (2.0 + 3.0im) * a' * a
        ma = m'
        @test ma.arg_c ≈ conj(2.0 + 3.0im)
    end

    @testset "Higher-order normal ordering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^2 * (a')^2 should produce multiple terms
        m = a * a * a' * a'
        result = simplify(normal_order(m))
        @test result isa QAdd
        @test length(result.arguments) >= 3
    end

    @testset "== consistency" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test QMul(1, QSym[a]) == QMul(1, QSym[a])
        @test (a + a') == (a + a')
    end
end
