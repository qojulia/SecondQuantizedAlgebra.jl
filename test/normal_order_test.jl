using SecondQuantizedAlgebra
using Test

@testset "normal_order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Already normal-ordered" begin
        m = ad * a
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
        @test isequal(result.arguments[1], m)
    end

    @testset "Single commutation: a*a† → a†a + 1" begin
        m = a * ad  # stored as QMul(1, [Create, Destroy]) due to canonical sort
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 2
    end

    @testset "QSym passthrough" begin
        result = normal_order(a)
        @test result isa QAdd
        @test length(result.arguments) == 1
    end

    @testset "QAdd — normal-orders each term" begin
        expr = a * ad + ad * a
        result = normal_order(expr)
        @test result isa QAdd
    end

    @testset "Different spaces don't commute" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result.arguments) == 1
    end

    @testset "Return type is always QAdd" begin
        @test normal_order(a) isa QAdd
        @test normal_order(2 * a) isa QAdd
        @test normal_order(a + ad) isa QAdd
    end
end
