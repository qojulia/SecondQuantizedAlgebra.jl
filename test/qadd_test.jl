using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QMul, QAdd, QSym, QField
using Test

@testset "QAdd" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QSym + QSym" begin
        s = a + ad
        @test s isa QAdd
        @test length(s.arguments) == 2
        @test all(x -> x isa QMul, s.arguments)
    end

    @testset "QMul + QMul" begin
        s = (2 * a) + (3 * ad)
        @test s isa QAdd
        @test length(s.arguments) == 2
    end

    @testset "QMul + QSym" begin
        s1 = (2 * a) + ad
        @test s1 isa QAdd
        @test length(s1.arguments) == 2

        s2 = ad + (2 * a)
        @test s2 isa QAdd
    end

    @testset "QAdd + QMul" begin
        s = (a + ad) + (3 * a)
        @test s isa QAdd
        @test length(s.arguments) == 3
    end

    @testset "QAdd + QAdd" begin
        s = (a + ad) + (a + ad)
        @test s isa QAdd
        @test length(s.arguments) == 4
    end

    @testset "QField + Number" begin
        s1 = a + 5
        @test s1 isa QAdd
        @test length(s1.arguments) == 2

        s2 = 5 + a
        @test s2 isa QAdd

        s3 = (a + ad) + 3
        @test s3 isa QAdd
        @test length(s3.arguments) == 3

        s4 = 3 + (a + ad)
        @test s4 isa QAdd
    end

    @testset "Subtraction" begin
        s = a - ad
        @test s isa QAdd
        @test length(s.arguments) == 2

        s2 = a - 3
        @test s2 isa QAdd
    end

    @testset "QAdd * Number" begin
        s = (a + ad) * 3
        @test s isa QAdd
        @test all(x -> x.arg_c == 3, s.arguments)
    end

    @testset "QAdd * QSym (distributes)" begin
        s = (a + ad) * a
        @test s isa QAdd
        @test length(s.arguments) == 2
        @test all(x -> length(x.args_nc) == 2, s.arguments)
    end

    @testset "QAdd * QMul (distributes)" begin
        s = (a + ad) * (2 * a)
        @test s isa QAdd
        @test length(s.arguments) == 2
    end

    @testset "QAdd * QAdd (distributes)" begin
        s = (a + ad) * (a + ad)
        @test s isa QAdd
        @test length(s.arguments) == 4
    end

    @testset "Equality and hashing" begin
        s1 = a + ad
        s2 = a + ad
        @test isequal(s1, s2)
        @test hash(s1) == hash(s2)
    end

    @testset "Adjoint" begin
        s = a + 2 * ad
        sd = s'
        @test sd isa QAdd
        @test length(sd.arguments) == 2
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(QAdd)
    end

    @testset "Type stability" begin
        # QSym + QSym
        @inferred a + ad
        # QMul + QMul
        @inferred (2 * a) + (3 * ad)
        # QMul + QSym / QSym + QMul
        @inferred (2 * a) + ad
        @inferred ad + (2 * a)
        # QAdd + QMul / QMul + QAdd
        @inferred (a + ad) + (3 * a)
        @inferred (3 * a) + (a + ad)
        # QAdd + QSym / QSym + QAdd
        @inferred (a + ad) + a
        @inferred a + (a + ad)
        # QAdd + QAdd
        @inferred (a + ad) + (a + ad)
        # QField + Number / Number + QField
        @inferred a + 5
        @inferred 5 + a
        @inferred (2 * a) + 5
        @inferred 5 + (2 * a)
        @inferred (a + ad) + 3
        @inferred 3 + (a + ad)
        # Subtraction
        @inferred a - ad
        @inferred a - 3
        @inferred 3 - a
        # Negation
        @inferred -(a + ad)
        # Distributive multiplication
        @inferred (a + ad) * 3
        @inferred 3 * (a + ad)
        @inferred (a + ad) * a
        @inferred a * (a + ad)
        @inferred (a + ad) * (2 * a)
        @inferred (2 * a) * (a + ad)
        @inferred (a + ad) * (a + ad)
        # Division
        @inferred (a + ad) / 2
        @inferred (a + ad) / 2.0
        # Adjoint
        @inferred adjoint(a + 2 * ad)
        # Equality / hashing
        @inferred isequal(a + ad, a + ad)
        @inferred hash(a + ad, UInt(0))
    end

    @testset "Allocations" begin
        s = a + ad
        m_2a = 2 * a
        m_3ad = 3 * ad

        # Warmup
        a + ad; m_2a + m_3ad; s + m_2a; s + s; s * 3; s * a; s * s

        # QSym + QSym
        @test @allocations(a + ad) <= 200
        # QMul + QMul
        @test @allocations(m_2a + m_3ad) <= 5
        # QAdd + QMul
        @test @allocations(s + m_2a) <= 100
        # QAdd + QAdd
        @test @allocations(s + s) <= 15
        # QAdd * Number
        @test @allocations(s * 3) <= 500
        # QAdd * QSym (distributive)
        @test @allocations(s * a) <= 15
        # QAdd * QAdd (distributive)
        @test @allocations(s * s) <= 10000
        # Negation
        @test @allocations(-s) <= 500
        # Adjoint
        @test @allocations(adjoint(s)) <= 200
    end
end
