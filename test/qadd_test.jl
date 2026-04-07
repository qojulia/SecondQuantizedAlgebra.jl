using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QAdd, QSym, QField, sorted_arguments, CNum
using Test

@testset "QAdd" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QAdd + QAdd (single-term)" begin
        s = (2 * a) + (3 * ad)
        @test s isa QAdd
        @test length(s) == 2
    end

    @testset "QAdd + QSym" begin
        s1 = (2 * a) + ad
        @test s1 isa QAdd
        @test length(s1) == 2

        s2 = ad + (2 * a)
        @test s2 isa QAdd
    end

    @testset "QAdd + QAdd (merge)" begin
        s = (a + ad) + (3 * a)
        @test s isa QAdd
        @test length(s) == 2  # a and 3a auto-collect to 4a, plus a†
    end

    @testset "QAdd + QAdd" begin
        s = (a + ad) + (a + ad)
        @test s isa QAdd
        @test length(s) == 2  # auto-collects to 2a + 2a†
    end

    @testset "QField + Number" begin
        s1 = a + 5
        @test s1 isa QAdd
        @test length(s1) == 2

        s2 = 5 + a
        @test s2 isa QAdd

        s3 = (a + ad) + 3
        @test s3 isa QAdd
        @test length(s3) == 3

        s4 = 3 + (a + ad)
        @test s4 isa QAdd
    end

    @testset "Subtraction" begin
        s = a - ad
        @test s isa QAdd
        @test length(s) == 2

        s2 = a - 3
        @test s2 isa QAdd
    end

    @testset "QAdd * Number" begin
        s = (a + ad) * 3
        @test s isa QAdd
        @test all(t -> prefactor(t) == 3, sorted_arguments(s))
    end

    @testset "QAdd * QSym (distributes)" begin
        s = (a + ad) * a
        @test s isa QAdd
        @test length(s) == 2
        @test all(t -> length(operators(t)) == 2, sorted_arguments(s))
    end

    @testset "QAdd * QAdd (distributes)" begin
        s = (a + ad) * (2 * a)
        @test s isa QAdd
        @test length(s) == 2
    end

    @testset "QAdd * QAdd (distributes)" begin
        s = (a + ad) * (a + ad)
        @test s isa QAdd
        @test length(s) == 4
    end

    @testset "Adjoint" begin
        s = a + 2 * ad
        sd = s'
        @test sd isa QAdd
        @test length(sd) == 2
    end

    @testset "Type stability" begin
        # QSym + QSym
        @inferred a + ad
        # QAdd + QAdd (single-term)
        @inferred (2 * a) + (3 * ad)
        # QAdd + QSym / QSym + QAdd
        @inferred (2 * a) + ad
        @inferred ad + (2 * a)
        # QAdd + QAdd (merge)
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

    @testset "sum and reduce" begin
        terms = [a + ad, 2 * a, 3 * ad, a * ad, ad * a]
        manual = terms[1] + terms[2] + terms[3] + terms[4] + terms[5]

        # sum matches manual chain of +
        @test sum(terms) == manual

        # reduce(+, ...) matches
        @test reduce(+, terms) == manual

        # Inputs are not mutated
        args_before = copy(terms[1].arguments)
        sum(terms)
        @test terms[1].arguments == args_before

        # Single element
        @test sum([a + ad]) == a + ad

        # Empty
        @test sum(QAdd[]) == zero(QAdd)

        # sum with a function
        ops = [a, ad]
        @test sum(x -> x * a, ops) == a * a + ad * a

        # zero
        @test iszero(zero(QAdd))
        @test zero(a + ad) == zero(QAdd)
    end

    @testset "Allocations" begin
        s = a + ad
        m_2a = 2 * a
        m_3ad = 3 * ad

        # Warmup
        a + ad; m_2a + m_3ad; s + m_2a; s + s; s * 3; s * a; s * s

        # QSym + QSym
        @test @allocations(a + ad) <= 200
        # QAdd + QAdd (single-term)
        @test @allocations(m_2a + m_3ad) <= 200
        # QAdd + QMul
        @test @allocations(s + m_2a) <= 200
        # QAdd + QAdd
        @test @allocations(s + s) <= 200
        # QAdd * Number
        @test @allocations(s * 3) <= 500
        # QAdd * QSym (distributive)
        @test @allocations(s * a) <= 200
        # QAdd * QAdd (distributive)
        @test @allocations(s * s) <= 10000
        # Negation
        @test @allocations(-s) <= 500
        # Adjoint
        @test @allocations(adjoint(s)) <= 200
    end
end
