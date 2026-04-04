using SecondQuantizedAlgebra
using Test

@testset "QAdd" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QSym + QSym" begin
        s = a + ad
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
        @test all(x -> x isa QMul{Int}, s.arguments)
    end

    @testset "QMul + QMul" begin
        s = (2 * a) + (3 * ad)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
    end

    @testset "QMul + QSym" begin
        s1 = (2 * a) + ad
        @test s1 isa QAdd{Int}
        @test length(s1.arguments) == 2

        s2 = ad + (2 * a)
        @test s2 isa QAdd{Int}
    end

    @testset "QAdd + QMul" begin
        s = (a + ad) + (3 * a)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 3
    end

    @testset "QAdd + QAdd" begin
        s = (a + ad) + (a + ad)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 4
    end

    @testset "QField + Number" begin
        s1 = a + 5
        @test s1 isa QAdd{Int}
        @test length(s1.arguments) == 2

        s2 = 5 + a
        @test s2 isa QAdd{Int}

        s3 = (a + ad) + 3
        @test s3 isa QAdd{Int}
        @test length(s3.arguments) == 3

        s4 = 3 + (a + ad)
        @test s4 isa QAdd{Int}
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
        @test s isa QAdd{Int}
        @test all(x -> x.arg_c == 3, s.arguments)
    end

    @testset "QAdd * QSym (distributes)" begin
        s = (a + ad) * a
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
        @test all(x -> length(x.args_nc) == 2, s.arguments)
    end

    @testset "QAdd * QMul (distributes)" begin
        s = (a + ad) * (2 * a)
        @test s isa QAdd{Int}
        @test length(s.arguments) == 2
    end

    @testset "QAdd * QAdd (distributes)" begin
        s = (a + ad) * (a + ad)
        @test s isa QAdd{Int}
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
        @inferred a + ad
        @inferred (2 * a) + (3 * ad)
        @inferred (a + ad) + a
        @inferred (a + ad) * 3
        @inferred (a + ad) * a
        @inferred (a + ad) * (a + ad)
    end
end
