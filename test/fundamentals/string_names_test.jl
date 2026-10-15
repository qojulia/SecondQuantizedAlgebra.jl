using SecondQuantizedAlgebra
using Test

@testset "string names rejected" begin
    @testset "Hilbert spaces" begin
        @test_throws ArgumentError FockSpace("c")
        @test_throws ArgumentError PauliSpace("s")
        @test_throws ArgumentError SpinSpace("S")
        @test_throws ArgumentError PhaseSpace("o")
        @test_throws ArgumentError NLevelSpace("atom", 3)
        @test_throws ArgumentError NLevelSpace("atom", 3, 2)
        @test_throws ArgumentError NLevelSpace("atom", (:g, :e))
    end

    @testset "Operators — single space" begin
        h = FockSpace(:c)
        @test_throws ArgumentError Destroy(h, "a")
        @test_throws ArgumentError Create(h, "a")

        hn = NLevelSpace(:atom, 3)
        @test_throws ArgumentError Transition(hn, "s", 1, 2)

        @test_throws ArgumentError Pauli(PauliSpace(:s), "σ", 1)
        @test_throws ArgumentError Spin(SpinSpace(:S), "S", 1)

        hx = PhaseSpace(:o)
        @test_throws ArgumentError Position(hx, "x")
        @test_throws ArgumentError Momentum(hx, "p")
    end

    @testset "Operators — product space" begin
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        @test_throws ArgumentError Destroy(hp, "a", 1)
        @test_throws ArgumentError Create(hp, "b", 2)
    end

    @testset "Indices and indexed variables" begin
        h = FockSpace(:c)
        @test_throws ArgumentError Index(h, "i", 3, h)
        @test_throws ArgumentError Index(h, "i", 3, 1)
        i = Index(h, :i, 3, h)
        j = Index(h, :j, 3, h)
        @test_throws ArgumentError IndexedVariable("g", i)
        @test_throws ArgumentError DoubleIndexedVariable("J", i, j)
        @test_throws ArgumentError DoubleIndexedVariable("J", i, j; identical = false)
    end

    @testset "Message guides to Symbol" begin
        err = try
            FockSpace("cavity")
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Symbol", err.msg)
        @test occursin(":cavity", err.msg)
    end
end
