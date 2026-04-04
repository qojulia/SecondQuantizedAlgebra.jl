using SecondQuantizedAlgebra
using QuantumOpticsBase
using Latexify
using Test
import SecondQuantizedAlgebra: simplify

@testset "phase_space" begin
    @testset "PhaseSpace" begin
        h = PhaseSpace(:q)
        @test h isa HilbertSpace
        @test h.name == :q
        @test PhaseSpace(:q) == PhaseSpace(:q)
        @test PhaseSpace(:q) != PhaseSpace(:r)
        @test hash(PhaseSpace(:q)) == hash(PhaseSpace(:q))
    end

    @testset "Position construction — single space" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        @test x isa Position
        @test x isa QSym
        @test x.name == :x
        @test x.space_index == 1
    end

    @testset "Momentum construction — single space" begin
        h = PhaseSpace(:q)
        p = Momentum(h, :p)
        @test p isa Momentum
        @test p isa QSym
        @test p.name == :p
        @test p.space_index == 1
    end

    @testset "Construction — product space" begin
        h = FockSpace(:c) ⊗ PhaseSpace(:q)
        x = Position(h, :x, 2)
        p = Momentum(h, :p, 2)
        @test x.space_index == 2
        @test p.space_index == 2
        @test_throws ArgumentError Position(h, :x, 1)  # FockSpace, not PhaseSpace
        @test_throws ArgumentError Momentum(h, :p, 1)
        @test_throws ArgumentError Position(h, :x, 3)  # out of range
    end

    @testset "Adjoint — Hermitian" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        @test x' == x
        @test p' == p
    end

    @testset "Equality and hashing" begin
        h = PhaseSpace(:q)
        x1 = Position(h, :x)
        x2 = Position(h, :x)
        p1 = Momentum(h, :p)
        @test isequal(x1, x2)
        @test !isequal(x1, Position(h, :y))
        @test hash(x1) == hash(x2)
        @test hash(x1) != hash(p1)
    end

    @testset "Lazy multiplication" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        m = x * p
        @test m isa QMul{Int}
        @test length(m.args_nc) == 2

        s = x + p
        @test s isa QAdd{Int}
    end

    @testset "Simplify: [X, P] = i" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        # p*x out of order → x*p - i
        result = simplify(p * x)
        @test result isa QAdd
        @test length(result.arguments) == 2  # x·p + (-i)

        # x*p already ordered → x·p (1 term)
        result2 = simplify(x * p)
        @test result2 isa QAdd
        @test length(result2.arguments) == 1
    end

    @testset "Simplify: mixed spaces don't interact" begin
        h = PhaseSpace(:q1) ⊗ PhaseSpace(:q2)
        x1 = Position(h, :x, 1)
        p2 = Momentum(h, :p, 2)

        # Different spaces — no commutation
        result = simplify(p2 * x1)
        @test length(result.arguments) == 1
    end

    @testset "Numeric conversion" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        b = FockBasis(10)

        x_num = to_numeric(x, b)
        p_num = to_numeric(p, b)
        @test x_num ≈ (destroy(b) + create(b)) / sqrt(2)
        @test p_num ≈ im * (create(b) - destroy(b)) / sqrt(2)
    end

    @testset "Numeric: commutator [X, P] = i" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        b = FockBasis(30)

        x_num = to_numeric(x, b)
        p_num = to_numeric(p, b)
        comm = x_num * p_num - p_num * x_num
        # Should be approximately i*I (exact in infinite basis, approximate in truncated)
        @test comm.data[1, 1] ≈ im atol = 1e-10
    end

    @testset "Numeric: composite basis" begin
        h = FockSpace(:c) ⊗ PhaseSpace(:q)
        x = Position(h, :x, 2)
        bf = FockBasis(3)
        bq = FockBasis(5)
        bc = bf ⊗ bq
        @test to_numeric(x, bc) isa LazyTensor
    end

    @testset "Printing" begin
        h = PhaseSpace(:q)
        @test sprint(show, h) == "ℋ(q)"

        x = Position(h, :x)
        p = Momentum(h, :p)
        @test sprint(show, x) == "x"
        @test sprint(show, p) == "p"
    end

    @testset "LaTeX" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        sx = latexify(x)
        sp = latexify(p)
        @test occursin("hat", sx)
        @test occursin("x", sx)
        @test occursin("hat", sp)
        @test occursin("p", sp)

        # MIME dispatch
        @test repr(MIME"text/latex"(), x) == sx
        @test repr(MIME"text/latex"(), p) == sp
    end

    @testset "@qnumbers" begin
        h = PhaseSpace(:q)
        @qnumbers x::Position(h)
        @qnumbers p::Momentum(h)
        @test x isa Position
        @test x.name == :x
        @test p isa Momentum
        @test p.name == :p
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(PhaseSpace)
        all_concrete(Position)
        all_concrete(Momentum)
    end
end
