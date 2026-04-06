using SecondQuantizedAlgebra
using QuantumOpticsBase
using Latexify
using Symbolics: @variables, Num
using Test
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym, HilbertSpace

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
        @test m isa QMul
        @test length(m.args_nc) == 2

        s = x + p
        @test s isa QAdd
    end

    @testset "Simplify: [X, P] = i" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        # p*x out of order → x*p - i
        result = simplify(p * x)
        @test result isa QAdd
        @test length(result) == 2  # x·p + (-i)

        # x*p already ordered → x·p (1 term)
        result2 = simplify(x * p)
        @test result2 isa QAdd
        @test length(result2) == 1
    end

    @testset "Simplify: mixed spaces don't interact" begin
        h = PhaseSpace(:q1) ⊗ PhaseSpace(:q2)
        x1 = Position(h, :x, 1)
        p2 = Momentum(h, :p, 2)

        # Different spaces — no commutation
        result = simplify(p2 * x1)
        @test length(result) == 1
    end

    @testset "Higher-order algebra" begin
        h = PhaseSpace(:motion)
        x = Position(h, :x)
        p = Momentum(h, :p)

        @test isequal(
            simplify(normal_order(p * x * x)), simplify(x * x * p - 2im * x)
        )
    end

    @testset "Multi-space algebra" begin
        hps1 = PhaseSpace(:motion1)
        hps2 = PhaseSpace(:motion2)
        h = hps1 ⊗ hps2

        x1 = Position(h, :x_1, 1)
        p1 = Momentum(h, :p_1, 1)
        x2 = Position(h, :x_2, 2)

        @test isequal(p1 * x2, x2 * p1)
        @test isequal(
            simplify(normal_order(p1 * x2 * x1)),
            simplify((x1 * p1 - 1im) * x2),
        )
    end

    @testset "Hamiltonian construction and commutators" begin
        hps1 = PhaseSpace(:motion1)
        hps2 = PhaseSpace(:motion2)
        h = hps1 ⊗ hps2

        x1 = Position(h, :x_1, 1)
        p1 = Momentum(h, :p_1, 1)
        x2 = Position(h, :x_2, 2)
        p2 = Momentum(h, :p_2, 2)

        @variables ω::Real m::Real miv::Real c::Real
        H = 0.5 * Num(miv) * p1^2 + 0.5Num(m) * Num(ω)^2 * x1^2 + Num(c) * x1^4

        @test iszero(simplify(commutator(H, x2)))
        @test iszero(simplify(commutator(H, p2)))
        @test isequal(simplify(1im * commutator(H, x1)), simplify(p1 * Num(miv)))
        @test isequal(
            simplify(1im * commutator(H, p1)),
            simplify(-x1 * Num(m) * Num(ω)^2 - 4 * Num(c) * x1^3),
        )
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
        @test comm.data[1, 1] ≈ im atol = 1.0e-10
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

    @testset "Type stability" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)

        @inferred Position(:x, 1)
        @inferred Momentum(:p, 1)
        @inferred adjoint(x)
        @inferred adjoint(p)
        @inferred isequal(x, x)
        @inferred isequal(p, p)
        @inferred hash(x, UInt(0))
        @inferred hash(p, UInt(0))
        @inferred x * p
        @inferred x + p
    end

    @testset "Allocations" begin
        h = PhaseSpace(:q)
        x = Position(h, :x)
        p = Momentum(h, :p)
        x2 = Position(h, :x)
        p2 = Momentum(h, :p)

        @test @allocations(Position(:x, 1)) == 0
        @test @allocations(Momentum(:p, 1)) == 0
        @test @allocations(adjoint(x)) == 0
        @test @allocations(adjoint(p)) == 0
        @test @allocations(isequal(x, x2)) == 0
        @test @allocations(isequal(p, p2)) == 0
        @test @allocations(hash(x, UInt(0))) == 0
        @test @allocations(hash(p, UInt(0))) == 0
    end
end
