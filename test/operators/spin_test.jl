using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, HilbertSpace, _single_qadd, _to_cnum, Index
using Test

@testset "spin" begin
    @testset "Spin construction — product space" begin
        h = FockSpace(:c) ⊗ SpinSpace(:s)
        Sx = Spin(h, :S, 1, 2)
        @test Sx.space_index == 2
        @test_throws ArgumentError Spin(h, :S, 1, 1)
    end

    @testset "Adjoint — Hermitian" begin
        h = SpinSpace(:s)
        Sx = Spin(h, :S, 1)
        @test Sx' == Sx
    end

    @testset "Arithmetic" begin
        h = SpinSpace(:s)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        m = Sx * Sy
        @test m isa QAdd

        s = Sx + Sy
        @test s isa QAdd
    end

    @testset "@qnumbers" begin
        h = SpinSpace(:s)
        @qnumbers Sx::Spin(h, 1)
        @test Sx isa Spin
        @test Sx.name == :Sx
        @test Sx.axis == 1
    end

    @testset "Pauli algebra — single space" begin
        hs = PauliSpace(:Spin1)
        s(axis) = Pauli(hs, :σ, axis)

        @test hash(s(1)) != hash(s(2))

        # Pauli multiplication rules (normal_order applies product rules)
        @test isequal(simplify(normal_order(s(1) * s(2))), simplify(1im * s(3)))
        @test !isequal(simplify(normal_order(s(1) * s(2))), simplify(1im * s(2)))
        @test isequal(simplify(normal_order(s(1) * s(3))), simplify(-1im * s(2)))
        @test isequal(simplify(normal_order(s(3) * s(3))), _single_qadd(_to_cnum(1), QSym[]))
        @test isequal(simplify(normal_order(s(3) * s(1))), simplify(1im * s(2)))
        @test isequal(
            simplify(normal_order(s(1) * s(2) * s(3))), _single_qadd(_to_cnum(1im), QSym[]),
        )

        # adjoint (Hermitian)
        for axis in 1:3
            @test adjoint(s(axis)) == s(axis)
        end
    end

    @testset "Multi-Pauli space operations" begin
        hs1 = PauliSpace(:Spin1)
        hs2 = PauliSpace(:Spin2)
        h = hs1 ⊗ hs2

        σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)

        # Different spins commute regardless of axis
        for (ax1, ax2) in ((1, 1), (1, 2), (3, 3))
            @test isequal(σ(2, ax2) * σ(1, ax1), σ(1, ax1) * σ(2, ax2))
        end
    end

    @testset "Collective spin commutation — [Sₐ, Sᵦ] = iε_{abc} S_c" begin
        hcs = SpinSpace(:Spin1)
        S(axis) = Spin(hcs, :S, axis)

        for (a, b, c, sign) in ((1, 2, 3, 1), (2, 3, 1, 1), (1, 3, 2, -1))
            @test isequal(simplify(commutator(S(a), S(b))), simplify(sign * 1im * S(c)))
        end
    end

    @testset "Spin operator properties" begin
        hcs = SpinSpace(:Spin1)
        S(axis) = Spin(hcs, :S, axis)

        @test isequal(S(1) * S(2), 1 * S(1) * S(2))
        @test !isequal(S(1) * S(2), S(1) * S(3))
        @test isequal(S(1) * S(1), (S(1))^2)
        @test isequal(average(S(1) + S(2)), average(S(2) + S(1)))
        @test isequal(simplify(S(2) + S(2)), simplify(2S(2)))
        @test isequal(2 * S(1) * S(2) * S(3), S(1) * S(2) * S(3) * 2)
    end

    @testset "Multi-spin operations" begin
        hcs1 = SpinSpace(:Spin1)
        hcs2 = SpinSpace(:Spin2)
        h = hcs1 ⊗ hcs2

        S(i, axis) = Spin(h, Symbol(:S_, i), axis, i)

        # Different spins commute regardless of axis
        for (ax1, ax2) in ((1, 1), (1, 2), (3, 3))
            @test isequal(S(2, ax2) * S(1, ax1), S(1, ax1) * S(2, ax2))
        end
    end

    @testset "Type stability" begin
        h = SpinSpace(:s)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        @inferred Spin(:S, 1, 1)
        @inferred adjoint(Sx)
        @inferred isequal(Sx, Sy)
        @inferred hash(Sx, UInt(0))
        @inferred Sx * Sy
        @inferred Sx + Sy
    end

    @static if VERSION >= v"1.12"
        @testset "Allocations" begin
            h = SpinSpace(:s)
            Sx = Spin(h, :S, 1)
            Sx2 = Spin(h, :S, 1)

            @test @allocations(Spin(:S, 1, 1)) == 0
            @test @allocations(adjoint(Sx)) == 0
            @test @allocations(isequal(Sx, Sx2)) == 0
            @test @allocations(hash(Sx, UInt(0))) == 0
        end
    end
end
