using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, HilbertSpace, QTermDict, _to_cnum, Index
using Test

@testset "spin" begin
    @testset "SpinSpace" begin
        h = SpinSpace(:s, 1 // 2)
        @test h isa HilbertSpace
        @test h.name == :s
        @test h.spin == 1 // 2
        @test SpinSpace(:s, 1 // 2) == SpinSpace(:s, 1 // 2)
        @test SpinSpace(:s, 1 // 2) != SpinSpace(:s, 1)
    end

    @testset "Spin construction — single space" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        @test Sx isa Spin
        @test Sx isa QSym
        @test Sx.axis == 1
        @test Sy.axis == 2
        @test Sz.axis == 3
        @test Sx.space_index == 1
        @test_throws ArgumentError Spin(h, :S, 4)
    end

    @testset "Spin construction — product space" begin
        h = FockSpace(:c) ⊗ SpinSpace(:s, 1)
        Sx = Spin(h, :S, 1, 2)
        @test Sx.space_index == 2
        @test_throws ArgumentError Spin(h, :S, 1, 1)
    end

    @testset "Adjoint — Hermitian" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        @test Sx' == Sx
    end

    @testset "Equality and hashing" begin
        h = SpinSpace(:s, 1 // 2)
        S1 = Spin(h, :S, 1)
        S2 = Spin(h, :S, 1)
        S3 = Spin(h, :S, 2)
        @test isequal(S1, S2)
        @test !isequal(S1, S3)
        @test hash(S1) == hash(S2)
        @test hash(S1) != hash(S3)
    end

    @testset "Arithmetic" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        m = Sx * Sy
        @test m isa QAdd

        s = Sx + Sy
        @test s isa QAdd
    end

    @testset "@qnumbers" begin
        h = SpinSpace(:s, 1 // 2)
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
        @test isequal(simplify(normal_order(s(3) * s(3))), QAdd(QTermDict(QSym[] => _to_cnum(1)), Index[], Tuple{Index, Index}[]))
        @test isequal(simplify(normal_order(s(3) * s(1))), simplify(1im * s(2)))
        @test isequal(
            simplify(normal_order(s(1) * s(2) * s(3))), QAdd(QTermDict(QSym[] => _to_cnum(1im)), Index[], Tuple{Index, Index}[])
        )

        # adjoint (Hermitian)
        @test adjoint(s(1)) == s(1)
        @test adjoint(s(2)) == s(2)
        @test adjoint(s(3)) == s(3)
    end

    @testset "Multi-Pauli space operations" begin
        hs1 = PauliSpace(:Spin1)
        hs2 = PauliSpace(:Spin2)
        h = hs1 ⊗ hs2

        σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)
        σx(i) = σ(i, 1)
        σy(i) = σ(i, 2)
        σz(i) = σ(i, 3)

        # Different spins commute
        @test isequal(σx(2) * σx(1), σx(1) * σx(2))
        @test isequal(σy(2) * σx(1), σx(1) * σy(2))
        @test isequal(σz(2) * σz(1), σz(1) * σz(2))
    end

    @testset "Collective spin commutation — [Sx, Sy] = iSz" begin
        hcs = SpinSpace(:Spin1, 1 // 2)
        S(axis) = Spin(hcs, :S, axis)

        @test isequal(simplify(commutator(S(1), S(2))), simplify(1im * S(3)))
        @test isequal(simplify(commutator(S(1), S(3))), simplify(-1im * S(2)))
        @test isequal(simplify(commutator(S(2), S(3))), simplify(1im * S(1)))
    end

    @testset "Spin operator properties" begin
        hcs = SpinSpace(:Spin1, 1 // 2)
        S(axis) = Spin(hcs, :S, axis)

        @test isequal(S(1) * S(2), 1 * S(1) * S(2))
        @test !isequal(S(1) * S(2), S(1) * S(3))
        @test isequal(S(1) * S(1), (S(1))^2)
        @test isequal(average(S(1) + S(2)), average(S(2) + S(1)))
        @test isequal(simplify(S(2) + S(2)), simplify(2S(2)))
        @test isequal(2 * S(1) * S(2) * S(3), S(1) * S(2) * S(3) * 2)
    end

    @testset "Multi-spin operations" begin
        hcs1 = SpinSpace(:Spin1, 1 // 2)
        hcs2 = SpinSpace(:Spin2, 1 // 2)
        h = hcs1 ⊗ hcs2

        S(i, axis) = Spin(h, Symbol(:S_, i), axis, i)
        Sx(i) = S(i, 1)
        Sy(i) = S(i, 2)
        Sz(i) = S(i, 3)

        # Different spins commute
        @test isequal(Sx(2) * Sx(1), Sx(1) * Sx(2))
        @test isequal(Sy(2) * Sx(1), Sx(1) * Sy(2))
        @test isequal(Sz(2) * Sz(1), Sz(1) * Sz(2))
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(SpinSpace)
        all_concrete(Spin)
    end

    @testset "Type stability" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)

        @inferred Spin(:S, 1, 1)
        @inferred adjoint(Sx)
        @inferred isequal(Sx, Sy)
        @inferred hash(Sx, UInt(0))
        @inferred Sx * Sy
        @inferred Sx + Sy
    end

    @testset "Allocations" begin
        h = SpinSpace(:s, 1 // 2)
        Sx = Spin(h, :S, 1)
        Sx2 = Spin(h, :S, 1)

        @test @allocations(Spin(:S, 1, 1)) == 0
        @test @allocations(adjoint(Sx)) == 0
        @test @allocations(isequal(Sx, Sx2)) == 0
        @test @allocations(hash(Sx, UInt(0))) == 0
    end
end
