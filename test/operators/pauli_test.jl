using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, HilbertSpace
using Test

@testset "pauli" begin
    @testset "Pauli construction — product space" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p)
        σx = Pauli(h, :σ, 1, 2)
        @test σx.space_index == 2
        @test_throws ArgumentError Pauli(h, :σ, 1, 1)
    end

    @testset "Adjoint — Hermitian" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        @test σx' == σx
    end

    @testset "Arithmetic" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)

        m = σx * σy
        @test m isa QAdd

        s = σx + σy
        @test s isa QAdd
    end

    @testset "@qnumbers" begin
        h = PauliSpace(:p)
        @qnumbers σx::Pauli(h, 1)
        @test is_pauli(σx)
        @test operator_name(σx) == :σx
        @test σx.l1 == 1
    end

    @testset "Type stability" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)

        @inferred Pauli(:σ, 1, 1)
        @inferred adjoint(σx)
        @inferred isequal(σx, σy)
        @inferred hash(σx, UInt(0))
        @inferred σx * σy
        @inferred σx + σy
    end

    @static if VERSION >= v"1.12"
        @testset "Allocations" begin
            h = PauliSpace(:p)
            σx = Pauli(h, :σ, 1)
            σx2 = Pauli(h, :σ, 1)

            @test @allocations(Pauli(:σ, 1, 1)) == 0
            @test @allocations(adjoint(σx)) == 0
            @test @allocations(isequal(σx, σx2)) == 0
            @test @allocations(hash(σx, UInt(0))) == 0
        end
    end

    @testset "Algebra spot checks" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)

        # σ_x σ_y = i σ_z (and cyclic permutations).
        for (a, b, c) in ((σx, σy, σz), (σy, σz, σx), (σz, σx, σy))
            @test iszero(simplify(a * b - 1im * c))
        end

        # σ_x² = 1.
        @test iszero(simplify(σx * σx - 1))

        # [σ_x, σ_y] = 2i σ_z.
        @test iszero(simplify(commutator(σx, σy) - 2im * σz))
    end
end
