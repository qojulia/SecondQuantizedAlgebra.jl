using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QTermDict, _CNUM_ONE, _to_cnum
using Test

@testset "fock operators" begin
    @testset "Construction — single space" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        @test a isa Destroy
        @test a isa QSym
        @test a.name == :a
        @test a.space_index == 1
    end

    @testset "Construction — product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        a = Destroy(h, :a, 1)
        b = Destroy(h, :b, 2)
        @test a.space_index == 1
        @test b.space_index == 2
        @test_throws ArgumentError Destroy(h, :c, 3)
    end

    @testset "Adjoint" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'
        @test ad isa Create
        @test ad.name == :a
        @test ad.space_index == 1
        @test ad' == a
        @test a'' == a
    end

    @testset "Equality and hashing" begin
        h = FockSpace(:c)
        a1 = Destroy(h, :a)
        a2 = Destroy(h, :a)
        b = Destroy(h, :b)
        @test isequal(a1, a2)
        @test !isequal(a1, b)
        @test hash(a1) == hash(a2)
        @test hash(a1) != hash(b)
        @test hash(a1) != hash(a1')
    end

    @testset "Canonical ordering helpers" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        import SecondQuantizedAlgebra: ladder
        @test ladder(a) == 1
        @test ladder(a') == 0
    end

    @testset "Algebraic simplification" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)

        @test isequal(simplify(a + a), simplify(2 * a))
        @test isequal(simplify(a / 2 + 0.5 * a), simplify(1.0 * a))
    end

    @testset "Commutation relations" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)

        # Eager ordering: a * a' is already normal-ordered
        @test isequal(a * a', 1 + a' * a)
        @test isequal(a * a' + 1, 2 + a' * a)
    end

    @testset "Hamiltonian evolution" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)

        ωc = 0.1313
        H = ωc * a' * a
        da = simplify(1.0im * commutator(H, a))
        @test isequal(da, simplify(-1.0im * ωc * a))
    end

    @testset "Multi-space commutators" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        h3 = FockSpace(:c3)
        h4 = FockSpace(:c4)
        h = h1 ⊗ h2 ⊗ h3 ⊗ h4

        a1 = Destroy(h, :a1, 1)
        a2 = Destroy(h, :a2, 2)

        @test isequal(
            simplify(commutator(a1 + a2, a1')), QAdd(QTermDict(QSym[] => _to_cnum(1)), Index[], Tuple{Index, Index}[])
        )
        @test isequal(
            simplify(commutator(a2', a1 + a2)), QAdd(QTermDict(QSym[] => _to_cnum(-1)), Index[], Tuple{Index, Index}[])
        )

        @test commutator(a1, 1) == commutator(1, a2)
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(Destroy)
        all_concrete(Create)
    end

    @testset "Type stability" begin
        import SecondQuantizedAlgebra: ladder
        h = FockSpace(:c)
        a = Destroy(h, :a)
        @inferred Destroy(:a, 1)
        @inferred Create(:a, 1)
        @inferred adjoint(a)
        @inferred adjoint(a')
        @inferred isequal(a, a)
        @inferred hash(a, UInt(0))
        @inferred hash(a', UInt(0))
        @inferred ladder(a)
        @inferred ladder(a')
    end

    @testset "Allocations" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)

        # Construction
        @test @allocations(Destroy(:a, 1)) == 0
        @test @allocations(Create(:a, 1)) == 0

        # Adjoint
        @test @allocations(adjoint(a)) == 0
        @test @allocations(adjoint(a')) == 0

        # Equality / hashing
        a2 = Destroy(h, :a)
        @test @allocations(isequal(a, a2)) == 0
        @test @allocations(hash(a, UInt(0))) == 0
    end
end
