using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: simplify, QAdd, QSym, _single_qadd, _CNUM_ONE, _to_cnum
using Test

@testset "fock operators" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Construction — product space" begin
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        @test a1.space_index == 1
        @test b1.space_index == 2
        @test_throws ArgumentError Destroy(hp, :c, 3)
    end

    @testset "Adjoint" begin
        @test is_create(ad)
        @test operator_name(ad) == :a
        @test ad.space_index == 1
        @test ad' == a
        @test a'' == a
    end

    @testset "Canonical ordering helpers" begin
        import SecondQuantizedAlgebra: ladder
        @test ladder(a) == 1
        @test ladder(ad) == 0
    end

    @testset "Algebraic simplification" begin
        @test isequal(simplify(a + a), simplify(2 * a))
        @test isequal(simplify(a / 2 + 0.5 * a), simplify(1.0 * a))
    end

    @testset "Commutation relations" begin
        # Eager ordering: a * a' is already normal-ordered
        @test isequal(a * ad, 1 + ad * a)
        @test isequal(a * ad + 1, 2 + ad * a)
    end

    @testset "Hamiltonian evolution" begin
        ωc = 0.1313
        H = ωc * ad * a
        da = simplify(1.0im * commutator(H, a))
        @test isequal(da, simplify(-1.0im * ωc * a))
    end

    @testset "Multi-space commutators" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        h3 = FockSpace(:c3)
        h4 = FockSpace(:c4)
        hm = h1 ⊗ h2 ⊗ h3 ⊗ h4

        a1 = Destroy(hm, :a1, 1)
        a2 = Destroy(hm, :a2, 2)

        @test isequal(
            simplify(commutator(a1 + a2, a1')), _single_qadd(_to_cnum(1), Op[])
        )
        @test isequal(
            simplify(commutator(a2', a1 + a2)), _single_qadd(_to_cnum(-1), Op[])
        )

        @test commutator(a1, 1) == commutator(1, a2)
    end

    @testset "Type stability" begin
        import SecondQuantizedAlgebra: ladder
        @inferred Destroy(:a, 1)
        @inferred Create(:a, 1)
        @inferred adjoint(a)
        @inferred adjoint(ad)
        @inferred isequal(a, a)
        @inferred hash(a, UInt(0))
        @inferred hash(ad, UInt(0))
        @inferred ladder(a)
        @inferred ladder(ad)
        @inferred a * ad
        @inferred ad * a
        @inferred a + ad
        @inferred a - ad
        @inferred a^2
        @inferred 2 * a
        @inferred commutator(a, ad)
    end

    @static if VERSION >= v"1.12"
        @testset "Allocations" begin
            a2 = Destroy(h, :a)
            @test @allocations(Destroy(:a, 1)) == 0
            @test @allocations(Create(:a, 1)) == 0
            @test @allocations(adjoint(a)) == 0
            @test @allocations(adjoint(ad)) == 0
            @test @allocations(isequal(a, a2)) == 0
            @test @allocations(hash(a, UInt(0))) == 0
        end
    end
end
