using SecondQuantizedAlgebra
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
