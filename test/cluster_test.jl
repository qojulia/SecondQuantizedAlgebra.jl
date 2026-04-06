using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: _unwrap_space, has_cluster, simplify, QAdd, QSym, sorted_arguments

@testset "ClusterSpace" begin
    @testset "Construction" begin
        ha = NLevelSpace(:atom, 2, 1)
        c = ClusterSpace(ha, 10, 2)
        @test c.original_space == ha
        @test c.N == 10
        @test c.order == 2
        @test _unwrap_space(c) === ha
        @test _unwrap_space(ha) === ha
    end

    @testset "Equality and hashing" begin
        ha = NLevelSpace(:atom, 2, 1)
        c1 = ClusterSpace(ha, 10, 2)
        c2 = ClusterSpace(ha, 10, 2)
        c3 = ClusterSpace(ha, 20, 2)
        @test c1 == c2
        @test c1 != c3
        @test hash(c1) == hash(c2)
    end

    @testset "has_cluster" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        @test has_cluster(c) == true
        @test has_cluster(ha) == false
        @test has_cluster(hf) == false

        h = hf ⊗ c
        @test has_cluster(h) == true
        h2 = hf ⊗ ha
        @test has_cluster(h2) == false
    end

    @testset "ProductSpace with ClusterSpace" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        # Operators on FockSpace position work normally
        a = Destroy(h, :a, 1)
        @test a isa Destroy
        @test a.space_index == 1
        @test a.copy_index == 1

        # Operators on ClusterSpace position validate against original_space
        σ = Transition(h, :σ, 1, 2, 2)
        @test σ isa Transition
        @test σ.space_index == 2
        @test σ.copy_index == 1
    end

    @testset "cluster_expand" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        σ = Transition(h, :σ, 1, 2, 2)
        copies = cluster_expand(σ, 2)
        @test length(copies) == 2
        @test copies[1].name == :σ_1
        @test copies[1].copy_index == 1
        @test copies[1].space_index == 2
        @test copies[2].name == :σ_2
        @test copies[2].copy_index == 2
        @test copies[2].space_index == 2

        # Expand via ProductSpace
        copies2 = cluster_expand(σ, h)
        @test length(copies2) == 2
        @test copies2[1].name == :σ_1
        @test copies2[2].name == :σ_2

        # Fock operator expand
        a = Destroy(h, :a, 1)
        @test_throws ArgumentError cluster_expand(a, h)  # FockSpace, not ClusterSpace
        fock_copies = cluster_expand(a, 3)
        @test length(fock_copies) == 3
        @test fock_copies[1].name == :a_1
        @test fock_copies[2].copy_index == 2
    end

    @testset "copy_index field on all QSym types" begin
        # Default copy_index = 1
        @test Destroy(:a, 1).copy_index == 1
        @test Create(:a, 1).copy_index == 1
        @test Transition(:σ, 1, 2, 1).copy_index == 1
        @test Pauli(:σ, 1, 1).copy_index == 1
        @test Spin(:S, 1, 1).copy_index == 1
        @test Position(:x, 1).copy_index == 1
        @test Momentum(:p, 1).copy_index == 1

        # Explicit copy_index
        @test Destroy(:a, 1, 2).copy_index == 2
        @test Create(:a, 1, 3).copy_index == 3
        @test Transition(:σ, 1, 2, 1, 2).copy_index == 2
    end

    @testset "Adjoint preserves copy_index" begin
        a = Destroy(:a, 1, 2)
        @test a'.copy_index == 2
        @test a'.name == :a

        σ = Transition(:σ, 1, 2, 1, 3)
        @test σ'.copy_index == 3
        @test σ'.i == 2 && σ'.j == 1
    end

    @testset "Equality includes copy_index" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        @test a1 != a2
        @test !isequal(a1, a2)
        @test hash(a1) != hash(a2)
    end

    @testset "Canonical ordering uses (space_index, copy_index)" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        m = a2 * a1  # should sort by copy_index
        ops = operators(only(sorted_arguments(m)))
        @test ops[1].copy_index == 1
        @test ops[2].copy_index == 2
    end

    @testset "Different copy_index operators commute" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)

        # Different copy → short-circuit to zero
        @test iszero(commutator(a1, a2'))

        # Same copy → non-trivial
        result = commutator(a1, a1')
        @test result isa QAdd
        @test !iszero(result)
    end

    @testset "Simplify respects copy_index" begin
        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)

        # Same copy: a·a† → a†·a + 1
        result = simplify(a1 * a1')
        @test length(result) == 2

        # Different copy: a₁·a₂† stays as-is (no commutation rule)
        result = simplify(a1 * a2')
        @test length(result) == 1
        @test length(operators(only(sorted_arguments(result)))) == 2
    end

    @testset "Symbolic N" begin
        @variables N_sym
        ha = NLevelSpace(:atom, 2, 1)
        c = ClusterSpace(ha, N_sym, 2)
        @test c.N === N_sym
        @test c.order == 2
    end

    @testset "Type stability" begin
        ha = NLevelSpace(:atom, 2, 1)
        hf = FockSpace(:c)
        c = ClusterSpace(ha, 10, 2)
        h = hf ⊗ c

        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)

        @inferred commutator(a, σ)
        @inferred simplify(a * a')

        a1 = Destroy(:a, 1, 1)
        a2 = Destroy(:a, 1, 2)
        @inferred commutator(a1, a2')
        @inferred commutator(a1, a1')
        @inferred simplify(a1 * a2')
    end
end
