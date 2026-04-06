using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: QMul, QAdd, QSym, CNum, NO_INDEX

@testset "insert_index" begin
    hf = FockSpace(:f)
    hn = NLevelSpace(:n, 2, 1)
    h = hf ⊗ hn

    a = Destroy(h, :a, 1)
    i = Index(h, :i, 10, hn)
    j = Index(h, :j, 10, hn)

    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

    @testset "Operators" begin
        # Inserts concrete copy_index, clears symbolic index
        σ1 = insert_index(σ(1, 2, i), i, 1)
        @test σ1 isa Transition
        @test σ1.copy_index == 1
        @test !has_index(σ1.index)
        @test σ1.space_index == 2

        # Equivalent to direct numbered construction
        @test isequal(σ1, IndexedOperator(Transition(h, :σ, 1, 2, 2), 1))

        # No-op when index doesn't match
        σ_noop = insert_index(σ(1, 2, i), j, 2)
        @test has_index(σ_noop.index)
        @test isequal(σ_noop, σ(1, 2, i))

        # Destroy/Create
        ai = IndexedOperator(a, i)
        a1 = insert_index(ai, i, 3)
        @test a1.copy_index == 3
        @test !has_index(a1.index)
    end

    @testset "Number passthrough" begin
        @test insert_index(1, i, 1) == 1
        @test insert_index(0, i, 5) == 0
        @test insert_index(3.14, i, 1) == 3.14
    end

    @testset "IndexedVariable" begin
        gi = IndexedVariable(:g, i)
        g1 = insert_index(gi, i, 1)
        @test g1 isa Symbolics.Num

        # No-op when index doesn't match
        g_noop = insert_index(gi, j, 2)
        @test isequal(g_noop, gi)
    end

    @testset "DoubleIndexedVariable" begin
        Γij = DoubleIndexedVariable(:Γ, i, j)
        # Insert one index
        Γ1j = insert_index(Γij, i, 1)
        @test Γ1j isa Symbolics.Num
        # Insert both indices
        Γ12 = insert_index(Γ1j, j, 2)
        @test Γ12 isa Symbolics.Num

        # identical=false → 0 when both become equal
        Ωij = DoubleIndexedVariable(:Ω, i, j; identical = false)
        Ω11 = insert_index(insert_index(Ωij, i, 1), j, 1)
        @test isequal(Ω11, Symbolics.Num(0))
    end

    @testset "QMul" begin
        gi = IndexedVariable(:g, i)
        m = gi * σ(1, 2, i)
        m1 = insert_index(m, i, 1)
        @test m1 isa QMul
        @test m1.args_nc[1].copy_index == 1
        @test !has_index(m1.args_nc[1].index)
    end

    @testset "QAdd" begin
        s = σ(1, 2, i) + σ(2, 1, j)
        s1 = insert_index(s, i, 1)
        # First term should be numbered, second unchanged
        @test s1 isa QAdd
    end

    @testset "Complex expressions" begin
        # Inserting into products of indexed operators
        σ1_numbered = IndexedOperator(Transition(h, :σ, 1, 2, 2), 2)
        expr = σ(1, 2, i) * σ1_numbered
        result = insert_index(expr, i, 1)
        @test result isa QMul
        @test all(!has_index(op.index) for op in result.args_nc)
    end
end
