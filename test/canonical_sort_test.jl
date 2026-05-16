using Test
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: _partial_sort!, _EMPTY_NE, Index, NO_INDEX, QSym

@testset "Partial sort: distinct sites" begin
    h = FockSpace(:c) ⊗ NLevelSpace(:a, 2)
    a = Destroy(h, :a, 1)
    σ = Transition(h, :σ, 1, 2)
    ops = QSym[σ, a]
    _partial_sort!(ops, _EMPTY_NE)
    # Destroy < Transition by type order
    @test ops[1] isa Destroy
    @test ops[2] isa Transition
end

@testset "Partial sort: same-site preserved" begin
    h = NLevelSpace(:a, 3)
    σ12 = Transition(h, :σ, 1, 2)
    σ21 = Transition(h, :σ, 2, 1)
    ops = QSym[σ12, σ21]
    pre = copy(ops)
    _partial_sort!(ops, _EMPTY_NE)
    @test ops == pre   # same site never reordered
end

@testset "Partial sort: undetermined preserved" begin
    ha = NLevelSpace(:a, 2)
    i = Index(ha, :i, 5, ha)
    j = Index(ha, :j, 5, ha)
    # Use IndexedOperator to attach symbolic indices to operators
    σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
    σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
    ops = QSym[σ_j, σ_i]   # physical order: j first
    _partial_sort!(ops, _EMPTY_NE)
    @test ops[1] === σ_j   # preserved because (i, j) relationship is Undetermined
    @test ops[2] === σ_i
end

@testset "Partial sort: ne resolves undetermined to distinct" begin
    ha = NLevelSpace(:a, 2)
    i = Index(ha, :i, 5, ha)
    j = Index(ha, :j, 5, ha)
    σ_i = IndexedOperator(Transition(ha, :σ, 1, 2), i)
    σ_j = IndexedOperator(Transition(ha, :σ, 1, 2), j)
    ops = QSym[σ_j, σ_i]
    _partial_sort!(ops, [(i, j)])
    # Now distinct: sort by index, i < j alphabetically
    @test ops[1] === σ_i
    @test ops[2] === σ_j
end
