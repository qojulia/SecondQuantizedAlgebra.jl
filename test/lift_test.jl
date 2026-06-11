using SecondQuantizedAlgebra
using Test
using SymbolicUtils: SymbolicUtils, symtype
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: AverageOperator, make_time_dependent

@testset "make_time_dependent lift" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables t
    leaf = average(a)
    lifted = make_time_dependent(leaf, t)
    @test SymbolicUtils.iscall(lifted)
    @test isequal(SymbolicUtils.arguments(lifted)[1], SymbolicUtils.unwrap(t))
    @test symtype(lifted) === Number
    @test is_average(lifted)
    @test isequal(undo_average(lifted), undo_average(leaf))   # operator recovered from metadata
    @test isequal(lifted, make_time_dependent(average(a), t))               # structural identity
    @test !isequal(lifted, make_time_dependent(average(a'), t))             # distinct operators distinct

    # the lift rewrites ONLY average leaves, leaving the surrounding structure intact:
    # lifting a compound expression equals the compound of its individually-lifted leaves
    expr = 2 * average(a) + average(a' * a)
    lifted_expr = make_time_dependent(expr, t)
    @test isequal(
        lifted_expr,
        2 * make_time_dependent(average(a), t) + make_time_dependent(average(a' * a), t)
    )

    # distinct operators on DIFFERENT subspaces must stay distinct unknowns: the lifted
    # name is structural (space_index + operator), not a hash digest that could collide
    hm = FockSpace(:m1) ⊗ FockSpace(:m2)
    a1 = Destroy(hm, :a, 1); a2 = Destroy(hm, :a, 2)
    @test !isequal(make_time_dependent(average(a1), t), make_time_dependent(average(a2), t))

    # acts_on reads the lifted-form metadata
    h2 = FockSpace(:a) ⊗ NLevelSpace(:b, 2)
    aa = Destroy(h2, :a, 1)
    σ = Transition(h2, :σ, 1, 2, 2)
    @test acts_on(make_time_dependent(average(aa' * σ), t)) == [1, 2]

    # inner_adjoint of a lifted average stays lifted (time-dependent) and recovers ⟨op'⟩
    aladj = inner_adjoint(make_time_dependent(average(aa), t))
    @test is_average(aladj)
    @test SymbolicUtils.hasmetadata(aladj, AverageOperator)               # still lifted
    @test isequal(undo_average(aladj), undo_average(average(aa')))

    # inner_adjoint preserves the Σ scope of a summed lifted average
    hi = FockSpace(:site); ai = Destroy(hi, :a); i = Index(hi, :i, 3, hi)
    summed = make_time_dependent(average(Σ(IndexedOperator(ai', i) * IndexedOperator(ai, i), i)), t)
    @test SecondQuantizedAlgebra.has_sum_metadata(inner_adjoint(summed))
end
