using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: get_sum_body, get_sum_indices, has_sum_metadata

@testset "Indexed time-dependent averages" begin
    h = FockSpace(:site)
    a = Destroy(h, :a)
    i = Index(h, :i, 3, h)
    @variables t

    summed = Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i)
    lifted = make_time_dependent(average(summed), t)

    @test is_indexed_sum(lifted)
    @test has_sum_metadata(lifted)
    @test get_sum_indices(lifted) == [i]
    @test is_average(get_sum_body(lifted))
    @test get_indices(undo_average(lifted)) == [i]
    @test has_sum_metadata(inner_adjoint(lifted))
end
