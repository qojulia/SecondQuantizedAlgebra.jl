using SecondQuantizedAlgebra
using Test
using SymbolicUtils: symtype
using Symbolics: @variables

@testset "Time-dependent averages" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables t

    @testset "lifting and recovery" begin
        leaf = average(a)
        lifted = make_time_dependent(leaf, t)

        @test is_average(lifted)
        @test symtype(lifted) === Number
        @test isequal(undo_average(lifted), undo_average(leaf))
        @test isequal(lifted, make_time_dependent(leaf, t))
        @test !isequal(lifted, make_time_dependent(average(a'), t))
        @test !isequal(make_time_dependent(average(a), t), make_time_dependent(average(a), 2t))
    end

    @testset "compound expressions preserve structure" begin
        compound = 2 * average(a) + average(a' * a)
        lifted = make_time_dependent(compound, t)
        expected = 2 * make_time_dependent(average(a), t) +
            make_time_dependent(average(a' * a), t)
        @test isequal(lifted, expected)
        @test !is_average(make_time_dependent(2a, t))
    end

    @testset "space identity and adjoints" begin
        hm = FockSpace(:left) ⊗ FockSpace(:right)
        a1 = Destroy(hm, :a, 1)
        a2 = Destroy(hm, :a, 2)
        @test !isequal(make_time_dependent(average(a1), t), make_time_dependent(average(a2), t))

        h2 = FockSpace(:mode) ⊗ NLevelSpace(:atom, 2)
        aa = Destroy(h2, :a, 1)
        σ = Transition(h2, :σ, 1, 2, 2)
        lifted = make_time_dependent(average(aa' * σ), t)
        @test acts_on(lifted) == [1, 2]

        aladj = inner_adjoint(lifted)
        @test is_average(aladj)
        @test isequal(undo_average(aladj), undo_average(average(σ' * aa)))
    end
end
