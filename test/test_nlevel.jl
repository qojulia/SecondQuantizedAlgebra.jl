using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "nlevel" begin

    @testset "Symbolic Level Systems" begin
        ha = NLevelSpace(:atom, (:g, :e))
        @test ha.GS == :g

        σ = Transition(ha, :σ, :g, :e)
        σee = Transition(ha, :σ, :e, :e)

        @testset "Basic Operator Properties" begin
            @test σ' == Transition(ha, :σ, :e, :g)
            @test σee' == σee
            @test σ'*σ == σee
        end

        @testset "Algebraic Relations" begin
            ex = σ*σ'
            @test isequal(simplify(σ*σ'), simplify(1 - σee))

            sz = σ'*σ - σ*σ'
            @test isequal(simplify(sz), simplify(2*σee - 1))
        end
    end

    @testset "Integer Level Systems" begin
        ha = NLevelSpace(:atom, 2)
        @test ha.GS == 1

        @test_throws AssertionError Transition(ha, :σ, :g, :e)

        σ = Transition(ha, :σ, 1, 2)
        σee = Transition(ha, :σ, :2, 2)
        σgg = Transition(ha, :σ, 1, 1)

        @testset "Basic Operator Properties" begin
            @test σ' == Transition(ha, :σ, 2, 1)
            @test σee' == σee
            @test σ'*σ == σee
        end

        @testset "Algebraic Relations" begin
            ex = σ*σ'
            @test isequal(ex, σgg)

            sz = σ'*σ - σ*σ'
            @test isequal(simplify(sz), simplify(2*σee - one(σgg)))
        end
    end

    @testset "Product Space Operations" begin
        ha1 = NLevelSpace(:atom, (:g, :e))
        ha2 = NLevelSpace(:atom, 2)
        @test ha1 != ha2

        hprod = ha1⊗ha2
        σ1 = Transition(hprod, :σ, :g, :e, 1)
        σ2 = Transition(hprod, :σ, 1, 2, 2)

        @testset "Basic Product Operations" begin
            @test σ1'*σ1 == Transition(hprod, :σ, :e, :e, 1)
            @test isequal(simplify(σ2*σ2'), simplify(1 - Transition(hprod, :σ, 2, 2, 2)))
            @test isequal(simplify(σ1*σ2), σ1*σ2)
        end

        @testset "Constructor Validation" begin
            @test_throws ErrorException Transition(hprod, :σ, :g, :e)
            @test Transition(hprod, :σ, :g, :e, 1) == σ1
            @test Transition(hprod, :σ, 1, 2, 2) == σ2
            @test_throws AssertionError Transition(hprod, :σ, 1, 2, 1)
        end

        @testset "Callable Constructor" begin
            @test_throws ErrorException Transition(hprod, :σ)
            s = Transition(hprod, :σ, 1)
            @test s isa SecondQuantizedAlgebra.CallableTransition
            @test s(:g, :e) == σ1
            @test acts_on(Transition(ha2, :s)) == 1
        end
    end

    @testset "Error Handling" begin
        @test_throws ArgumentError NLevelSpace(:atom, (:g, :e), 1)
    end
end
