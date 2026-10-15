using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: Coeff, coefficient, sprepost, spre, spost, to_num

same(c, x) = isequal(to_num(c), to_num(x * one(Coeff)))

@testset "operator coefficient" begin
    @variables Δ F g κ
    h = FockSpace(:a) ⊗ FockSpace(:b)
    a = Destroy(h, :a, 1)
    b = Destroy(h, :b, 2)

    @testset "driven oscillator" begin
        H = Δ * a' * a + F * a' + conj(F) * a + 7
        @test coefficient(H, a' * a) isa Coeff
        @test same(coefficient(H, a' * a), Δ)
        @test same(coefficient(H, a'), F)
        @test same(coefficient(H, a), conj(F))
        @test same(coefficient(H, 1), 7)
        @test iszero(coefficient(H, a' * a' * a))
        @test iszero(coefficient(H, b))
        @test coefficient(H - 7, 1) isa Coeff
        @test iszero(coefficient(H - 7, 1))
    end

    @testset "canonical order of the term" begin
        H = g * (a' * b + b' * a) + κ * a * b
        @test same(coefficient(H, b' * a), g)
        @test same(coefficient(H, a * b'), g)
        @test same(coefficient(H, b * a), κ)
        @test same(coefficient(H, a * b), κ)
    end

    @testset "rebuilds the expression" begin
        H = Δ * a' * a + 2 // 3 * a' * a' * b + F * b' + 1
        terms = (a' * a, a' * a' * b, b', 1)
        @test sum(coefficient(H, t) * t for t in terms) == H
    end

    @testset "operand kinds" begin
        @test same(coefficient(a, a), 1)
        @test iszero(coefficient(a, 1))
        @test same(coefficient(3, 1), 3)
        @test same(coefficient(g, 1), g)
        @test iszero(coefficient(g, a))
    end

    @testset "same-site composition" begin
        ha = NLevelSpace(:atom, 2, 1)
        σ(i, j) = Transition(ha, :σ, i, j)
        H = Δ * σ(2, 2) + g * σ(1, 1)
        @test same(coefficient(H, σ(2, 1) * σ(1, 2)), Δ)
        @test same(coefficient(H, σ(1, 1)), g)
    end

    @testset "invalid terms" begin
        H = Δ * a' * a + F * a'
        @test_throws ArgumentError coefficient(H, 2a)
        @test_throws ArgumentError coefficient(H, g * a)
        @test_throws ArgumentError coefficient(H, a + a')
        @test_throws ArgumentError coefficient(H, a * a')
        @test_throws ArgumentError coefficient(H, 0)
        @test_throws ArgumentError coefficient(H, 2)
        @test_throws ArgumentError coefficient(cos(a + a'), a)
        @test_throws ArgumentError coefficient(H, cos(a + a'))
    end

    @testset "index constraints" begin
        hs = FockSpace(:s)
        i = Index(hs, :i, 3, hs)
        j = Index(hs, :j, 3, hs)
        c = Destroy(hs, :c)
        ci = IndexedOperator(c, i)
        cj = IndexedOperator(c, j)
        @test same(coefficient(2 * ci' * cj, ci' * cj), 2)
        @test iszero(coefficient(2 * ci' * cj, ci' * ci))

        X = assume_distinct_index(3 * ci' * cj, [(i, j)])
        constrained = assume_distinct_index(ci' * cj, [(i, j)])
        @test same(coefficient(X, constrained), 3)
        @test_throws ArgumentError coefficient(X, ci' * cj)
        @test_throws ArgumentError coefficient(2 * ci' * cj, constrained)
        @test_throws ArgumentError coefficient(X + ci' * cj, constrained)

        @test_throws ArgumentError coefficient(Σ(ci' * ci, i), ci' * ci)
        @test_throws ArgumentError coefficient(ci' * ci, Σ(ci' * ci, i))
    end

    @testset "superoperator sandwiches" begin
        S = 2 * sprepost(a, a') - spre(a' * a) + Δ * spost(a' * a) + κ * sprepost(b, 1)
        @test coefficient(S, sprepost(a, a')) isa Coeff
        @test same(coefficient(S, sprepost(a, a')), 2)
        @test same(coefficient(S, spre(a' * a)), -1)
        @test same(coefficient(S, spost(a' * a)), Δ)
        @test same(coefficient(S, spre(b)), κ)
        @test iszero(coefficient(S, sprepost(a', a)))
        @test iszero(coefficient(S, one(S)))
        @test_throws ArgumentError coefficient(S, 2 * sprepost(a, a'))
        @test_throws ArgumentError coefficient(S, sprepost(a, a') + spre(a))
        @test_throws ArgumentError coefficient(S, zero(S))
    end
end
