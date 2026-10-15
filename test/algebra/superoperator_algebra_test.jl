using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperTerm, sprepost, spre, spost, superterm,
    left_operator, right_operator, compose, QAdd

@testset "SuperOperator algebra" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    @variables g Δ

    A = sprepost(a, a')
    B = spre(a' * a) + spost(a)
    C = g * sprepost(a', a) - spre(a)
    ρ = a' + a * a + 2

    @testset "constructors" begin
        @test @inferred(sprepost(a, a')) isa SuperOperator
        @test @inferred(spre(a)) == sprepost(a, 1)
        @test @inferred(spost(a)) == sprepost(1, a)
        @test spre(a)(ρ) == a * ρ
        @test spost(a)(ρ) == ρ * a
        @test A(ρ) == a * ρ * a'
        @test sprepost(2a, 3a') == 6 * sprepost(a, a')
        @test sprepost(a + a', 1) == spre(a) + spre(a')
        @test iszero(sprepost(0, a))
        @test iszero(zero(SuperOperator))
        Z = zero(SuperOperator)
        @test zero(A) == Z && Z == 0 * spre(1) && Z == A - A
        @test hash(Z) == hash(A - A)
        @test iszero(Z(ρ)) && Z + A == A
        @test sum(SuperOperator[]) == Z
        @test sum((sprepost(a, a') for _ in 1:0); init = zero(SuperOperator)) == Z
        @test sum(c * sprepost(a, a') for c in 1:2; init = zero(SuperOperator)) == 3A
        acc = zero(SuperOperator)
        for c in 1:2
            acc += c * sprepost(a, a')
        end
        @test acc == 3A
        @test isone(one(SuperOperator))
        @test one(SuperOperator)(ρ) == ρ
        t = superterm(Op[a], Op[a'])
        @test t isa SuperTerm
        @test SuperOperator([t => 1], Index[]) == A
        @test left_operator(t) == a + 0
        @test right_operator(t) == a' + 0
        @test length(A) == 1
        @test only(collect(A)).first == t
    end

    @testset "composition" begin
        @test (A ∘ B)(ρ) == A(B(ρ))
        @test (A * B)(ρ) == A(B(ρ))
        @test compose(A, B) == A ∘ B
        @test ((A ∘ B) ∘ C) == (A ∘ (B ∘ C))
        @test (A ∘ one(SuperOperator)) == A
        @test spre(a) ∘ spre(a') == spre(a * a')
        @test spost(a) ∘ spost(a') == spost(a' * a)
    end

    @testset "adjoint and transpose" begin
        @test (A ∘ B)' == B' ∘ A'
        @test transpose(A ∘ B) == transpose(B) ∘ transpose(A)
        @test A'' == A
        @test transpose(transpose(C)) == C
        @test transpose(sprepost(a, a')) == sprepost(a', a)
        @test sprepost(a, a')' == sprepost(a', a)
    end

    @testset "powers and scalars" begin
        @test A^0 == one(SuperOperator)
        @test A^1 == A
        @test A^2 == A ∘ A
        @test_throws ArgumentError A^-1
        @test (2 + A)(ρ) == 2ρ + A(ρ)
        @test (A - 1)(ρ) == A(ρ) - ρ
        @test (A + B)(ρ) == A(ρ) + B(ρ)
        @test (Δ * A)(ρ) == Δ * A(ρ)
        @test A / 2 == (1 // 2) * A
        @test iszero(A - A)
    end

    @testset "commutators" begin
        @test commutator(A, B) == A ∘ B - B ∘ A
        @test anticommutator(A, B) == A ∘ B + B ∘ A
    end

    @testset "tensor product of disjoint maps" begin
        hb = FockSpace(:b)
        hp = h ⊗ hb
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        X = sprepost(a1, a1') ⊗ spre(b1)
        @test X == sprepost(a1 * b1, a1')
        @test acts_on(X) == [1, 2]
    end

    @testset "structural passes" begin
        S = g * sprepost(a * a', a)
        @test normal_order(S) == S
        @test simplify(S) == S
        @test substitute(S, Dict(g => 2)) == 2 * sprepost(a * a', a)
        @test isequal(get_variables(S), get_variables(g * a))
        @test isequal(hash(S), hash(g * sprepost(a' * a + 1, a)))
    end
end
