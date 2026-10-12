using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperExpr, sprepost, spre, spost,
    poly, formal_terms, trigonometric_form, liouvillian, get_sum_indices, constraint_pairs

@testset "Superoperator structure" begin
    hc = FockSpace(:c)
    ha = NLevelSpace(:atom, 2, 1)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    @variables N g κ γ Δ φ
    σ(x, y, q) = IndexedOperator(Transition(h, :σ, x, y, 2), q)
    i = Index(h, :i, N, 2)
    k = Index(h, :k, N, 2)
    l = Index(h, :l, N, 2)
    ρ = a' * a + σ(2, 1, k) + 2

    @testset "SuperOperator queries" begin
        S = g * sprepost(a, a') + κ * spre(a' * a)
        @test Set(get_operators(S)) == Set([a, a'])
        @test Set(string.(Symbolics.get_variables(S))) == Set(["g", "κ"])
        @test isempty(get_indices(S)) && isempty(get_sum_indices(S))
        @test length(keys(S)) == 2
        @test isequal(real(get_prefactor(g * spre(a))), g)
        @test_throws ArgumentError get_prefactor(S)

        Si = Σ(sprepost(σ(1, 2, i), σ(2, 1, i)), i, [k])
        @test get_sum_indices(Si) == [i]
        @test i in get_indices(Si)
        @test (i, k) in constraint_pairs(Si) || (k, i) in constraint_pairs(Si)
        two = Σ(sprepost(σ(1, 2, i), σ(2, 1, l)), i, l)
        @test two == Σ(Σ(sprepost(σ(1, 2, i), σ(2, 1, l)), i), l)
    end

    @testset "SuperOperator arithmetic" begin
        S = g * sprepost(a, a') + spre(a)
        @test +S == S
        @test zero(S) == zero(SuperOperator) && one(S) == one(SuperOperator)
        @test (S / κ)(ρ) == S(ρ) * (1 / κ)
        @test (2 - S)(ρ) == 2 * ρ - S(ρ)
        @test dagger(S) == adjoint(S)
        @test iszero(simplify(expand((g + κ)^2 * spre(a)) - (g^2 + 2g * κ + κ^2) * spre(a)))
    end

    @testset "SuperOperator scope passes" begin
        P = spre(σ(1, 1, k)) + sprepost(σ(1, 1, k), σ(2, 2, l))
        @test expand_completeness(P)(ρ) == expand_completeness(P(ρ))
        D = spre(σ(1, 2, k) * σ(2, 1, l))
        @test expand_completeness(assume_distinct_index(D, [(k, l)])(ρ)) ==
            assume_distinct_index(D(ρ), [(k, l)])
    end

    @testset "Lindbladian queries" begin
        L = hamiltonian_action(Δ * a' * a) + κ * dissipator(a) +
            Σ(γ * dissipator(σ(1, 2, i)), i, [k])
        @test get_sum_indices(L) == [i]
        @test i in get_indices(L)
        @test Set(get_operators(L)) == Set([a, a', σ(1, 2, i)])
        @test acts_on(L) == [1, 2]
        @test issubset(["Δ", "κ", "γ"], string.(Symbolics.get_variables(L)))
        @test (i, k) in constraint_pairs(L) || (k, i) in constraint_pairs(L)
        @test hash(L) == hash(
            hamiltonian_action(Δ * a' * a) + κ * dissipator(a) +
                Σ(γ * dissipator(σ(1, 2, i)), i, [k])
        )
        @test_throws ArgumentError one(L)
    end

    @testset "Lindbladian passes match the flattened map" begin
        L = hamiltonian_action(Δ * a' * a) + κ * dissipator(a + a') + g * dissipator(σ(1, 2, k))
        S = SuperOperator(L)
        @test SuperOperator(normal_order(L)) == normal_order(S)
        @test iszero(simplify(SuperOperator(expand(L)) - expand(S)))
        rules = Dict(κ => 2, Δ => 3)
        @test SuperOperator(substitute(L, rules)) == substitute(S, rules)
        @test SuperOperator(expand_completeness(L)) == expand_completeness(S)

        Ld = dissipator(σ(1, 2, k) * σ(2, 1, l))
        @test expand_completeness(SuperOperator(assume_distinct_index(Ld, [(k, l)]))) ==
            assume_distinct_index(SuperOperator(Ld), [(k, l)])

        Lφ = hamiltonian_action(exp(im * φ) * a + exp(-im * φ) * a')
        @test iszero(simplify(SuperOperator(trigonometric_form(Lφ)) - trigonometric_form(SuperOperator(Lφ))))
    end

    @testset "Lindbladian arithmetic" begin
        L = hamiltonian_action(Δ * a' * a) + dissipator(a)
        @test +L == L
        @test SuperOperator(L * g) == SuperOperator(g * L)
        @test L + 2 == SuperOperator(L) + 2
        @test dagger(L) == adjoint(L)
        @test liouvillian(Δ * a' * a) == hamiltonian_action(Δ * a' * a)
        @test SuperOperator(Σ(κ * dissipator(a), i)) == Σ(SuperOperator(κ * dissipator(a)), i)
    end

    @testset "SuperExpr" begin
        hφ = PhaseSpace(:t)
        x = Position(hφ, :x)
        X = g * spre(cos(x)) + spost(x)
        @test X isa SuperExpr
        @test SuperExpr(X) === X
        @test poly(X) == spost(x)
        @test length(collect(formal_terms(X))) == 1
        @test isequal(only(keys(X.formal)), only(keys(X.formal)))
        @test x in get_operators(X)
        @test acts_on(X) == [1]
        @test isempty(get_indices(X)) && isempty(get_sum_indices(X))
        @test "g" in string.(Symbolics.get_variables(X))
        @test hash(X) == hash(g * spre(cos(x)) + spost(x))
        @test iszero(zero(SuperExpr)) && iszero(zero(X))

        @test normal_order(X) == X
        @test simplify(expand(X)) == simplify(X)
        @test substitute(X, Dict(g => 2)) == 2 * spre(cos(x)) + spost(x)
        @test +X == X
        @test X * κ == κ * X
        @test X / 2 == X * (1 // 2)
        @test dagger(X) == adjoint(X)
        @test (X - spre(x))(1) == X(1) - x
        @test_throws ArgumentError inv(X)
        @test spre(x)(cos(x)) == x * cos(x)
    end
end
