using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, map_coefficients, to_cnum

@testset "Coefficient maps" begin
    hc = FockSpace(:c)
    ha = NLevelSpace(:atom, 2, 1)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    @variables N g κ x
    σ(p, q, k) = IndexedOperator(Transition(h, :σ, p, q, 2), k)
    i = Index(h, :i, N, 2)
    k = Index(h, :k, N, 2)
    gi = IndexedVariable(:g, i)
    sub(rules) = c -> substitute(c, rules)

    @testset "QAdd" begin
        q = x * a' * a + κ * a + 3
        @test isequal(map_coefficients(sub(Dict(x => 2)), q), 2 * a' * a + κ * a + 3)
        @test isequal(map_coefficients(sub(Dict(x => 0)), q), κ * a + 3)
        @test iszero(map_coefficients(c -> 0, q))
        @test isequal(map_coefficients(identity, q), q)
        @test isequal(map_coefficients(c -> 1 // 3, a), (1 // 3) * a)
        @test isequal(map_coefficients(c -> g, a), g * a)
        @test isequal(map_coefficients(c -> to_cnum(2), a' * a), 2 * a' * a)
        @test isequal(map_coefficients(c -> c^2, x * a + 2 * a'), x^2 * a + 4 * a')
        @test iszero(map_coefficients(simplify, (g^2 - g * g) * a))
    end

    @testset "summation scopes" begin
        rules = Dict(gi => g)
        q = Σ(gi * a' * a, i) + Σ(gi * a' * σ(1, 2, i), i)
        r = map_coefficients(sub(rules), q)
        @test isequal(r, N * g * a' * a + Σ(g * a' * σ(1, 2, i), i))
        @test isequal(r, substitute(q, rules))
        qc = Σ(gi * a' * a * σ(2, 2, k), i, [k])
        @test isequal(map_coefficients(sub(rules), qc), substitute(qc, rules))
        @test isequal(map_coefficients(sub(rules), qc), (N - 1) * g * a' * a * σ(2, 2, k))
        unit = cos(gi)^2 + sin(gi)^2
        @test isequal(simplify(Σ(unit * a' * a, i)), N * a' * a)
        @test simplify(Σ(unit * spre(a), i)) == N * spre(a)
        q2 = Σ(a' * σ(1, 2, i), i) + a' * a
        @test_throws ArgumentError map_coefficients(c -> c * gi, q2)
    end

    @testset "SuperOperator" begin
        S = x * sprepost(a, a') + κ * spre(a' * a) - spost(a)
        @test map_coefficients(sub(Dict(x => 2)), S) ==
            2 * sprepost(a, a') + κ * spre(a' * a) - spost(a)
        @test map_coefficients(sub(Dict(κ => 0)), S) == x * sprepost(a, a') - spost(a)
        @test iszero(map_coefficients(c -> 0, S))
        @test iszero(map_coefficients(simplify, (g^2 - g * g) * spre(a)))
        Si = Σ(sprepost(gi * a, a'), i) + Σ(sprepost(σ(1, 2, i), σ(2, 1, i)), i)
        rules = Dict(gi => g)
        @test map_coefficients(sub(rules), Si) == substitute(Si, rules)
        @test map_coefficients(sub(rules), Si) ==
            N * g * sprepost(a, a') + Σ(sprepost(σ(1, 2, i), σ(2, 1, i)), i)
        @test_throws ArgumentError map_coefficients(c -> c * gi, Σ(spre(σ(1, 2, i)), i) + spre(a))
    end

    @testset "Lindbladian" begin
        L = hamiltonian_action(x * a' * a) + g * dissipator(a + x * a') + x * dissipator(a) +
            Σ(gi * dissipator(σ(1, 2, i)), i)
        rules = Dict(x => 2, gi => g)
        Lm = map_coefficients(sub(rules), L)
        @test Lm isa Lindbladian
        @test Lm == substitute(L, rules)
        @test SuperOperator(Lm) == map_coefficients(sub(rules), SuperOperator(L))
        @test iszero(map_coefficients(c -> 0, L))
        @test map_coefficients(sub(Dict(x => 0)), L) ==
            g * dissipator(a) + Σ(gi * dissipator(σ(1, 2, i)), i)
    end
end
