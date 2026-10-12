using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperExpr, QAdd, QExpr, sprepost, spre, spost,
    liouvillian, hamiltonian, channels, jump, partner, scope, LindbladChannel

# 𝒟[J, K] by hand: JρK† - ½K†Jρ - ½ρK†J.
raw_dissipator(J, K = J) = sprepost(J, K') - (1 // 2) * spre(K' * J) - (1 // 2) * spost(K' * J)

@testset "Lindbladian" begin
    h = FockSpace(:c)
    hp = h ⊗ FockSpace(:b)
    @qnumbers a::Destroy(h)
    a1 = Destroy(hp, :a, 1)
    b1 = Destroy(hp, :b, 2)
    @variables γ::Real κ::Real Δ::Real g::Real

    @testset "normalization" begin
        @test dissipator(2a) == 4 * dissipator(a)
        @test dissipator(-a) == dissipator(a)
        @test dissipator(im * a) == dissipator(a)
        @test iszero(dissipator(1))
        @test iszero(dissipator(0))
        L = dissipator(a + 3)
        @test only(channels(L)).first == only(channels(dissipator(a))).first
        # 𝒟[a + α] = 𝒟[a] - i[G, ·] with G = (i/2)(ᾱ a - α a†).
        @test hamiltonian(L) == (3im // 2) * (a - a')
        @test SuperOperator(L) == raw_dissipator(a + 3)
        @test iszero(hamiltonian_action(1))
        @test hamiltonian_action(a' * a + 5) == hamiltonian_action(a' * a)
        ch = only(channels(dissipator(a))).first
        @test ch isa LindbladChannel{QAdd}
        @test jump(ch) == a + 0
        @test partner(ch) === jump(ch)
        @test isempty(scope(ch))
    end

    @testset "cross terms" begin
        J = a1 + b1
        L = dissipator(J)
        Lsplit = dissipator(a1) + dissipator(b1) + dissipator(a1, b1) + dissipator(b1, a1)
        @test L != Lsplit
        @test SuperOperator(L) == SuperOperator(Lsplit)
        @test SuperOperator(dissipator(a1, b1)) == raw_dissipator(a1, b1)
        @test (L == SuperOperator(L)) == false
    end

    @testset "lowering" begin
        L = hamiltonian_action(Δ * a' * a) + γ * dissipator(a) + κ * dissipator(a')
        S = SuperOperator(L)
        Hs = Δ * a' * a
        expected = -im * spre(Hs) + im * spost(Hs) + γ * raw_dissipator(a) + κ * raw_dissipator(a')
        @test S == expected
        ρ = a' * a + a + a'
        @test L(ρ) == S(ρ)
    end

    @testset "trace preservation" begin
        models = (
            hamiltonian_action(Δ * a' * a + g * (a + a')) + γ * dissipator(a),
            κ * dissipator(a' * a) + dissipator(a + 2a'),
            dissipator(a1, b1) + dissipator(b1, a1) + dissipator(a1 * b1'),
        )
        for L in models
            @test iszero(transpose(L)(1))
            @test iszero(transpose(SuperOperator(L))(1))
        end
    end

    @testset "Hermiticity preservation" begin
        L = hamiltonian_action(Δ * a' * a + g * (a + a')) + γ * dissipator(a) + κ * dissipator(a' * a)
        ρ = a' * a + a + a' + a' * a' * a * a
        @test L(ρ) == L(ρ)'
    end

    @testset "liouvillian" begin
        H = Δ * a1' * a1
        @test liouvillian(H, [a1, b1]) == hamiltonian_action(H) + dissipator(a1) + dissipator(b1)
        @test liouvillian(H, [a1, b1]; rates = [γ, κ]) ==
            hamiltonian_action(H) + γ * dissipator(a1) + κ * dissipator(b1)
        Γ = [1 2im; -2im 3]
        Lm = liouvillian(H, [a1, b1]; rates = Γ)
        @test SuperOperator(Lm) == SuperOperator(
            hamiltonian_action(H) + dissipator(a1) + 3 * dissipator(b1) +
                2im * dissipator(a1, b1) - 2im * dissipator(b1, a1),
        )
        @test liouvillian(H, QAdd[]) == hamiltonian_action(H)
        @test_throws ArgumentError liouvillian(H, [a1, b1]; rates = [γ])
        @test_throws ArgumentError liouvillian(H, [a1, b1]; rates = ones(3, 3))
    end

    @testset "identity and promotion" begin
        L = hamiltonian_action(a' * a) + dissipator(a)
        @test_throws ArgumentError one(Lindbladian)
        @test L^0 == one(SuperOperator)
        @test 1 + L isa SuperOperator
        @test (L + sprepost(a, a)) isa SuperOperator
        @test (L ∘ L) isa SuperOperator
        @test transpose(L) isa SuperOperator
        @test L' isa SuperOperator
        @test (L ∘ L)(a) == L(L(a))
        @test 2 * L == L + L
        @test -L + L == zero(L)
    end

    @testset "return types" begin
        hφ = PhaseSpace(:t)
        φ = Position(hφ, :φ)
        Lp = hamiltonian_action(φ * φ) + dissipator(φ)
        Lf = hamiltonian_action(cos(φ)) + dissipator(φ)
        Sp = SuperOperator(Lp)
        X = SuperExpr(Lf)
        @test Lp isa Lindbladian{QAdd}
        @test Lf isa Lindbladian{QExpr}
        @test @inferred(Lp + Lp) isa Lindbladian{QAdd}
        @test @inferred(Lf + Lf) isa Lindbladian{QExpr}
        @test @inferred(Lp + Lf) isa Lindbladian{QExpr}
        @test @inferred(Lp + Sp) isa SuperOperator
        @test @inferred(Lf + Sp) isa SuperExpr
        @test @inferred(2 * Lp) isa Lindbladian{QAdd}
        @test @inferred(2 * Lf) isa Lindbladian{QExpr}
        @test @inferred(2 * Sp) isa SuperOperator
        @test @inferred(2 * X) isa SuperExpr
        @test @inferred(2 + Lp) isa SuperOperator
        @test @inferred(2 + Lf) isa SuperExpr
        @test @inferred(Lp ∘ Lp) isa SuperOperator
        @test @inferred(Lf ∘ Lp) isa SuperExpr
        @test @inferred(Sp ∘ X) isa SuperExpr
        @test @inferred(transpose(Lp)) isa SuperOperator
        @test @inferred(transpose(Lf)) isa SuperExpr
        @test @inferred(adjoint(Lf)) isa SuperExpr
        @test @inferred(Lp(φ + 0)) isa QAdd
        @test @inferred(Lf(φ + 0)) isa QExpr
        @test @inferred(Sp(cos(φ))) isa QExpr
        @test @inferred(X(φ + 0)) isa QExpr
        @test isone(one(SuperExpr))
    end
end
