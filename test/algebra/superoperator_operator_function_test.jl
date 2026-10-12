using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperExpr, QAdd, QExpr, sprepost, spre, spost,
    formal_terms, poly

@testset "Superoperators with formal operator functions" begin
    hφ = PhaseSpace(:t)
    φ = Position(hφ, :φ)
    π_ = Momentum(hφ, :π)
    @variables EJ::Real EC::Real γ::Real

    @testset "transmon Hamiltonian" begin
        H = 4EC * π_ * π_ - EJ * cos(φ)
        L = hamiltonian_action(H)
        @test L isa Lindbladian{QExpr}
        X = SuperExpr(L)
        @test !isempty(collect(formal_terms(X)))
        @test taylor(X, 0:4) == SuperOperator(hamiltonian_action(taylor(H, 0:4)))
        @test taylor(L, 0:4) == hamiltonian_action(taylor(H, 0:4))
        @test L(φ) isa QExpr
    end

    @testset "formal jump" begin
        D = γ * dissipator(cos(φ))
        @test D isa Lindbladian{QExpr}
        @test taylor(SuperExpr(D), 0:2) == SuperOperator(γ * dissipator(taylor(cos(φ), 0:2)))
        @test iszero(transpose(taylor(D, 0:3))(1))
    end

    @testset "key canonicality" begin
        @test sprepost(2cos(φ), φ) == 2 * sprepost(cos(φ), φ)
        @test sprepost(cos(φ), 3φ) == 3 * sprepost(cos(φ), φ)
        @test sprepost(cos(φ), φ) isa SuperExpr
        @test iszero(sprepost(cos(φ), φ) - sprepost(cos(φ), φ))
        X = sprepost(cos(φ), 1) + spre(φ)
        @test poly(X) == spre(φ)
    end

    @testset "promotion" begin
        Lf = hamiltonian_action(cos(φ))
        Lp = dissipator(φ)
        @test (Lf + Lp) isa Lindbladian{QExpr}
        @test (Lf ∘ Lp) isa SuperExpr
        @test (Lf ∘ Lp)(φ) == Lf(Lp(φ))
        @test SuperOperator(SuperExpr(Lp)) == SuperOperator(Lp)
        @test_throws ArgumentError SuperOperator(SuperExpr(Lf))
    end
end
