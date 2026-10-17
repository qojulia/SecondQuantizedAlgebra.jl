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
        # Loss terms are split into the same additive atoms as the gain term.
        @test iszero(simplify(transpose(dissipator(cos(φ) + sin(φ)))(1)))
        @test iszero(simplify(transpose(dissipator(cos(φ) + 2sin(φ), sin(φ) + φ))(1)))
        shifted = dissipator(cos(φ), φ + 1)
        @test taylor(SuperExpr(shifted), 0:2) == SuperOperator(dissipator(taylor(cos(φ), 0:2), φ + 1))
        @test iszero(simplify(transpose(shifted)(1)))
    end

    @testset "key canonicality" begin
        @test sprepost(2cos(φ), φ) == 2 * sprepost(cos(φ), φ)
        @test sprepost(cos(φ), 3φ) == 3 * sprepost(cos(φ), φ)
        @test sprepost(cos(φ), φ) isa SuperExpr
        @test iszero(sprepost(cos(φ), φ) - sprepost(cos(φ), φ))
        X = sprepost(cos(φ), 1) + spre(φ)
        @test poly(X) == spre(φ)
    end

    @testset "mixed sides split into atoms" begin
        X = sprepost(cos(φ) + φ, π_ + 2)
        @test X isa SuperExpr
        @test poly(X) == sprepost(φ, π_) + 2 * spre(φ)
        @test X == sprepost(cos(φ), π_) + 2 * sprepost(cos(φ), 1) + sprepost(φ, π_) + 2 * spre(φ)
        @test X(φ) == cos(φ) * φ * π_ + 2 * cos(φ) * φ + φ * φ * π_ + 2 * φ * φ
    end

    @testset "adjoint and scalar division" begin
        X = sprepost(cos(φ), π_)
        @test qadjoint(im * X) == -im * X
        @test qadjoint(X + im * spre(φ)) == X - im * spre(φ)
        @test X / γ == inv(γ) * X
        @test X // 3 == (1 // 3) * X
        @test X // γ == inv(γ) * X
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
