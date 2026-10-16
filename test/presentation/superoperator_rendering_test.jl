using SecondQuantizedAlgebra
using Latexify: latexify
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperExpr, sprepost, spre, spost, channels

@testset "Superoperator rendering" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    @variables γ

    @testset "Lindbladian" begin
        L = hamiltonian_action(a' * a) + dissipator(a)
        @test repr(L) == "-i[a' * a, ·] + 𝒟[a]"
        @test repr(γ * dissipator(a) - 2 * dissipator(a')) == "γ * 𝒟[a] - 2 * 𝒟[a']"
        @test repr(zero(L)) == "0"
        @test repr(dissipator(a, a')) == "𝒟[a, a']"
        @test repr(only(channels(dissipator(a))).first) == "𝒟[a]"
        @test repr("text/latex", L) == "\$-i\\left[a^{\\dagger}a,\\cdot\\right] + \\mathcal{D}\\left[a\\right]\$"
        @test repr("text/latex", γ * dissipator(a)) == "\$\\gamma\\,\\mathcal{D}\\left[a\\right]\$"
    end

    @testset "SuperOperator" begin
        @test repr(sprepost(a, a')) == "a ⋅ (·) ⋅ a'"
        @test repr(spre(a)) == "a ⋅ (·)"
        @test repr(spost(a)) == "(·) ⋅ a"
        @test repr(one(SuperOperator)) == "(·)"
        @test repr(zero(SuperOperator)) == "0"
        @test repr(SuperOperator(hamiltonian_action(a' * a) + dissipator(a))) ==
            "(-1//2 + 1im) * (·) ⋅ a' * a + a ⋅ (·) ⋅ a' - (1//2 + 1im) * a' * a ⋅ (·)"
        @test repr(sprepost(a, a') - (1 // 2) * spre(a' * a)) == "a ⋅ (·) ⋅ a' - 1//2 * a' * a ⋅ (·)"
        @test repr("text/latex", sprepost(a, a')) == "\$a\\,\\cdot\\,a^{\\dagger}\$"
        @test latexify(spre(a)) isa AbstractString
        @test repr((γ + im * γ^2) * spre(a)) == "(γ + (γ^2)im) * a ⋅ (·)"
        @test startswith(repr("text/latex", [spre(a), spost(a)]), "\$\$ \\left[")
    end

    @testset "indexed" begin
        ha = NLevelSpace(:atom, 2)
        @variables N
        i = Index(ha, :i, N, ha)
        j = Index(ha, :j, N, ha)
        σ(x, y, q) = IndexedOperator(Transition(ha, :σ, x, y), q)
        L = dissipator(σ(1, 2, i), DoubleIndexedVariable(:Γ, i, j), i => j)
        @test repr(L) == "Σ(i=1:N) Γ(i, i) * 𝒟[σ_i₁₂] + Σ(i=1:N)Σ(j=1:N)(i≠j) Γ(i, j) * 𝒟[σ_i₁₂, σ_j₁₂]"
        @test occursin("\\mathcal{D}\\left[{\\sigma}_{i}^{{12}},{\\sigma}_{j}^{{12}}\\right]", repr("text/latex", L))
        @test startswith(repr(SuperOperator(Σ(dissipator(σ(1, 2, i)), i))), "Σ(i=1:N)")
    end

    @testset "SuperExpr" begin
        φ = Position(PhaseSpace(:t), :φ)
        X = SuperExpr(hamiltonian_action(cos(φ)))
        @test repr(X) == "im * (·) ⋅ cos(φ) - im * cos(φ) ⋅ (·)"
        @test repr("text/latex", X) ==
            "\$\\mathit{i} \\cdot\\,\\cos\\left( \\hat{\\varphi} \\right) - \\mathit{i} \\cos\\left( \\hat{\\varphi} \\right)\\,\\cdot\$"
        @test repr(hamiltonian_action(cos(φ))) == "-i[cos(φ), ·]"
    end
end
