using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, liouvillian, channels,
    scope, get_sum_indices

@testset "Indexed superoperators" begin
    ha = NLevelSpace(:atom, 2)
    @variables N
    i = Index(ha, :i, N, ha)
    j = Index(ha, :j, N, ha)
    k = Index(ha, :k, N, ha)
    l = Index(ha, :l, N, ha)
    σ(x, y, q) = IndexedOperator(Transition(ha, :σ, x, y), q)

    @testset "independent decay" begin
        L = Σ(dissipator(σ(1, 2, i)), i)
        @test only(scope(only(channels(L)).first)) == i
        @test transpose(L)(σ(2, 2, k)) == -σ(2, 2, k)
        @test transpose(L)(σ(1, 2, k)) == -(1 // 2) * σ(1, 2, k)
        @test iszero(transpose(L)(1))
        # Pumping 𝒟[σ21ᵢ] keeps the trace.
        @test iszero(transpose(Σ(dissipator(σ(2, 1, i)), i))(1))
        @test SuperOperator(L) == Σ(SuperOperator(dissipator(σ(1, 2, i))), i)
    end

    @testset "collective decay" begin
        Γ = DoubleIndexedVariable(:Γ, i, j)
        L = dissipator(σ(1, 2, i), Γ, i => j)
        @test length(channels(L)) == 2
        @test iszero(transpose(L)(1))
        Lnd = dissipator(σ(1, 2, i), DoubleIndexedVariable(:Γ, i, j; identical = false), i => j)
        @test length(channels(Lnd)) == 1
        @test all(length(scope(ch)) == 2 for (ch, _) in channels(Lnd))
        @test iszero(transpose(Lnd)(1))
        @test_throws ArgumentError dissipator(σ(1, 2, i), Γ, j => i)
        @test_throws ArgumentError dissipator(σ(1, 2, i) * σ(1, 2, j), Γ, i => j)
    end

    @testset "multi-index families" begin
        hc = NLevelSpace(:atom, 2) ⊗ NLevelSpace(:mode, 2)
        p = Index(hc, :p, N, 1)
        q = Index(hc, :q, N, 1)
        r = Index(hc, :r, N, 2)
        s = Index(hc, :s, N, 2)
        τ(x, y, n) = IndexedOperator(Transition(hc, :τ, x, y, 1), n)
        μ(x, y, n) = IndexedOperator(Transition(hc, :μ, x, y, 2), n)
        L = dissipator(τ(1, 2, p) * μ(1, 2, r), 1, p => q, r => s)
        @test length(channels(L)) == 4
        @test iszero(transpose(L)(1))
    end

    @testset "scope rules" begin
        A = Σ(spre(σ(1, 2, i)), i)
        @test_throws ArgumentError A ∘ A
        @test_throws ArgumentError spre(σ(2, 1, i)) ∘ A
        @test_throws ArgumentError A ∘ spre(σ(2, 1, i))
        @test_throws ArgumentError sprepost(Σ(σ(1, 2, i), i), 1)
        @test_throws ArgumentError dissipator(Σ(σ(1, 2, i), i))
        @test_throws ArgumentError A^2
        L = Σ(dissipator(σ(1, 2, i)), i)
        @test_throws ArgumentError L(σ(2, 2, i))
        @test_throws ArgumentError L ∘ L
        @test A ∘ change_index(A, i, j) == Σ(Σ(spre(σ(1, 2, i) * σ(1, 2, j)), i), j)
    end

    @testset "range scaling and constrained cardinality" begin
        @test iszero(transpose(Σ(dissipator(σ(2, 1, i)), i))(1))
        S = Σ(spre(σ(1, 2, j)), i, [j])
        @test S == (N - 1) * spre(σ(1, 2, j))
        @test isempty(get_sum_indices(S))
        S2 = Σ(Σ(one(SuperOperator), i, [j]), j)
        @test S2 == N * (N - 1) * one(SuperOperator)
        # Two free partners not known to be distinct: the count is undecidable.
        S3 = Σ(spre(σ(1, 2, j) * σ(1, 2, k)), i, [j, k])
        @test i in get_sum_indices(S3)
    end

    @testset "matrix-rate orientation" begin
        # QuantumCumulants' kernel Γ(i,j) Jᵢ† A Jⱼ is SQA's Γ(j,i) 𝒟[Jᵢ, Jⱼ].
        h2 = NLevelSpace(:a, 2) ⊗ NLevelSpace(:b, 2)
        J = [Transition(h2, :s, 1, 2, 1), Transition(h2, :t, 1, 2, 2)]
        Γ = [1 (1 // 2 + 1im // 3); (1 // 2 - 1im // 3) 2]
        A = Transition(h2, :s, 2, 2, 1) + J[1]' * J[2]
        qc = sum(
            Γ[x, y] * (J[x]' * A * J[y] - (1 // 2) * (J[x]' * J[y] * A + A * J[x]' * J[y]))
                for x in 1:2, y in 1:2
        )
        L = liouvillian(0, J; rates = transpose(Γ))
        @test transpose(L)(A) == qc
    end

    @testset "index relabelling" begin
        L = Σ(dissipator(σ(1, 2, i)), i)
        @test change_index(L, i, k) == Σ(dissipator(σ(1, 2, k)), k)
        S = spre(σ(1, 2, i)) + spost(σ(2, 1, j))
        @test change_index(S, Dict(i => j, j => i)) == spre(σ(1, 2, j)) + spost(σ(2, 1, i))
    end
end
