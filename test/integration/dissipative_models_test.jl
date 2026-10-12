# Dissipative models crossing the symbolic, indexed and numeric layers.

using SecondQuantizedAlgebra
import QuantumOpticsBase as QOB
using SymbolicUtils: SymbolicUtils
using Symbolics: @variables, Num
using Test

import SecondQuantizedAlgebra: SuperOperator, liouvillian, to_numeric_master, QSym

# Numeric Heisenberg drift Lᵀ(A) from the Schrödinger matrix M: vec(Bᵀ) = Mᵀ vec(Aᵀ).
heisenberg(M, A) = transpose(reshape(transpose(M) * vec(transpose(A)), size(A)))

@testset "Dissipative models" begin
    @testset "damped Jaynes-Cummings" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        a = Destroy(h, :a, 1)
        σ(x, y) = Transition(h, :σ, x, y, 2)
        @variables ωc ωa g κ γ
        H = ωc * a' * a + ωa * σ(2, 2) + g * (a' * σ(1, 2) + a * σ(2, 1))
        L = hamiltonian_action(H) + κ * dissipator(a) + γ * dissipator(σ(1, 2))
        @test transpose(L)(a' * a) == -κ * a' * a + im * g * (a * σ(2, 1) - a' * σ(1, 2))
        @test transpose(L)(σ(2, 2)) == -γ * σ(2, 2) + im * g * (a' * σ(1, 2) - a * σ(2, 1))
        @test isequal(
            average(transpose(L)(a' * a)),
            average(-κ * a' * a + im * g * (a * σ(2, 1) - a' * σ(1, 2))),
        )
        @test iszero(transpose(L)(1))

        par = Dict(ωc => 1.0, ωa => 1.2, g => 0.1, κ => 0.3, γ => 0.05)
        bh = QOB.FockBasis(3) ⊗ QOB.NLevelBasis(2)
        Ln = to_numeric(L, bh; parameter = par)
        Hn, c_ops = to_numeric_master(L, h, (3, nothing); backend = QuantumOpticsBackend(), parameter = par)
        @test QOB.liouvillian(Hn, c_ops).data ≈ Ln.data
        A = to_numeric(a' * a, bh).data
        D = to_numeric(transpose(L)(a' * a), bh; parameter = par).data
        @test heisenberg(Ln.data, Matrix(A)) ≈ D
    end

    @testset "superradiant ensemble" begin
        ha = NLevelSpace(:atom, 2)
        @variables N
        i = Index(ha, :i, N, ha)
        j = Index(ha, :j, N, ha)
        k = Index(ha, :k, N, ha)
        σ(x, y, q) = IndexedOperator(Transition(ha, :σ, x, y), q)
        b1 = QOB.NLevelBasis(2)
        bb = b1 ⊗ b1 ⊗ b1
        sites = Dict{Int, Vector{Int}}(1 => [1, 2, 3])
        Jn = [QOB.embed(bb, q, QOB.transition(b1, 1, 2)) for q in 1:3]
        Z = 0 * QOB.dagger(Jn[1]) * Jn[1]

        collective = dissipator(σ(1, 2, i), 1, i => j)
        independent = Σ(dissipator(σ(1, 2, i)), i)
        @test to_numeric(collective, bb, sites).data ≈ QOB.liouvillian(Z, [sum(Jn)]).data
        @test to_numeric(independent, bb, sites).data ≈ QOB.liouvillian(Z, Jn).data
        H, c_ops = to_numeric_master(collective, ha ⊗ ha ⊗ ha, (2, 2, 2); backend = QuantumOpticsBackend(), sites)
        @test length(c_ops) == 1
        @test QOB.liouvillian(H, c_ops).data ≈ QOB.liouvillian(Z, [sum(Jn)]).data

        # QuantumCumulants-style kernel with a complex Hermitian Γ(i,j), migrated as Γ(j,i).
        Γ = DoubleIndexedVariable(:Γ, i, j)
        L = dissipator(σ(1, 2, i), change_index(Γ, Dict(i => j, j => i)), i => j)
        Γf = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
            :Γ; type = SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing}, shape = UnitRange{Int}[],
        )
        Γv = [1.0 0.2 + 0.1im 0.05im; 0.2 - 0.1im 0.8 0.1; -0.05im 0.1 0.6]
        subs = Dict{Num, ComplexF64}(Num(Γf(p, q)) => Γv[p, q] for p in 1:3, q in 1:3)
        Ln = to_numeric(L, bb, sites, Dict{QSym, Any}(), subs).data
        qc = QOB.liouvillian(Z, Jn; rates = transpose(Γv)).data
        @test Ln ≈ qc
        excitation = Σ(σ(2, 2, k), k)
        drift = transpose(L)(excitation)
        Dn = to_numeric(drift, bb, sites, Dict{QSym, Any}(), subs).data
        An = to_numeric(excitation, bb, sites).data
        @test heisenberg(Ln, Matrix(An)) ≈ Dn
    end
end
