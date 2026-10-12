# Numeric superoperators and master-equation solver inputs on both backends.

using SecondQuantizedAlgebra
import QuantumOpticsBase as QOB
import QuantumOptics as QO
import QuantumToolbox as QTB
import SciMLOperators as SO
using SymbolicUtils: SymbolicUtils
using Symbolics: @variables, Num
import LinearAlgebra: tr, I, exp
using Test

import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, liouvillian, to_numeric_master

const QOBB = QuantumOpticsBackend()
const QTBB = QuantumToolboxBackend()

apply_numeric(S, ρ::AbstractMatrix) = reshape(S * vec(ρ), size(ρ))

@testset "Superoperator numeric conversion" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    @variables Δ γ
    b = QOB.FockBasis(4)
    an = QOB.destroy(b)
    n = length(b)

    @testset "QuantumOptics superoperators" begin
        L = hamiltonian_action(2 * a' * a) + (1 // 2) * dissipator(a)
        ref = QOB.liouvillian(2 * QOB.dagger(an) * an, [an]; rates = [0.5])
        Ln = to_numeric(L, b)
        @test Ln isa QOB.SparseSuperOpType
        @test to_numeric(L, b; op_type = identity) isa QOB.SparseSuperOpType
        @test Ln.data ≈ ref.data
        @test to_numeric(SuperOperator(L), b).data ≈ ref.data
        @test to_numeric(sprepost(a, a'), b).data ≈ QOB.sprepost(an, QOB.dagger(an)).data
        @test iszero(to_numeric(zero(SuperOperator), b).data)
        @test to_numeric(one(SuperOperator), b).data ≈ I(n^2)
        @test to_numeric(L, QOB.fockstate(b, 1)).data ≈ ref.data
        @test to_numeric(L, h, 4; backend = QOBB).data ≈ ref.data
        @test to_numeric(spre(a), b, Dict(a => 2an)).data ≈ QOB.spre(2an).data
    end

    @testset "complex Hermitian rate matrix" begin
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        bp = QOB.FockBasis(2) ⊗ QOB.FockBasis(2)
        J = [QOB.embed(bp, 1, QOB.destroy(QOB.FockBasis(2))), QOB.embed(bp, 2, QOB.destroy(QOB.FockBasis(2)))]
        Γ = [1.0 0.3 + 0.4im; 0.3 - 0.4im 0.8]
        L = liouvillian(a1' * a1, [a1, b1]; rates = Γ)
        ref = QOB.liouvillian(QOB.dagger(J[1]) * J[1], J; rates = Γ)
        @test to_numeric(L, bp).data ≈ ref.data
        H, c_ops = to_numeric_master(L, hp, (2, 2); backend = QOBB)
        @test length(c_ops) == 2
        @test QOB.liouvillian(H, c_ops).data ≈ ref.data
    end

    @testset "parameters" begin
        L = hamiltonian_action(Δ * a' * a) + γ * dissipator(a)
        ref = QOB.liouvillian(1.5 * QOB.dagger(an) * an, [an]; rates = [0.25])
        @test to_numeric(L, b; parameter = Dict(Δ => 1.5, γ => 0.25)).data ≈ ref.data
        H, c_ops = to_numeric_master(L, h, 4; backend = QOBB, parameter = Dict(Δ => 1.5, γ => 0.25))
        @test QOB.liouvillian(H, c_ops).data ≈ ref.data
    end

    @testset "duality" begin
        S = sprepost(a, a') + 2 * spre(a' * a) - spost(a * a) + (1 // 3) * sprepost(a', a * a)
        Sn = to_numeric(S, b).data
        Stn = to_numeric(transpose(S), b).data
        ρ = rand(ComplexF64, n, n)
        A = rand(ComplexF64, n, n)
        @test tr(A * apply_numeric(Sn, ρ)) ≈ tr(apply_numeric(Stn, A) * ρ)
    end

    @testset "indexed sites" begin
        ha = NLevelSpace(:atom, 2)
        @variables N
        i = Index(ha, :i, N, ha)
        j = Index(ha, :j, N, ha)
        σ(x, y, q) = IndexedOperator(Transition(ha, :σ, x, y), q)
        Γsym = DoubleIndexedVariable(:Γ, i, j)
        L = dissipator(σ(1, 2, i), Γsym, i => j)
        Γf = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
            :Γ; type = SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing}, shape = UnitRange{Int}[],
        )
        Γv = [1.0 0.3 0.1; 0.3 0.7 0.2; 0.1 0.2 0.9]
        subs = Dict{Num, Float64}(Num(Γf(p, q)) => Γv[p, q] for p in 1:3, q in 1:3)
        b1 = QOB.NLevelBasis(2)
        bb = b1 ⊗ b1 ⊗ b1
        sites = Dict{Int, Vector{Int}}(1 => [1, 2, 3])
        J = [QOB.embed(bb, k, QOB.transition(b1, 1, 2)) for k in 1:3]
        ref = QOB.liouvillian(0 * QOB.dagger(J[1]) * J[1], J; rates = Γv)
        @test to_numeric(L, bb, sites, Dict{SecondQuantizedAlgebra.QSym, Any}(), subs).data ≈ ref.data
        H, c_ops = to_numeric_master(L, ha ⊗ ha ⊗ ha, (2, 2, 2); backend = QOBB, sites, scalar_subs = subs)
        @test length(c_ops) == 3
        @test QOB.liouvillian(H, c_ops).data ≈ ref.data

        # Sides that differ only in their unrolled site compare equal as `Vector{Op}`, so the
        # side-matrix cache must distinguish slots.
        Sx = Σ(sprepost(σ(1, 2, i), σ(2, 1, i)), i)
        refx = sum(QOB.sprepost(Jn, QOB.dagger(Jn)) for Jn in J)
        @test to_numeric(Sx, bb, sites).data ≈ refx.data

        # Independent decay against finite-site enumeration.
        Li = Σ(dissipator(σ(1, 2, i)), i)
        refi = QOB.liouvillian(0 * QOB.dagger(J[1]) * J[1], J)
        @test to_numeric(Li, bb, sites).data ≈ refi.data
    end

    @testset "collapse-operator factorization" begin
        # A rank-one block: one eigenvalue is zero up to roundoff and is clipped.
        hp = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hp, :a, 1)
        b1 = Destroy(hp, :b, 2)
        L = liouvillian(0, [a1, b1]; rates = [0.5 0.5; 0.5 0.5 + 1.0e-17])
        H, c_ops = to_numeric_master(L, hp, (2, 2); backend = QOBB)
        @test length(c_ops) == 1
        @test to_numeric(L, hp, (2, 2); backend = QOBB).data ≈ QOB.liouvillian(H, c_ops).data
        # Degenerate rates: any eigenbasis gives the same map.
        Ld = liouvillian(0, [a1, b1]; rates = [1.0 0.0; 0.0 1.0])
        Hd, cd = to_numeric_master(Ld, hp, (2, 2); backend = QOBB)
        @test length(cd) == 2
        @test to_numeric(Ld, hp, (2, 2); backend = QOBB).data ≈ QOB.liouvillian(Hd, cd).data
    end

    @testset "QuantumOptics master equation" begin
        # Driven damped cavity: compare timeevolution.master on the solver inputs with the
        # exact propagator exp(Lt) of the converted superoperator. A heating jump a† is left
        # out: the solver forms a a† in the truncated space, where it differs from the
        # normal-ordered a†a + 1 at the cutoff.
        par = Dict(Δ => 0.7, γ => 0.4)
        L = hamiltonian_action(Δ * a' * a + (1 // 5) * (a + a')) + γ * dissipator(a) + (1 // 10) * dissipator(a' * a)
        H, c_ops = to_numeric_master(L, h, 4; backend = QOBB, parameter = par)
        ρ0 = QOB.dm(QOB.fockstate(b, 2))
        tout, ρt = QO.timeevolution.master([0.0, 1.5], ρ0, H, c_ops; abstol = 1.0e-10, reltol = 1.0e-10)
        Ln = to_numeric(L, b; parameter = par).data
        exact = reshape(exp(1.5 * Matrix(Ln)) * vec(Matrix(ρ0.data)), n, n)
        @test ρt[end].data ≈ exact rtol = 1.0e-6
        # Pure loss: ⟨n⟩ = n₀ e^{-γt}.
        Hl, cl = to_numeric_master(γ * dissipator(a), h, 4; backend = QOBB, parameter = Dict(γ => 0.4))
        _, ρl = QO.timeevolution.master([0.0, 1.5], ρ0, Hl, cl; abstol = 1.0e-10, reltol = 1.0e-10)
        @test real(QOB.expect(QOB.dagger(an) * an, ρl[end])) ≈ 2 * exp(-0.4 * 1.5) rtol = 1.0e-6
    end

    @testset "QuantumToolbox" begin
        L = hamiltonian_action(2 * a' * a) + (1 // 2) * dissipator(a)
        aq = QTB.destroy(5)
        ref = QTB.liouvillian(2 * aq' * aq, [sqrt(0.5) * aq])
        Lq = to_numeric(L, h, 4; backend = QTBB)
        @test Lq isa QTB.QuantumObject{QTB.SuperOperator}
        @test Lq.data ≈ ref.data
        @test to_numeric(L, 5).data ≈ ref.data
        H, c_ops = to_numeric_master(L, h, 4; backend = QTBB)
        @test QTB.liouvillian(H, c_ops).data ≈ ref.data
        @test_throws ArgumentError to_numeric(SecondQuantizedAlgebra.SuperExpr(hamiltonian_action(cos(a + a'))), 5)
    end

    @testset "QuantumToolbox time dependence" begin
        @variables t κ(t)
        L = hamiltonian_action(Δ * a' * a) + κ * dissipator(a)
        tp = Dict(κ => (t -> 0.5 + 0.1t))
        par = Dict(Δ => 2.0)
        Lq = to_numeric(L, h, 4; backend = QTBB, parameter = par, time_parameter = tp)
        @test Lq isa QTB.QuantumObjectEvolution
        aq = QTB.destroy(5)
        M1 = SO.concretize(SO.update_coefficients(Lq.data, nothing, nothing, 1.0))
        @test M1 ≈ QTB.liouvillian(2 * aq' * aq, [sqrt(0.6) * aq]).data
        H, c_ops = to_numeric_master(L, h, 4; backend = QTBB, parameter = par, time_parameter = tp)
        ψ0 = QTB.fock(5, 3)
        e_ops = [aq' * aq]
        sol = QTB.mesolve(H, ψ0, [0.0, 1.0], c_ops; e_ops, progress_bar = Val(false))
        sol2 = QTB.mesolve(Lq, ψ0, [0.0, 1.0]; e_ops, progress_bar = Val(false))
        exact = 3 * exp(-(0.5 + 0.05))
        @test real(sol.expect[1, end]) ≈ exact rtol = 1.0e-5
        @test real(sol2.expect[1, end]) ≈ exact rtol = 1.0e-5
    end
end
