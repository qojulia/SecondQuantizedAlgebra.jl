# QuantumToolbox backend for `to_numeric`/`numeric_average`.
#
# `QuantumToolbox` is imported (not `using`) so its `⊗`/`tensor`/`expect` do not collide
# with the QuantumInterface ones that `SecondQuantizedAlgebra` extends; `⊗` here is the
# shared QuantumInterface operator (works on Hilbert spaces and QuantumOptics bases). Parity
# is checked through `expect`/`numeric_average` scalars (never raw matrices, whose kron
# convention differs across backends), with states built per-backend in slot order.

using SecondQuantizedAlgebra
import QuantumToolbox as QTB
import QuantumOpticsBase as QOB
import SciMLOperators as SO
using Symbolics: @variables
import LinearAlgebra: tr, norm
using Test

const QTBB = QuantumToolboxBackend()

@testset "QuantumToolbox backend" begin
    @testset "lazy + n-stable + static==td" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        H2 = to_numeric(a' * a, h, 5; backend = QTBB)
        H3 = to_numeric(a' * a + a + a', h, 5; backend = QTBB)
        @variables E::Number
        Htd = to_numeric(E * a + conj(E) * a', h, 5; backend = QTBB, time_parameter = Dict(E => t -> 1 + im * t))

        # The vector-backed VecSum gives ONE concrete type for any term count and for static
        # vs time-dependent (the property the custom VecSum buys over SciMLOperators'
        # Tuple-locked AddedOperator).
        @test H2 isa QTB.QuantumObjectEvolution
        @test typeof(H2) === typeof(H3)
        @test typeof(H2) === typeof(Htd)
        # The positional dims form of a single operator is the bare eager QuantumObject
        # (the HilbertSpace / keyword form always wraps in a one-term lazy QobjEvo).
        @test to_numeric(a, 5) isa QTB.QuantumObject
        @test to_numeric(a' * a, h, 5; backend = QTBB) isa QTB.QuantumObjectEvolution
    end

    @testset "Fock parity vs QuantumOptics" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        ψq = QOB.coherentstate(QOB.FockBasis(7), 0.4 + 0.2im)
        ψt = QTB.coherent(8, 0.4 + 0.2im)
        for op in (a, a', a' * a, a + a', 2 * a' * a + a)
            @test numeric_average(op, ψq) ≈ numeric_average(op, ψt) atol = 1.0e-6
        end
    end

    @testset "Transition parity (1-indexed vs 0-indexed)" begin
        h = NLevelSpace(:atom, 3, 1)
        σ(i, j) = Transition(h, :σ, i, j)
        ψq = QOB.nlevelstate(QOB.NLevelBasis(3), 1) +
            2 * QOB.nlevelstate(QOB.NLevelBasis(3), 2) +
            3 * QOB.nlevelstate(QOB.NLevelBasis(3), 3)
        ψq = ψq / norm(ψq.data)
        ψt = QTB.basis(3, 0) + 2 * QTB.basis(3, 1) + 3 * QTB.basis(3, 2)
        ψt = ψt / norm(ψt)
        for i in 1:3, j in 1:3
            @test numeric_average(σ(i, j), ψq) ≈ numeric_average(σ(i, j), ψt) atol = 1.0e-8
        end
    end

    @testset "Pauli / Spin parity" begin
        hp = PauliSpace(:p)
        ψq = (QOB.spindown(QOB.SpinBasis(1 // 2)) + QOB.spinup(QOB.SpinBasis(1 // 2))) / sqrt(2)
        ψt = (QTB.basis(2, 1) + QTB.basis(2, 0)) / sqrt(2)
        for ax in 1:3
            @test numeric_average(Pauli(hp, :σ, ax), ψq) ≈ numeric_average(Pauli(hp, :σ, ax), ψt) atol = 1.0e-8
        end

        hs = SpinSpace(:s)
        for s in (1 // 2, 1)
            d = Int(2s + 1)
            bq = QOB.SpinBasis(s)
            ψq = QOB.spinup(bq)
            ψt = QTB.basis(d, 0)
            for ax in 1:3
                @test numeric_average(Spin(hs, :S, ax), ψq) ≈ numeric_average(Spin(hs, :S, ax), ψt) atol = 1.0e-8
            end
        end
    end

    @testset "Composite parity" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        a = Destroy(h, :a, 1)
        σ(i, j) = Transition(h, :σ, i, j, 2)
        ba = QOB.NLevelBasis(3)
        ψq = QOB.coherentstate(QOB.FockBasis(4), 0.3) ⊗
            ((QOB.nlevelstate(ba, 1) + QOB.nlevelstate(ba, 2)) / sqrt(2))
        ψt = QTB.tensor(QTB.coherent(5, 0.3), (QTB.basis(3, 0) + QTB.basis(3, 1)) / sqrt(2))
        for op in (a' * a, a' * σ(1, 2) + a * σ(2, 1), a' * a * σ(2, 2))
            @test numeric_average(op, ψq) ≈ numeric_average(op, ψt) atol = 1.0e-6
        end
    end

    @testset "Time-dependent (native QobjEvo)" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        @variables E::Number
        N = 6
        f = to_numeric(E * a + conj(E) * a', h, N; backend = QTBB, time_parameter = Dict(E => t -> 1 + im * t))
        @test f isa QTB.QuantumObjectEvolution
        ψt = QTB.coherent(N, 0.3 - 0.1im)
        aT = QTB.destroy(N)
        @test QTB.expect(f(nothing, 0.5), ψt) ≈ QTB.expect((1 + 0.5im) * aT + (1 - 0.5im) * aT', ψt) atol = 1.0e-8
    end

    @testset "sesolve / mesolve consume the lazy VecSum" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        N = 6
        H = to_numeric(2.0 * a' * a + 0.5 * (a + a'), h, N; backend = QTBB)
        Hc = QTB.Qobj(SO.concretize(H.data); dims = N)   # eager reference
        ψ = QTB.coherent(N, 0.3 + 0.1im)
        tl = collect(range(0, 1.0, length = 11))
        eop = QTB.create(N) * QTB.destroy(N)

        rs = QTB.sesolve(H, ψ, tl; e_ops = [eop], progress_bar = Val(false))
        rsc = QTB.sesolve(Hc, ψ, tl; e_ops = [eop], progress_bar = Val(false))
        @test maximum(abs, rs.expect[1, :] .- rsc.expect[1, :]) < 1.0e-6

        c_ops = [QTB.destroy(N)]
        rm = QTB.mesolve(H, ψ, tl, c_ops; e_ops = [eop], progress_bar = Val(false))
        rmc = QTB.mesolve(Hc, ψ, tl, c_ops; e_ops = [eop], progress_bar = Val(false))
        @test maximum(abs, rm.expect[1, :] .- rmc.expect[1, :]) < 1.0e-6
        @test abs(tr(rm.states[end]) - 1) < 1.0e-6
    end

    @testset "numeric_average: averages, scalars, vector of states" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        α = 0.4 + 0.2im
        ψt = QTB.coherent(8, α)
        @test numeric_average(average(a), ψt) ≈ α atol = 1.0e-4
        @test numeric_average(average(a) + average(a'), ψt) ≈ α + conj(α) atol = 1.0e-4
        @test numeric_average(3, ψt) === ComplexF64(3)

        ψs = [QTB.coherent(8, β) for β in (0.1 + 0.2im, -0.2 + 0.1im)]
        @test numeric_average(a, ψs) ≈ [QTB.expect(QTB.destroy(8), ψ) for ψ in ψs] atol = 1.0e-8
    end

    @testset "materialize on demand" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        H = to_numeric(a' * a, h, 4; backend = QTBB)
        M = SO.concretize(H.data)
        @test M ≈ QTB.create(4).data * QTB.destroy(4).data
    end

    @testset "empty QAdd assembles to zero" begin
        # A `QAdd` that cancels to zero terms still assembles to a valid zero operator on
        # the full space (the empty-terms branch of `numeric_assemble`).
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        z = 1 * a - 1 * a
        M = to_numeric(z, h, 5; backend = QTBB)
        @test SO.concretize(M.data) ≈ zeros(ComplexF64, 5, 5)
    end

    @testset "time-dependent with constant term" begin
        # A TD assembly mixing a constant term (`Δ*a'a`) with a genuinely time-dependent
        # term (`f*a`) exercises both `_td_scalar` branches.
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        @variables Δ f
        H = to_numeric(
            Δ * a' * a + f * a, h, 5;
            backend = QTBB, parameter = Dict(Δ => 1.0), time_parameter = Dict(f => t -> 2.0 + 0im),
        )
        M = H(nothing, 0.0).data
        Href = QTB.create(5).data * QTB.destroy(5).data + 2.0 * QTB.destroy(5).data
        @test M ≈ Href
    end

    @testset "unsupported operator role errors" begin
        # A collective transition has no QuantumToolbox realisation: the closed operator
        # ladder falls through to the throwing `Val` extension point.
        hc = CollectiveNLevelSpace(:atom, 3)
        S12 = CollectiveTransition(hc, :S, 1, 2)
        @test_throws ArgumentError to_numeric(S12, 3)
    end

    @testset "default backend errors when both are loaded" begin
        # With QuantumOpticsBase and QuantumToolbox both loaded, the `HilbertSpace` form
        # cannot pick a backend and must ask for an explicit `backend =`.
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        @test_throws ArgumentError to_numeric(a' * a, h, 5)
    end
end
