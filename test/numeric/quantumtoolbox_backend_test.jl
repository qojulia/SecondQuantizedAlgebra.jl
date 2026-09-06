# QuantumToolbox backend for `to_numeric`/`numeric_average`.
#
# `QuantumToolbox` is imported (not `using`) so its `⊗`/`tensor`/`expect` do not collide
# with the shared TensorCore ones that `SecondQuantizedAlgebra` extends; `⊗` here is the
# shared TensorCore operator (works on Hilbert spaces and QuantumOptics bases). Parity
# is checked through `expect`/`numeric_average` scalars (never raw matrices, whose kron
# convention differs across backends), with states built per-backend in slot order.

using SecondQuantizedAlgebra
import QuantumToolbox as QTB
import QuantumOpticsBase as QOB
import SciMLOperators as SO
using Symbolics: @variables
import LinearAlgebra: tr, norm, diag
using Test

const QTBB = QuantumToolboxBackend()

function inferred_qtb_uniform(op::Op, h::FockSpace, n::Int)
    return to_numeric(op, h, n; backend = QTBB)
end

function inferred_qtb_product(op::Op, h::ProductSpace, dims::NTuple{2, Int})
    return to_numeric(op, h, dims; backend = QTBB)
end

function inferred_qtb_raw_product(op::Op, dims::NTuple{2, Int})
    return to_numeric(op, dims)
end

@testset "QuantumToolbox backend" begin
    @testset "eager default; lazy assembly is n-stable" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        H2 = to_numeric(a' * a, h, 5; backend = QTBB, op_type = identity)
        H3 = to_numeric(a' * a + a + a', h, 5; backend = QTBB, op_type = identity)
        @variables E::Number
        Htd = to_numeric(E * a + conj(E) * a', h, 5; backend = QTBB, time_parameter = Dict(E => t -> 1 + im * t))

        # The vector-backed VecSum gives ONE concrete type for any term count and for static
        # vs time-dependent (the property the custom VecSum buys over SciMLOperators'
        # Tuple-locked AddedOperator).
        @test H2 isa QTB.QuantumObjectEvolution
        @test typeof(H2) === typeof(H3)
        @test typeof(H2) === typeof(Htd)
        # Static defaults are eager and sparse on both direct and uniform forms; lazy is
        # an explicit representation choice.
        @test to_numeric(a, 5) isa QTB.QuantumObject
        @test typeof(to_numeric(a, 5).data) === typeof(QTB.destroy(5).data)
        @test to_numeric(a' * a, h, 5; backend = QTBB) isa QTB.QuantumObject
        @test typeof(to_numeric(a' * a, h, 5; backend = QTBB).data) ===
            typeof(QTB.destroy(6).data)
    end

    @testset "concrete materialization and expectations" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        ψ = QTB.coherent(8, 0.4 + 0.2im)
        bra = ψ'
        ρ = QTB.ket2dm(ψ)
        n = a' * a

        @test @inferred(to_numeric(a, 8)) isa QTB.QuantumObject
        @test @inferred(to_numeric(n, 8)) isa QTB.QuantumObject
        @test @inferred(inferred_qtb_uniform(a, h, 7)) isa QTB.QuantumObject

        hp = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        ap = Destroy(hp, :a, 1)
        @test @inferred(inferred_qtb_product(ap, hp, (7, 3))) isa QTB.QuantumObject
        @test @inferred(inferred_qtb_raw_product(ap, (8, 3))) isa QTB.QuantumObject

        lazy = to_numeric(n, 8; op_type = identity)
        @test @inferred(SO.concretize(lazy.data)) isa AbstractMatrix
        @test @inferred(QTB.expect(lazy, ψ)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), ψ)
        @test @inferred(QTB.expect(lazy, bra)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), bra)
        @test @inferred(QTB.expect(lazy, ρ)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), ρ)
        @test numeric_basis(bra) == 8
        @test @inferred(numeric_average(a, ψ)) ≈ QTB.expect(QTB.destroy(8), ψ)
        @test @inferred(numeric_average(n, ψ)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), ψ)
        @test @inferred(numeric_average(n, bra)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), bra)
        @test @inferred(numeric_average(n, ρ)) ≈ QTB.expect(QTB.create(8) * QTB.destroy(8), ρ)

        qop = to_numeric(n, 8)
        @test numeric_expect(QTBB, qop, ψ) ≈ QTB.expect(qop, ψ)
        @test numeric_expect(QTBB, qop, bra) ≈ QTB.expect(qop, bra)
        @test numeric_expect(QTBB, qop, ρ) ≈ QTB.expect(qop, ρ)

        dense_qop = QTB.Qobj(Matrix(qop.data), QTB.Operator(), 8)
        @test numeric_expect(QTBB, dense_qop, ρ) ≈ QTB.expect(dense_qop, ρ)

        # Materialization accepts custom dense leaves but preserves the sparse public form.
        dense_a = QTB.Qobj(Matrix(QTB.destroy(8).data), QTB.Operator(), 8)
        dense_lazy = to_numeric(a + a', 8; operators = Dict(a => dense_a), op_type = identity)
        dense_data = @inferred SO.concretize(dense_lazy.data)
        @test typeof(dense_data) === typeof(QTB.destroy(8).data)
        @test dense_data ≈ QTB.destroy(8).data + QTB.create(8).data
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
        ψt = QTB.coherent(N + 1, 0.3 - 0.1im)
        aT = QTB.destroy(N + 1)
        @test QTB.expect(f(nothing, 0.5), ψt) ≈ QTB.expect((1 + 0.5im) * aT + (1 - 0.5im) * aT', ψt) atol = 1.0e-8

        @test to_numeric(
            E * a, h, N;
            backend = QTBB, time_parameter = Dict(E => t -> 1), op_type = identity,
        ) isa QTB.QuantumObjectEvolution
        @test_throws ArgumentError to_numeric(
            E * a, h, N;
            backend = QTBB, time_parameter = Dict(E => t -> 1), op_type = QTB.to_dense,
        )
    end

    @testset "sesolve / mesolve consume the lazy VecSum" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        N = 6
        H = to_numeric(
            2.0 * a' * a + 0.5 * (a + a'), h, N;
            backend = QTBB, op_type = identity,
        )
        Hc = QTB.Qobj(SO.concretize(H.data); dims = N + 1)   # eager reference
        ψ = QTB.coherent(N + 1, 0.3 + 0.1im)
        tl = collect(range(0, 1.0, length = 11))
        eop = QTB.create(N + 1) * QTB.destroy(N + 1)

        rs = QTB.sesolve(H, ψ, tl; e_ops = [eop], progress_bar = Val(false))
        rsc = QTB.sesolve(Hc, ψ, tl; e_ops = [eop], progress_bar = Val(false))
        @test maximum(abs, rs.expect[1, :] .- rsc.expect[1, :]) < 1.0e-6

        c_ops = [QTB.destroy(N + 1)]
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
        H = to_numeric(a' * a, h, 4; backend = QTBB, op_type = identity)
        M = SO.concretize(H.data)
        @test M ≈ QTB.create(5).data * QTB.destroy(5).data
    end

    @testset "empty QAdd assembles to zero" begin
        # A `QAdd` that cancels to zero terms still assembles to a valid zero operator on
        # the full space (the empty-terms branch of `numeric_assemble`).
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        z = 1 * a - 1 * a
        M = to_numeric(z, h, 5; backend = QTBB, op_type = identity)
        @test SO.concretize(M.data) ≈ zeros(ComplexF64, 6, 6)
    end

    @testset "time-dependent with constant term" begin
        # A TD assembly mixing a constant term (`Δ*a'a`) with a genuinely time-dependent
        # term (`f*a`) exercises both `td_scalar` branches.
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        @variables Δ f
        H = to_numeric(
            Δ * a' * a + f * a, h, 5;
            backend = QTBB, parameter = Dict(Δ => 1.0), time_parameter = Dict(f => t -> 2.0 + 0im),
        )
        M = H(nothing, 0.0).data
        Href = QTB.create(6).data * QTB.destroy(6).data + 2.0 * QTB.destroy(6).data
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

    @testset "Indexed path works for QuantumToolbox" begin
        # Sum of excited-state projectors over two atoms is the excited-atom-count operator.
        ha = NLevelSpace(:atom, 2)
        i = Index(ha, :i, 2, ha)
        σ(a, b, k) = IndexedOperator(Transition(ha, :σ, a, b), k)
        H = Σ(σ(2, 2, i), i)
        sites = Dict{Int, Vector{Int}}(1 => [1, 2])

        Mt = to_numeric(H, [2, 2], sites)               # QTB (Vector{Int} dims)
        Mt_tuple = to_numeric(H, (2, 2), sites)         # tuple dims use the same public path
        @test Mt isa QTB.QuantumObject
        @test Mt_tuple.data ≈ Mt.data

        P = QTB.projection(2, 1, 1)                     # |1><1|, 0-indexed excited level
        ref = QTB.tensor(P, QTB.qeye(2)) + QTB.tensor(QTB.qeye(2), P)
        @test Mt.data ≈ ref.data

        # Cross-check the eager expectation via a scalar (matrices differ by kron convention
        # across backends, so compare scalars, not raw matrices).
        ψt = QTB.tensor(QTB.basis(2, 1), QTB.basis(2, 1))   # both excited -> eigenvalue 2
        @test real(QTB.expect(Mt, ψt)) ≈ 2 atol = 1.0e-10

        @test_throws ArgumentError to_numeric(H, 2, sites)
        @test_throws ArgumentError to_numeric(H, [2, 2], Dict(1 => [1, 1]))
    end

    @testset "Fock/PhaseSpace dims match QuantumOptics" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        # cutoff 5 gives dimension 6 (occupations 0:5) on BOTH backends.
        @test SecondQuantizedAlgebra.numeric_basis(QuantumToolboxBackend(), FockSpace(:f), 5) == 6
        @test SecondQuantizedAlgebra.numeric_basis(QuantumToolboxBackend(), PhaseSpace(:x), 5) == 6

        Ht = to_numeric(a' * a, h, 5; backend = QuantumToolboxBackend())
        @test size(Ht.data) == (6, 6)
        # highest occupation is 5 (would be 4 under the old off-by-one convention).
        @test maximum(real, diag(Matrix(Ht.data))) ≈ 5

        # Direct QTB dimensions are raw matrix dimensions; the HilbertSpace form receives
        # a symbolic Fock cutoff.
        @test size(to_numeric(a, 5).data) == (5, 5)
        @test size(to_numeric(a, h, 5; backend = QTBB).data) == (6, 6)

        # Bare QAdd on raw QTB dimensions (positional QAdd-on-dims form).
        qadd_raw = to_numeric(a' + a, 5)
        @test size(qadd_raw.data) == (5, 5)
        @test Matrix(qadd_raw.data) ≈ Matrix((to_numeric(a', 5) + to_numeric(a, 5)).data)
    end

    @testset "explicit op_type materializes eager" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)
        ref = QTB.create(6).data * QTB.destroy(6).data     # cutoff 5 gives dim 6

        sp = to_numeric(a' * a, h, 5; backend = QTBB, op_type = QTB.to_sparse)
        @test sp isa QTB.QuantumObject
        @test sp.data ≈ ref

        de = to_numeric(a' * a, h, 5; backend = QTBB, op_type = QTB.to_dense)
        @test de isa QTB.QuantumObject
        @test de.data ≈ ref

        co = to_numeric(a' * a, h, 5; backend = QTBB, op_type = SO.concretize)
        @test co isa AbstractMatrix
        @test co ≈ ref

        # identity keeps the lazy assembly.
        lz = to_numeric(a' * a, h, 5; backend = QTBB, op_type = identity)
        @test lz isa QTB.QuantumObjectEvolution
    end

    @testset "vector API for QuantumToolbox" begin
        h = FockSpace(:f)
        @qnumbers a::Destroy(h)

        # HilbertSpace form (backend-neutral entry point).
        vhs = to_numeric([a, a'], h, 5; backend = QTBB)
        @test vhs isa AbstractVector
        @test length(vhs) == 2
        @test all(x -> x isa QTB.QuantumObject, vhs)

        # Direct QTB-dims form.
        vd = to_numeric([a, a'], 5)
        @test vd isa AbstractVector
        @test length(vd) == 2
        @test all(x -> x isa QTB.QuantumObject, vd)
    end

    @testset "public backend hooks and product dims validation" begin
        ψ = QTB.basis(3, 0)
        @test numeric_backend(ψ) isa QuantumToolboxBackend
        @test numeric_basis(ψ) == 3
        @test SecondQuantizedAlgebra.numeric_num_subsystems(QTBB, 3) == 1
        @test SecondQuantizedAlgebra.numeric_num_subsystems(QTBB, [2, 3]) == 2

        h = FockSpace(:a) ⊗ FockSpace(:b)
        a = Destroy(h, :a, 1)
        @test_throws ArgumentError to_numeric(a, h, [2]; backend = QTBB)
        @test_throws ArgumentError to_numeric(a, h, [2, 3, 99]; backend = QTBB)
        # dims that are not an indexable collection are rejected up front.
        @test_throws ArgumentError to_numeric(a, h, nothing; backend = QTBB)

        # A valid ProductSpace HilbertSpace conversion builds per-subspace dims (Fock
        # cutoffs 2 and 3 give dimensions 3 and 4, so the composite is 12×12).
        Mp = to_numeric(a, h, [2, 3]; backend = QTBB)
        @test size(Mp.data) == (12, 12)
    end

    @testset "p-aware time_parameter (differentiable control)" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        @variables Δ ε κ
        N = 5
        nop = to_numeric(a' * a, h, N; backend = QTBB)
        ψ0 = QTB.fock(N + 1, 0)
        tl = collect(range(0, 2.0, length = 21))

        # A two-argument (p, t) value function makes the QobjEvo p-aware.
        H = to_numeric(
            Δ * a' * a + ε * (a' + a), h, N;
            backend = QTBB,
            time_parameter = Dict(Δ => (p, t) -> p[1], ε => (p, t) -> p[2] * cos(t)),
        )
        @test H isa QTB.QuantumObjectEvolution

        finaln(pv) = real(QTB.expect(nop, QTB.sesolve(H, ψ0, tl; params = pv, progress_bar = Val(false)).states[end]))
        n1 = finaln([1.0, 0.5])
        n2 = finaln([1.0, 0.6])
        # p actually flows to the coefficient (forward finite-difference sanity, no AD dep).
        @test !isapprox(n1, n2; atol = 1.0e-8)

        # Equivalence with a hand-assembled (p, t) QobjEvo.
        adaga = to_numeric(a' * a, h, N; backend = QTBB)
        xop = to_numeric(a' + a, h, N; backend = QTBB)
        Hhand = QTB.QobjEvo(((adaga, (p, t) -> p[1]), (xop, (p, t) -> p[2] * cos(t))))
        nhand = real(QTB.expect(nop, QTB.sesolve(Hhand, ψ0, tl; params = [1.0, 0.5], progress_bar = Val(false)).states[end]))
        @test n1 ≈ nhand atol = 1.0e-8

        # Mixed arity: one t-only and one (p, t) symbol; the t-only symbol is lifted.
        Hmix = to_numeric(
            Δ * a' * a + ε * (a' + a), h, N;
            backend = QTBB,
            time_parameter = Dict(Δ => t -> 1.0, ε => (p, t) -> p[1] * cos(t)),
        )
        @test Hmix.data isa SO.AbstractSciMLOperator
        smix = QTB.sesolve(Hmix, ψ0, tl; params = [0.5], progress_bar = Val(false))
        @test real(QTB.expect(nop, smix.states[end])) isa Float64

        # mesolve superoperator path with a p-aware collapse operator: no lazy-tensor warning,
        # a VecSum superoperator, and κ (from p) flows into the dissipation.
        c1 = to_numeric(sqrt(κ) * a, h, N; backend = QTBB, time_parameter = Dict(κ => (p, t) -> p[1]))
        @test c1.data isa SO.AbstractSciMLOperator
        L = QTB.liouvillian(to_numeric(2.0 * a' * a, h, N; backend = QTBB), [c1])
        ρ0 = QTB.ket2dm(QTB.fock(N + 1, 2))
        tl2 = collect(range(0, 1.5, length = 16))
        fme(pv) = real(QTB.expect(nop, QTB.mesolve(L, ρ0, tl2; params = pv, progress_bar = Val(false)).states[end]))
        @test !isapprox(fme([0.3]), fme([0.6]); atol = 1.0e-8)

        # QuantumToolbox p-aware collapse decays faster with a larger rate.
        @test fme([0.6]) < fme([0.3])

        # Value-function arity is validated with actionable errors, not silently misdetected:
        # variadic, conflicting-arity, and out-of-range (0 or 3+) functions are rejected.
        bad_va(args...) = args[1]
        bad_multi(t) = t
        bad_multi(p, t) = p[1]
        bad_zero() = 1.0
        bad_three(x, y, z) = x
        for badf in (bad_va, bad_multi, bad_zero, bad_three)
            @test_throws ArgumentError to_numeric(ε * a, h, N; backend = QTBB, time_parameter = Dict(ε => badf))
        end
    end

end
