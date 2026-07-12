using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QAdd, QSym, _single_qadd, _to_cnum, _to_complex, _fold_const
using QuantumOpticsBase
using Symbolics: @variables, substitute
import SymbolicUtils
using Test

@testset "numeric conversion" begin
    @testset "Single space — basic" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
    end

    @testset "Products" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a' * a, b) == create(b) * destroy(b)
        @test to_numeric(2 * a, b) == 2 * destroy(b)
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        expected = destroy(b) + create(b)
        @test result == expected
    end

    @testset "Scalar" begin
        b = FockBasis(7)
        @test to_numeric(3, b) == 3 * one(b)
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)

        # Default materialises to a concrete sparse operator; the lazy
        # representation is opt-in via `op_type = identity`.
        a1_lazy = to_numeric(a1, b; op_type = identity)
        @test a1_lazy isa LazyTensor
        a1_num = to_numeric(a1, b)
        @test a1_num isa Operator
        @test a1_num == sparse(a1_lazy)
        @test to_numeric(a1, b; op_type = dense) == dense(a1_lazy)
    end

    @testset "numeric_average" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end

    @testset "NLevel numeric" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        b = NLevelBasis(3)
        @test to_numeric(σ12, b) == transition(b, 1, 2)
        @test to_numeric(σ12', b) == transition(b, 2, 1)
    end

    @testset "Pauli numeric" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)
        b = SpinBasis(1 // 2)
        @test to_numeric(σx, b) == sigmax(b)
        @test to_numeric(σy, b) == sigmay(b)
        @test to_numeric(σz, b) == sigmaz(b)
    end

    @testset "Spin numeric" begin
        h = SpinSpace(:s)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        b = SpinBasis(5 // 2)
        @test to_numeric(Sx, b) == 0.5 * sigmax(b)
        @test to_numeric(Sy, b) == 0.5 * sigmay(b)
        @test to_numeric(Sz, b) == 0.5 * sigmaz(b)
    end

    @testset "Composite NLevel + Fock" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        @qnumbers a::Destroy(h, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        bf = FockBasis(3)
        bn = NLevelBasis(3)
        bc = bf ⊗ bn
        @test to_numeric(σ12, bc; op_type = identity) isa LazyTensor
        @test to_numeric(σ12, bc) isa Operator
    end

    @testset "Type stability" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        ψ = fockstate(b, 2)

        # to_numeric: leaf + QAdd + dict-substitution paths
        @test @inferred(to_numeric(a, b)) isa AbstractOperator
        @test @inferred(to_numeric(a' * a, b)) isa AbstractOperator
        @test @inferred(to_numeric(2 * a + 3 * a' + 5, b)) isa AbstractOperator
        @test @inferred(to_numeric(a, b, Dict{QSym, Union{}}())) isa AbstractOperator

        # numeric_average: every BasicSymbolic branch infers to ComplexF64
        @test @inferred(numeric_average(a' * a, ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a) + average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(2 * average(a) * average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a)^2, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3.5 + 1im, ψ)) isa ComplexF64

        # `_to_complex(::Any) → ComplexF64` is the keystone; canary against a
        # 7th overload tripping Julia's union-split budget.
        @test Base.return_types(_to_complex, (Any,))[1] === ComplexF64
    end

    @testset "numeric_average: Average expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        # numeric_average on average(a) should give ⟨a⟩ = α
        avg_a = average(a)
        @test numeric_average(avg_a, ψ) ≈ α

        # Sum of averages
        avg_sum = average(a) + average(a')
        @test numeric_average(avg_sum, ψ) ≈ α + conj(α)

        # Scalar times average
        avg_scaled = 2 * average(a)
        @test numeric_average(avg_scaled, ψ) ≈ 2α
    end

    @testset "to_numeric — Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        custom_op = 2 * destroy(b)
        d = Dict(a => custom_op)

        # Dict substitution replaces operator
        @test to_numeric(a, b, d) == custom_op
        # Adjoint not in dict → falls back to normal
        @test to_numeric(a', b, d) == create(b)

        # QAdd with Dict
        result_mul = to_numeric(a' * a, b, d)
        expected_mul = create(b) * custom_op
        @test result_mul == expected_mul
    end

    @testset "numeric_average — Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        d = Dict{QSym, Any}()  # empty dict — same as no dict
        @test numeric_average(a, ψ, d) ≈ α

        # Average expression with dict
        avg_a = average(a)
        @test numeric_average(avg_a, ψ, d) ≈ α

        # Number passthrough — coerced to ComplexF64 like every other branch
        @test numeric_average(3, ψ, d) === ComplexF64(3)
    end

    @testset "Composite basis with gaps" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3, 1)
        hprod_gap = hfock ⊗ hnlevel ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod_gap = bfock ⊗ bnlevel ⊗ bnlevel

        a = Destroy(hprod_gap, :a, 1)
        σprod_gap(i, j) = Transition(hprod_gap, :σ, i, j, 3)

        for i in 1:3, j in 1:3
            i == j == 1 && continue
            op1 = a * σprod_gap(i, j)
            op2 = a' * σprod_gap(i, j)
            @test to_numeric(op1, bprod_gap; op_type = identity) == LazyTensor(
                bprod_gap,
                [1, 3],
                (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            @test to_numeric(op2, bprod_gap; op_type = identity) == LazyTensor(
                bprod_gap,
                [1, 3],
                (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            # Default returns the same operator, materialised sparse.
            @test to_numeric(op1, bprod_gap) == sparse(
                to_numeric(op1, bprod_gap; op_type = identity),
            )
        end
    end

    @testset "Large Hilbert space" begin
        hfock = FockSpace(:fock)
        @qnumbers a::Destroy(hfock)
        bfock = FockBasis(100)

        @test isequal(
            2 * create(bfock) + 2 * destroy(bfock),
            to_numeric(2 * a + 2 * a', bfock),
        )
        @test iszero(
            (2 * create(bfock) + 2 * destroy(bfock)) -
                to_numeric(2 * a + 2 * a', bfock),
        )
        @test isequal(to_numeric(2 * a, bfock), 2 * to_numeric(a, bfock))
    end

    @testset "numeric_average — product state" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3, 1)
        hprod = hfock ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel

        α = 0.1 + 0.2im
        ψ = coherentstate(bfock, α)
        ψprod = ψ ⊗ nlevelstate(bnlevel, 1)

        σprod(i, j) = Transition(hprod, :σ, i, j, 2)

        idfock = one(bfock)
        for i in 1:3, j in 1:3
            op = σprod(i, j)
            op_num = idfock ⊗ QuantumOpticsBase.transition(bnlevel, i, j)
            @test numeric_average(op, ψprod) ≈ expect(op_num, ψprod)
        end
    end

    @testset "numeric_average — LazyKet state" begin
        # `numeric_average` must contract against a lazy operator here: `expect`
        # has no method for a sparse operator paired with a `LazyKet`, and a
        # `LazyKet` exists precisely to avoid materialising the full state.
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, (:a, :b, :c))
        hprod = hfock ⊗ hnlevel
        bfock = FockBasis(10)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel

        @qnumbers a::Destroy(hprod, 1)
        σ(i, j) = Transition(hprod, :σ, i, j, 2)

        α = 0.3 + 0.0im
        ket_fock = coherentstate(bfock, α)
        ket_nlevel = (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2)
        ψ_lazy = LazyKet(bprod, (ket_fock, ket_nlevel))
        ψ_dense = ket_fock ⊗ ket_nlevel

        for op in (a, a' * σ(:a, :c), a + a' * σ(:a, :c))
            op_num = to_numeric(op, bprod)
            @test numeric_average(op, ψ_lazy) ≈ expect(op_num, ψ_dense)
        end
        # average-expression entry point (undo_average path) on a LazyKet
        @test numeric_average(average(a' * a), ψ_lazy) ≈
            expect(to_numeric(a' * a, bprod), ψ_dense)
    end

    @testset "numeric_average — comprehensive expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a + a'a, ψ) ≈ α + abs2(α)
        @test numeric_average(average(a) + average(a'a), ψ) ≈ α + abs2(α)
        @test numeric_average(average(a + a'a), ψ) ≈ α + abs2(α)
        @test numeric_average(average(a) * average(a'a), ψ) ≈ α * α' * α
        @test numeric_average(average(a)^2, ψ) ≈ α^2
        @test numeric_average(3, ψ) ≈ 3
    end

    @testset "numeric_average — vector of states" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im, 0.5 + 0.0im)
        ψs = [coherentstate(b, α) for α in αs]

        # Vector input → vector output, same as broadcasting the scalar form
        @test numeric_average(a, ψs) ≈ [α for α in αs]
        @test numeric_average(a' * a, ψs) ≈ [abs2(α) for α in αs]
        @test numeric_average(average(a), ψs) ≈ [α for α in αs]

        # expect alias matches
        @test expect(a, ψs) ≈ numeric_average(a, ψs)
        @test expect(a' * a, ψs[1]) ≈ numeric_average(a' * a, ψs[1])
        @test numeric_average(average(a), ψs[1]) ≈ αs[1]

        # Empty vector errors
        empty_ψs = typeof(ψs[1])[]
        @test_throws ArgumentError numeric_average(a, empty_ψs)
    end

    @testset "numeric_average — vector of states with Dict" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im)
        ψs = [coherentstate(b, α) for α in αs]
        d = Dict{QSym, Any}()       # empty dict — passthrough to base case

        @test numeric_average(a, ψs, d) ≈ numeric_average(a, ψs)
        @test numeric_average(average(a' * a), ψs, d) ≈ [abs2(α) for α in αs]
        @test expect(a, ψs, d) ≈ numeric_average(a, ψs, d)
    end

    @testset "numeric_average — Dict comprehensive" begin
        nQDs = 2
        h_qc1 = FockSpace(:ada)
        h_qc2 = FockSpace(:n)
        h_qc = h_qc1 ⊗ h_qc2
        a = Destroy(h_qc, :a, 1)
        n = Destroy(h_qc, :n, 2)
        ad = a'

        bs = NLevelBasis(2)
        b_all = tensor([bs for i in 1:nQDs]...)
        s(α, i, j) = embed(b_all, α, transition(bs, i, j))
        b_test = FockBasis(2) ⊗ FockBasis(3)

        dd = Dict([ad, a] .=> [s(2, 2, 1), s(2, 1, 2)])

        @test to_numeric(a, b_test, Dict{QSym, Any}()) == to_numeric(a, b_test)
        @test to_numeric(ad * n, b_test, Dict{QSym, Any}()) == to_numeric(ad * n, b_test)
        @test to_numeric(2 * ad * a, b_test, Dict{QSym, Any}()) == to_numeric(2 * ad * a, b_test)
        @test to_numeric(ad, b_all, dd) == s(2, 2, 1)
        @test to_numeric(2 * ad * a, b_all, dd) == 2 * s(2, 2, 1) * s(2, 1, 2)
        @test dense(to_numeric(3, b_all, dd)) == one(b_all) * 3

        ψ0 = tensor([nlevelstate(bs, 2) for i in 1:nQDs]...)
        @test numeric_average(average(ad * a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0)
        @test numeric_average(average(ad) * average(ad * a) + 3, ψ0, dd) ==
            expect(s(2, 2, 1), ψ0) * expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + 3
        @test numeric_average(3 * average(ad)^2, ψ0, dd) ==
            3 * expect(s(2, 2, 1), ψ0)^2
        @test numeric_average(average(ad * a) + average(a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + expect(s(2, 1, 2), ψ0)
    end

    @testset "Allocations — to_numeric" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        # Warmup
        to_numeric(a, b)
        to_numeric(a', b)
        to_numeric(a' * a, b)

        # Single operators should be bounded
        @test @allocations(to_numeric(a, b)) < 50
        @test @allocations(to_numeric(a', b)) < 50
        @test @allocations(to_numeric(a' * a, b)) < 1500
    end

    @testset "Round-trip with coherent state" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
        @test to_numeric(a' * a, b) == create(b) * destroy(b)

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)

        # numeric_average on QAdd
        expr = a + a' * a
        @test numeric_average(expr, ψ) ≈ α + abs2(α)

        # to_numeric with scalar QAdd (empty operators)
        @test to_numeric(_single_qadd(_to_cnum(3), Op[]), b) == 3 * one(b)
    end

    @testset "Number-symtype coefficient round-trip" begin
        @variables g::Number
        b = FockBasis(3)
        a = Destroy(FockSpace(:c), :a)

        # g†g must reduce to |g|², i.e. conj(g)*g, after substituting a complex value.
        for v in (2 + 3im, 1 + 1im, 0.5 - 0.25im)
            op = substitute((g * a)' * (g * a), Dict(g => v))
            @test to_numeric(op, b) ≈ abs2(v) * to_numeric(a' * a, b)
        end

        # A bare complex coupling carries through (no conj): g·a → value·a.
        @test to_numeric(substitute(g * a, Dict(g => 2 + 3im)), b) ≈ (2 + 3im) * to_numeric(a, b)
    end

    @testset "Public coefficient lowering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        @variables x::Real

        coeff = first(x * a).second
        lowered = SecondQuantizedAlgebra.to_num(coeff)
        @test isequal(real(lowered), x)
        @test iszero(imag(lowered))
    end

    @testset "to_numeric keyword translation" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        Ad = create(b)
        @variables x::Real ϕ::Real E::Number

        @test to_numeric(substitute(sqrt(x) * a, Dict(x => 2.0)), b) ≈ sqrt(2.0) * A
        @test to_numeric(substitute(exp(im * ϕ) * a, Dict(ϕ => 0.5)), b) ≈ exp(0.5im) * A

        @test to_numeric(sqrt(x) * a, b; parameter = Dict(x => 2.0)) ≈ sqrt(2.0) * A
        @test to_numeric(exp(im * ϕ) * a, b; parameter = Dict(ϕ => 0.5)) ≈ exp(0.5im) * A

        custom = 2 * A
        @test to_numeric(a', b; operators = Dict(a => custom)) == custom'
        @test to_numeric(a', b; operators = Dict(a => custom), adjoint_ops = false) == Ad
        @test to_numeric([a, a'], b; operators = Dict(a => custom)) == [custom, custom']

        f = to_numeric(E * a + conj(E) * a', b; time_parameter = Dict(E => t -> 1 + im * t))
        @test f(0.5) ≈ (1 + 0.5im) * A + (1 - 0.5im) * Ad
    end

    @testset "keyword to_numeric on composite basis" begin
        hf = FockSpace(:c) ⊗ FockSpace(:d)
        a = Destroy(hf, :a, 1)
        @variables Δ::Real
        b = FockBasis(2) ⊗ FockBasis(3)
        Ia = identityoperator(b)
        an = to_numeric(a, b)

        # Regression: a scalar/constant term must build the full-system identity, not
        # the identity of a single subsystem (the `_lazy_one(b)` fix). Before, the
        # `+ 5` term produced a wrongly-sized identity on a composite basis.
        H = to_numeric(Δ * a' * a + 5, b; parameter = Dict(Δ => 2.0))
        @test size(dense(H).data) == (length(b), length(b))
        @test dense(H).data ≈ dense(2.0 * an' * an + 5 * Ia).data

        # Custom operator on one subsystem; the constant term is still full-sized.
        custom = 2 * an
        H2 = to_numeric(a + 3, b; operators = Dict(a => custom))
        @test dense(H2).data ≈ dense(custom + 3 * Ia).data

        # A term that cancels to zero still yields a full-sized zero operator.
        z = to_numeric(a - a, b)
        @test size(dense(z).data) == (length(b), length(b))
        @test dense(z).data ≈ zeros(length(b), length(b))

        # Missing adjoint rule is auto-added for custom operators.
        r = to_numeric(a', b; operators = Dict(a => custom))
        @test dense(r).data ≈ dense(custom').data
        r2 = to_numeric(a', b; operators = Dict(a => custom), adjoint_ops = false)
        @test dense(r2).data ≈ dense(to_numeric(a', b)).data
    end

    @testset "time_parameter variants" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        @variables E::Number

        # A plain number value becomes a constant-in-time closure.
        f0 = to_numeric(E * a, b; time_parameter = Dict(E => 3.0))
        @test f0(0.0) ≈ 3.0 * A
        @test f0(10.0) ≈ 3.0 * A

        # A constant-coefficient term still returns a callable when time_parameter is set.
        f1 = to_numeric(2.0 * a, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test f1(7.0) ≈ 2.0 * A

        # `conj(v)` is accepted as a time_parameter key.
        f2 = to_numeric(conj(E) * a, b; time_parameter = Dict(conj(E) => t -> 2 + im * t))
        @test f2(1.0) ≈ (2 + 1im) * A
    end

    @testset "keyword to_numeric, scalar argument" begin
        b = FockBasis(4)
        Ib = one(b)
        @variables x::Real E::Number

        # Plain number scalar routes through the keyword path to a scaled identity.
        @test to_numeric(3, b; parameter = Dict(x => 1.0)) == 3 * Ib
        @test to_numeric(2.0 + 0im, b; parameter = Dict{Any, Any}()) == (2.0 + 0im) * Ib

        # A symbolic scalar resolved by `parameter` becomes a constant identity.
        @test to_numeric(2 * x, b; parameter = Dict(x => 1.5)) ≈ 3.0 * Ib

        # A symbolic scalar with no value cannot be translated.
        @test_throws ArgumentError to_numeric(x, b)

        # With `time_parameter`, a constant scalar yields a constant-in-time closure.
        fconst = to_numeric(2.0, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test fconst(0.0) == 2.0 * Ib
        @test fconst(9.0) == 2.0 * Ib

        # A genuinely time-dependent scalar yields a time-varying closure.
        ft = to_numeric(E, b; time_parameter = Dict(E => t -> 1 + im * t))
        @test ft(0.0) ≈ (1 + 0im) * Ib
        @test ft(2.0) ≈ (1 + 2im) * Ib
    end

    @testset "keyword to_numeric, state argument" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @variables Δ::Real

        # The state form derives the basis from the state and forwards keywords.
        op_state = to_numeric(Δ * a, ψ; parameter = Dict(Δ => 2.0))
        @test op_state == 2.0 * destroy(b)
    end

    @testset "constant symbolic coefficient reduction" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(3)
        A = destroy(b)
        @variables z::Number

        # Coefficients that stay symbolic-but-constant (built from `real`/`imag`/`conj`,
        # which are not folded into the native coefficient tier) are reduced to a
        # concrete number when lowered for `to_numeric`. Each shape exercises a distinct
        # arithmetic node in the constant folder.
        c0 = 1.0 + 2.0im
        cases = (
            ("real", real(conj(z)), real(conj(c0))),
            ("imag", imag(conj(z)), imag(conj(c0))),
            ("conj", conj(real(z) + im * imag(z)), conj(c0)),
            ("plus", real(z) + real(conj(z)), real(c0) + real(conj(c0))),
            ("times", real(z) * real(conj(z)), real(c0) * real(conj(c0))),
            ("div", real(z) / real(conj(z)), real(c0) / real(conj(c0))),
            ("pow", real(z)^3, real(c0)^3),
        )
        for (label, coeff, expected) in cases
            op = substitute(coeff * a, Dict(z => c0))
            @test to_numeric(op, b) ≈ expected * A
        end

        # An unrecognized constant function (`sin`) is not handled by the folder and
        # falls back to the compile-based reduction path.
        op_sin = substitute(sin(z) * a, Dict(z => 0.5))
        @test to_numeric(op_sin, b) ≈ sin(0.5) * A

        # The folder also handles raw subtraction nodes. Symbolics canonicalizes
        # subtraction and negation to `+`/`*`, so a literal `-` node never reaches
        # the folder through normal arithmetic; build it directly to cover both arms.
        neg_unary = SymbolicUtils.term(-, 5.0 + 0im; type = Number)
        neg_binary = SymbolicUtils.term(-, 7.0 + 0im, 2.0 + 0im; type = Number)
        @test _fold_const(neg_unary) == -(5.0 + 0im)
        @test _fold_const(neg_binary) == (7.0 + 0im) - (2.0 + 0im)
    end

    @testset "keyword to_numeric: op_type, vector, complex params, errors" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        Ad = create(b)
        @variables x::Real z::Number E::Number

        # op_type transforms each emitted operator.
        Hd = to_numeric(2.0 * a' * a, b; op_type = dense)
        @test Hd isa QuantumOpticsBase.Operator
        @test Hd ≈ dense(2.0 * Ad * A)

        # Vector form forwards keywords to each element.
        @test to_numeric([a, a'], b; parameter = Dict(x => 1.0)) == [A, Ad]

        # A complex-valued parameter key is split into real/imaginary substitutions.
        @test to_numeric(z * a, b; parameter = Dict(z => 2 + 3im)) ≈ (2 + 3im) * A

        # Error paths.
        @test_throws ArgumentError to_numeric(x * a, b)                         # symbolic, no value
        @test_throws ArgumentError to_numeric(a, b; operators = Dict(x => A))   # non-QSym key
        @test_throws ArgumentError to_numeric(
            x * E * a, b; time_parameter = Dict(E => t -> 1.0 + 0im),
        )                                                                       # untimed variable
    end

    @testset "op_type is shape-independent" begin
        # A bare operator, a product, and a sum must all return the same
        # operator type for a given op_type — the representation depends only on
        # op_type, never on the shape of the expression.
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)
        exprs = (a1, a1' * a1, a1' * a1 + a2' * a2)

        for expr in exprs
            @test to_numeric(expr, b) isa Operator                     # sparse default
            @test to_numeric(expr, b) == to_numeric(expr, b; op_type = sparse)
            @test to_numeric(expr, b; op_type = dense) isa Operator
        end
        # The lazy opt-in yields lazy types across every shape.
        @test to_numeric(exprs[1], b; op_type = identity) isa LazyTensor
        @test to_numeric(exprs[3], b; op_type = identity) isa LazySum
        # Values agree across representations.
        for expr in exprs
            @test to_numeric(expr, b) == sparse(to_numeric(expr, b; op_type = identity))
            @test dense(to_numeric(expr, b)) ≈ dense(to_numeric(expr, b; op_type = dense))
        end
    end

end
