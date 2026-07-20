using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QAdd, QSym, _single_qadd, _to_cnum, _to_complex, _fold_const,
    _to_numeric_static, NumericContext
using QuantumOpticsBase
using Symbolics: @variables, substitute
import SymbolicUtils
using Test

# `to_numeric` materialises via `op_type` (default `sparse`), so the return type depends only
# on `op_type`, never on the shape of the expression: a bare operator, a product, and a sum
# all return a concrete `Operator`. `op_type = identity` opts into the natural lazy form
# (`LazyTensor`/`LazyProduct`/`LazySum`). The time-dependent form returns a native
# `TimeDependentSum`, evaluated via `H(t)` and compared through `expect`.
_dat(x) = dense(x).data

# Tiny backend used as a conformance test for the documented third-party static interface.
struct MockNumericBackend <: SecondQuantizedAlgebra.NumericBackend end
struct MockEagerOperator
    data::Matrix{ComplexF64}
end
_mock_identity(n) = ComplexF64[i == j for i in 1:n, j in 1:n]
SecondQuantizedAlgebra.numeric_basis(::MockNumericBackend, ::FockSpace, cutoff) = Int(cutoff) + 1
SecondQuantizedAlgebra.numeric_num_subsystems(::MockNumericBackend, ::Int) = 1
SecondQuantizedAlgebra.numeric_subbasis(::MockNumericBackend, n::Int, slot::Int) =
    slot == 1 ? n : throw(ArgumentError("mock basis has no subsystem slot $slot"))
SecondQuantizedAlgebra.numeric_operator(
    ::MockNumericBackend, ::SecondQuantizedAlgebra.Op, n::Int,
) = _mock_identity(n)
SecondQuantizedAlgebra.numeric_embed(::MockNumericBackend, ::Int, ::Int, leaf) = leaf
SecondQuantizedAlgebra.numeric_identity(::MockNumericBackend, n::Int) = _mock_identity(n)
function SecondQuantizedAlgebra.numeric_assemble(::MockNumericBackend, n::Int, terms)
    result = zeros(ComplexF64, n, n)
    for (coefficient, factors) in terms
        term = isempty(factors) ? _mock_identity(n) : foldl(*, factors)
        result .+= coefficient .* term
    end
    return result
end
SecondQuantizedAlgebra.numeric_materialize(::MockNumericBackend, assembled, ::Nothing) =
    MockEagerOperator(assembled)
SecondQuantizedAlgebra.numeric_materialize(
    ::MockNumericBackend, assembled, ::typeof(identity),
) = assembled

@testset "numeric conversion" begin
    @testset "third-party static backend contract" begin
        h = FockSpace(:mock)
        @qnumbers a::Destroy(h)

        eager = to_numeric(2 * a, h, 3; backend = MockNumericBackend())
        lazy = to_numeric(2 * a, h, 3; backend = MockNumericBackend(), op_type = identity)
        @test eager isa MockEagerOperator
        @test eager.data == 2 .* _mock_identity(4)
        @test lazy == eager.data
    end

    @testset "Single space basic" begin
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

        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(a' * a, b; op_type = identity) isa LazySum
        @test _dat(to_numeric(a' * a, b)) == _dat(create(b) * destroy(b))
        @test _dat(to_numeric(2 * a, b)) == _dat(2 * destroy(b))
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        @test result isa Operator
        @test to_numeric(a + a', b; op_type = identity) isa LazySum
        @test _dat(result) == _dat(destroy(b) + create(b))
    end

    @testset "Scalar" begin
        b = FockBasis(7)
        @test to_numeric(3, b) == 3 * one(b)
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)

        a1_num = to_numeric(a1, b)
        @test a1_num isa Operator
        # The lazy form is the vector-backed LazySum for every shape (single op included).
        @test to_numeric(a1, b; op_type = identity) isa LazySum
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
        @test to_numeric(σ12, bc) isa Operator
        @test to_numeric(σ12, bc; op_type = identity) isa LazySum
    end

    @testset "Type stability" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        ψ = fockstate(b, 2)

        # The lazy assembly is inference-stable and one concrete type for any term count: the
        # property the 5-arg vector-backed `LazySum` constructor buys. `sparse` materialization
        # (the default) is a top-level convenience and NOT on the `@inferred` contract, since
        # `sparse` of the abstract-eltype `LazySum` widens.
        ctx = NumericContext(QuantumOpticsBackend(), b, Dict{QSym, Union{}}())
        @test @inferred(_to_numeric_static(a' * a, ctx)) isa AbstractOperator
        @test @inferred(_to_numeric_static(2 * a + 3 * a' + 5, ctx)) isa AbstractOperator
        @test typeof(_to_numeric_static(a + a', ctx)) === typeof(_to_numeric_static(a + a' + a' * a, ctx))

        # to_numeric: leaf + QAdd + dict-substitution paths materialise a concrete operator.
        @test to_numeric(a, b) isa AbstractOperator
        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(2 * a + 3 * a' + 5, b) isa Operator
        @test to_numeric(a, b, Dict{QSym, Union{}}()) isa AbstractOperator

        # numeric_average: every BasicSymbolic branch infers to ComplexF64
        @test @inferred(numeric_average(a' * a, ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a) + average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(2 * average(a) * average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a)^2, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3.5 + 1im, ψ)) isa ComplexF64

        # `_to_complex(::Any) -> ComplexF64` is the keystone; canary against a
        # 7th overload tripping Julia's union-split budget.
        @test Base.return_types(_to_complex, (Any,))[1] === ComplexF64
        @test length(methods(_to_complex)) == 2
    end

    @testset "op_type materialization" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        # The default materialises `sparse`; `op_type = identity` keeps the lazy assembly.
        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(2 * a + 3 * a', b) isa Operator
        @test to_numeric(a, b) isa Operator
        @test to_numeric(a' * a, b; op_type = identity) isa LazySum
        @test to_numeric(2 * a + 3 * a', b; op_type = identity) isa LazySum
        @test to_numeric(a, b) == sparse(to_numeric(a, b; op_type = identity))
        @test _dat(to_numeric(a' * a, b; op_type = dense)) == _dat(to_numeric(a' * a, b))
    end

    @testset "numeric_average: Average expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        avg_a = average(a)
        @test numeric_average(avg_a, ψ) ≈ α

        avg_sum = average(a) + average(a')
        @test numeric_average(avg_sum, ψ) ≈ α + conj(α)

        avg_scaled = 2 * average(a)
        @test numeric_average(avg_scaled, ψ) ≈ 2α
    end

    @testset "to_numeric: Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        custom_op = 2 * destroy(b)
        d = Dict(a => custom_op)

        # Dict substitution replaces operator
        @test to_numeric(a, b, d) == custom_op
        # Adjoint not in dict, falls back to normal
        @test to_numeric(a', b, d) == create(b)

        # QAdd with Dict
        result_mul = to_numeric(a' * a, b, d)
        @test _dat(result_mul) == _dat(create(b) * custom_op)
    end

    @testset "numeric_average: Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        d = Dict{QSym, Any}()  # empty dict, same as no dict
        @test numeric_average(a, ψ, d) ≈ α

        avg_a = average(a)
        @test numeric_average(avg_a, ψ, d) ≈ α

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
            ref1 = LazyTensor(
                bprod_gap, [1, 3],
                (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            ref2 = LazyTensor(
                bprod_gap, [1, 3],
                (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            @test _dat(to_numeric(op1, bprod_gap)) == _dat(ref1)
            @test _dat(to_numeric(op2, bprod_gap)) == _dat(ref2)
        end
    end

    @testset "Large Hilbert space" begin
        hfock = FockSpace(:fock)
        @qnumbers a::Destroy(hfock)
        bfock = FockBasis(100)

        ref = 2 * create(bfock) + 2 * destroy(bfock)
        got = to_numeric(2 * a + 2 * a', bfock)
        @test isequal(_dat(ref), _dat(got))
        @test iszero(_dat(ref) - _dat(got))
        @test isequal(_dat(to_numeric(2 * a, bfock)), _dat(2 * destroy(bfock)))
    end

    @testset "numeric_average: product state" begin
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

    @testset "numeric_average: comprehensive expressions" begin
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

    @testset "numeric_average: vector of states" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im, 0.5 + 0.0im)
        ψs = [coherentstate(b, α) for α in αs]

        @test numeric_average(a, ψs) ≈ [α for α in αs]
        @test numeric_average(a' * a, ψs) ≈ [abs2(α) for α in αs]
        @test numeric_average(average(a), ψs) ≈ [α for α in αs]

        @test expect(a, ψs) ≈ numeric_average(a, ψs)
        @test expect(a' * a, ψs[1]) ≈ numeric_average(a' * a, ψs[1])
        @test numeric_average(average(a), ψs[1]) ≈ αs[1]

        empty_ψs = typeof(ψs[1])[]
        @test_throws ArgumentError numeric_average(a, empty_ψs)
    end

    @testset "numeric_average: vector of states with Dict" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im)
        ψs = [coherentstate(b, α) for α in αs]
        d = Dict{QSym, Any}()

        @test numeric_average(a, ψs, d) ≈ numeric_average(a, ψs)
        @test numeric_average(average(a' * a), ψs, d) ≈ [abs2(α) for α in αs]
        @test expect(a, ψs, d) ≈ numeric_average(a, ψs, d)
    end

    @testset "numeric_average: Dict comprehensive" begin
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
        @test _dat(to_numeric(ad * n, b_test, Dict{QSym, Any}())) == _dat(to_numeric(ad * n, b_test))
        @test _dat(to_numeric(2 * ad * a, b_test, Dict{QSym, Any}())) == _dat(to_numeric(2 * ad * a, b_test))
        @test to_numeric(ad, b_all, dd) == s(2, 2, 1)
        @test _dat(to_numeric(2 * ad * a, b_all, dd)) == _dat(2 * s(2, 2, 1) * s(2, 1, 2))
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

    @testset "Allocations: to_numeric" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        # Warmup
        to_numeric(a, b)
        to_numeric(a', b)
        to_numeric(a' * a, b)

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
        @test _dat(to_numeric(a' * a, b)) == _dat(create(b) * destroy(b))

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)

        expr = a + a' * a
        @test numeric_average(expr, ψ) ≈ α + abs2(α)

        # to_numeric with scalar QAdd (empty operators), now a lazy identity
        @test dense(to_numeric(_single_qadd(_to_cnum(3), Op[]), b)) == 3 * one(b)
    end

    @testset "Number-symtype coefficient round-trip" begin
        @variables g::Number
        b = FockBasis(3)
        a = Destroy(FockSpace(:c), :a)

        for v in (2 + 3im, 1 + 1im, 0.5 - 0.25im)
            op = substitute((g * a)' * (g * a), Dict(g => v))
            @test _dat(to_numeric(op, b)) ≈ abs2(v) * _dat(to_numeric(a' * a, b))
        end

        @test _dat(to_numeric(substitute(g * a, Dict(g => 2 + 3im)), b)) ≈ (2 + 3im) * _dat(to_numeric(a, b))
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
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables x::Real ϕ::Real E::Number

        @test _dat(to_numeric(substitute(sqrt(x) * a, Dict(x => 2.0)), b)) ≈ sqrt(2.0) * _dat(A)
        @test _dat(to_numeric(substitute(exp(im * ϕ) * a, Dict(ϕ => 0.5)), b)) ≈ exp(0.5im) * _dat(A)

        @test _dat(to_numeric(sqrt(x) * a, b; parameter = Dict(x => 2.0))) ≈ sqrt(2.0) * _dat(A)
        @test _dat(to_numeric(exp(im * ϕ) * a, b; parameter = Dict(ϕ => 0.5))) ≈ exp(0.5im) * _dat(A)

        custom = 2 * A
        @test _dat(to_numeric(a', b; operators = Dict(a => custom))) == _dat(custom')
        @test _dat(to_numeric(a', b; operators = Dict(a => custom), adjoint_ops = false)) == _dat(Ad)
        @test [_dat(x) for x in to_numeric([a, a'], b; operators = Dict(a => custom))] == [_dat(custom), _dat(custom')]

        # Time-dependent: native TimeDependentSum, evaluated at a time, compared via expect.
        f = to_numeric(E * a + conj(E) * a', b; time_parameter = Dict(E => t -> 1 + im * t))
        @test f isa TimeDependentSum
        @test expect(f(0.5), ψ) ≈ expect((1 + 0.5im) * A + (1 - 0.5im) * Ad, ψ)
    end

    @testset "keyword to_numeric on composite basis" begin
        hf = FockSpace(:c) ⊗ FockSpace(:d)
        a = Destroy(hf, :a, 1)
        @variables Δ::Real
        b = FockBasis(2) ⊗ FockBasis(3)
        Ia = identityoperator(b)
        an = to_numeric(a, b)

        # Regression: a scalar/constant term must build the full-system identity, not the
        # identity of a single subsystem (the `numeric_identity(composite)` fix).
        H = to_numeric(Δ * a' * a + 5, b; parameter = Dict(Δ => 2.0))
        @test size(dense(H).data) == (length(b), length(b))
        @test dense(H).data ≈ dense(2.0 * an' * an + 5 * Ia).data

        custom = 2 * an
        H2 = to_numeric(a + 3, b; operators = Dict(a => custom))
        @test dense(H2).data ≈ dense(custom + 3 * Ia).data

        z = to_numeric(a - a, b)
        @test size(dense(z).data) == (length(b), length(b))
        @test dense(z).data ≈ zeros(length(b), length(b))

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
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables E::Number

        # A plain number value becomes a constant-in-time term.
        f0 = to_numeric(E * a, b; time_parameter = Dict(E => 3.0))
        @test expect(f0(0.0), ψ) ≈ expect(3.0 * A, ψ)
        @test expect(f0(10.0), ψ) ≈ expect(3.0 * A, ψ)

        # A constant-coefficient term still returns a native TD operator.
        f1 = to_numeric(2.0 * a, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test expect(f1(7.0), ψ) ≈ expect(2.0 * A, ψ)

        # `conj(v)` is accepted as a time_parameter key.
        f2 = to_numeric(conj(E) * a, b; time_parameter = Dict(conj(E) => t -> 2 + im * t))
        @test expect(f2(1.0), ψ) ≈ expect((2 + 1im) * A, ψ)

        @test to_numeric(E * a, b; time_parameter = Dict(E => t -> 1), op_type = identity) isa TimeDependentSum
        @test_throws ArgumentError to_numeric(E * a, b; time_parameter = Dict(E => t -> 1), op_type = dense)
    end

    @testset "public backend hooks and product dims validation" begin
        b = FockBasis(3)
        ψ = fockstate(b, 0)
        @test numeric_backend(ψ) isa QuantumOpticsBackend
        @test numeric_basis(ψ) == b
        @test SecondQuantizedAlgebra.numeric_num_subsystems(QuantumOpticsBackend(), b) == 1

        h = FockSpace(:a) ⊗ FockSpace(:b)
        a = Destroy(h, :a, 1)
        @test_throws ArgumentError to_numeric(a, h, [2]; backend = QuantumOpticsBackend())
        @test_throws ArgumentError to_numeric(a, h, [2, 3, 99]; backend = QuantumOpticsBackend())
    end

    @testset "keyword to_numeric, scalar argument" begin
        b = FockBasis(4)
        Ib = one(b)
        ψ = coherentstate(b, 0.3 - 0.1im)
        @variables x::Real E::Number

        @test to_numeric(3, b; parameter = Dict(x => 1.0)) == 3 * Ib
        @test to_numeric(2.0 + 0im, b; parameter = Dict{Any, Any}()) == (2.0 + 0im) * Ib

        @test dense(to_numeric(2 * x, b; parameter = Dict(x => 1.5))) ≈ 3.0 * Ib

        @test_throws ArgumentError to_numeric(x, b)

        # With `time_parameter`, a constant scalar yields a constant-in-time TD operator.
        fconst = to_numeric(2.0, b; time_parameter = Dict(E => t -> 1.0 + 0im))
        @test expect(fconst(0.0), ψ) ≈ expect(2.0 * Ib, ψ)
        @test expect(fconst(9.0), ψ) ≈ expect(2.0 * Ib, ψ)

        ft = to_numeric(E, b; time_parameter = Dict(E => t -> 1 + im * t))
        @test expect(ft(0.0), ψ) ≈ expect((1 + 0im) * Ib, ψ)
        @test expect(ft(2.0), ψ) ≈ expect((1 + 2im) * Ib, ψ)
    end

    @testset "keyword to_numeric, state argument" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @variables Δ::Real

        op_state = to_numeric(Δ * a, ψ; parameter = Dict(Δ => 2.0))
        @test _dat(op_state) == _dat(2.0 * destroy(b))
    end

    @testset "constant symbolic coefficient reduction" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(3)
        A = destroy(b)
        @variables z::Number

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
            @test _dat(to_numeric(op, b)) ≈ expected * _dat(A)
        end

        op_sin = substitute(sin(z) * a, Dict(z => 0.5))
        @test _dat(to_numeric(op_sin, b)) ≈ sin(0.5) * _dat(A)

        neg_unary = SymbolicUtils.term(-, 5.0 + 0im; type = Number)
        neg_binary = SymbolicUtils.term(-, 7.0 + 0im, 2.0 + 0im; type = Number)
        @test _fold_const(neg_unary) == -(5.0 + 0im)
        @test _fold_const(neg_binary) == (7.0 + 0im) - (2.0 + 0im)
    end

    @testset "keyword to_numeric: vector, complex params, errors" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(4)
        A = destroy(b)
        Ad = create(b)
        @variables x::Real z::Number E::Number

        # `op_type` selects the materialization; `dense`/`sparse` on the default result work too.
        Hd = to_numeric(2.0 * a' * a, b; op_type = dense)
        @test Hd isa QuantumOpticsBase.Operator
        @test Hd ≈ dense(2.0 * Ad * A)
        @test dense(to_numeric(2.0 * a' * a, b)) ≈ Hd

        # Vector form forwards keywords to each element.
        @test [_dat(x) for x in to_numeric([a, a'], b; parameter = Dict(x => 1.0))] == [_dat(A), _dat(Ad)]

        # A complex-valued parameter key is split into real/imaginary substitutions.
        @test _dat(to_numeric(z * a, b; parameter = Dict(z => 2 + 3im))) ≈ (2 + 3im) * _dat(A)

        # Error paths.
        @test_throws ArgumentError to_numeric(x * a, b)                         # symbolic, no value
        @test_throws ArgumentError to_numeric(a, b; operators = Dict(x => A))   # non-QSym key
        @test_throws ArgumentError to_numeric(
            x * E * a, b; time_parameter = Dict(E => t -> 1.0 + 0im),
        )                                                                       # untimed variable
    end

    @testset "HilbertSpace form (uniform entry)" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        ψ = coherentstate(b, 0.2 + 0.1im)

        # Backend defaults to the single loaded backend (QuantumOpticsBase here).
        @test _dat(to_numeric(a' * a, h, 7)) == _dat(to_numeric(a' * a, b))
        @test numeric_average(a' * a, ψ) ≈ expect(to_numeric(a' * a, h, 7), ψ)

        # Composite HilbertSpace with per-subspace dims.
        hp = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        ap = Destroy(hp, :a, 1)
        σp = Transition(hp, :σ, 1, 2, 2)
        bp = FockBasis(4) ⊗ NLevelBasis(3)
        @test _dat(to_numeric(ap' * σp, hp, (4, 3))) == _dat(to_numeric(ap' * σp, bp))
    end

    @testset "op_type is shape-independent" begin
        # The return type depends only on op_type, not on the expression shape.
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)
        exprs = (a1, a1' * a1, a1' * a1 + a2' * a2)

        for expr in exprs
            @test to_numeric(expr, b) isa Operator
            @test to_numeric(expr, b) == to_numeric(expr, b; op_type = sparse)
            @test to_numeric(expr, b; op_type = dense) isa Operator
        end
        @test to_numeric(exprs[1], b; op_type = identity) isa LazySum
        @test to_numeric(exprs[3], b; op_type = identity) isa LazySum
        for expr in exprs
            @test to_numeric(expr, b) == sparse(to_numeric(expr, b; op_type = identity))
            @test dense(to_numeric(expr, b)) ≈ dense(to_numeric(expr, b; op_type = dense))
        end
    end

    @testset "numeric_average — LazyKet state" begin
        # `expect` of a sparse operator has no method for a `LazyKet`; numeric_average
        # assembles the lazy form, so a `LazyKet` state works without materializing.
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
            @test numeric_average(op, ψ_lazy) ≈ expect(to_numeric(op, bprod), ψ_dense)
        end
        @test numeric_average(average(a' * a), ψ_lazy) ≈
            expect(to_numeric(a' * a, bprod), ψ_dense)
    end

    @testset "numeric_average — unsupported symbolic operation" begin
        # A symbolic average expression whose top operation is neither `+`, `*`, nor `^`
        # (here `sqrt`) is not reducible to an expectation value and must error.
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(5)
        ψ = fockstate(b, 2)
        @test_throws ArgumentError numeric_average(sqrt(average(a' * a)), ψ)
        @test_throws ArgumentError numeric_average(average(a' * a) / (1 + average(a)), ψ)
    end

    @testset "complex parameter key" begin
        # A `Complex` parameter key `κ` splits into real/imag substitutions, so both `κ`
        # and `conj(κ)` in the expression resolve from a single complex value.
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(6)
        @variables κ::Complex
        H = κ * a + conj(κ) * a'
        M = to_numeric(H, b; parameter = Dict(κ => 1.0 + 2.0im))
        Mref = to_numeric((1.0 + 2.0im) * a + (1.0 - 2.0im) * a', b)
        @test _dat(M) ≈ _dat(Mref)
    end

    @testset "time_parameter — unsupported key" begin
        # A single-variable `time_parameter` key that is neither a bare variable nor
        # `conj(v)` (here `2*E`) is rejected.
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(6)
        @variables E::Number
        @test_throws ArgumentError to_numeric(E * a, b; time_parameter = Dict(2 * E => t -> 1.0 + 0im))
    end

end
