using SecondQuantizedAlgebra
import SecondQuantizedAlgebra as SQA
using SecondQuantizedAlgebra:
    OperatorKind, CustomOperator, operator_kind, operator_label, optype, canonical_isless,
    QAdd, to_qadd
using Test
using Symbolics: @variables
using Latexify: latexify
import QuantumOpticsBase as QOB
import QuantumToolbox

const SU11 = OperatorKind(
    :su11, FockSpace, 3;
    adjoint = (3, 2, 1),
    relations = function (K)
        Kp, K0, Km = K
        return ((K0, Kp) => Kp * K0 + Kp, (Km, K0) => K0 * Km + Km, (Km, Kp) => Kp * Km + 2 * K0)
    end,
    realization = let a = Destroy(FockSpace(:mode), :a)
        (a' * a' / 2, (a' * a + 1 // 2) / 2, a * a / 2)
    end,
    suffixes = ("₊", "₀", "₋"), latex = ("_{+}", "_{0}", "_{-}"),
)

const FREE = OperatorKind(:free, FockSpace, 2; adjoint = (1, 2))

const PROJECTOR = OperatorKind(
    :projector, NLevelSpace, 2;
    adjoint = (1, 2),
    relations = P -> ((P[1], P[1]) => P[1], (P[2], P[2]) => P[2], (P[1], P[2]) => 0, (P[2], P[1]) => 0),
    realization = (Transition(NLevelSpace(:atom, 2), :σ, 1, 1), Transition(NLevelSpace(:atom, 2), :σ, 2, 2)),
)

ladder(dim) = [j == k + 1 ? sqrt(k) : 0.0 for k in 1:dim, j in 1:dim]

@testset "operator kinds defined outside the package" begin
    h = FockSpace(:mode)
    Kp, K0, Km = (CustomOperator(SU11, h, :K, g) for g in 1:3)

    @testset "construction and accessors" begin
        @test Kp isa Op
        @test optype(Kp) === SQA.OP_CUSTOM
        @test operator_kind(Kp) === SU11
        @test operator_label(K0) == 2
        @test repr(SU11) == "su11"
        @test isequal(Kp', Km) && isequal(K0', K0)
        @test operator_name(Kp) === :K
        a = Destroy(h, :a)
        @test_throws ArgumentError operator_kind(a)
        @test_throws ArgumentError operator_label(a)
        @test_throws ArgumentError CustomOperator(SU11, h, "K", 1)
        @test_throws ArgumentError CustomOperator(SU11, h, :K, 4)
        @test_throws ArgumentError CustomOperator(SU11, NLevelSpace(:atom, 2), :K, 1)
        @test (@inferred CustomOperator(SU11, h, :K, 1)) isa Op
    end

    @testset "su(1,1) commutation relations" begin
        @test isequal(commutator(K0, Kp), to_qadd(Kp))
        @test isequal(commutator(K0, Km), -Km)
        @test isequal(commutator(Km, Kp), 2 * K0)
        @test isequal(Km * Kp, Kp * Km + 2 * K0)
        @test isequal((Km * K0) * Kp, Km * (K0 * Kp))
        for x in (Kp, K0, Km), y in (Kp, K0, Km), z in (Kp, K0, Km)
            jacobi = commutator(x, commutator(y, z)) + commutator(y, commutator(z, x)) +
                commutator(z, commutator(x, y))
            @test iszero(jacobi)
        end
        C = K0 * K0 - (Kp * Km + Km * Kp) / 2
        @test all(iszero(commutator(C, K)) for K in (Kp, K0, Km))
        @test isequal((Km * Kp)', Km * Kp)
        @test (@inferred Km * Kp) isa QAdd
    end

    @testset "independence from other operators" begin
        a = Destroy(h, :a)
        L = CustomOperator(SU11, h, :L, 3)
        @test isequal(L * Kp, Kp * L)
        @test isequal(a * Km, Km * a)
        F, G = CustomOperator(FREE, h, :F, 1), CustomOperator(FREE, h, :F, 2)
        @test !isequal(F * G, G * F)
        @test isequal(F * Kp, Kp * F)
    end

    @testset "local identities" begin
        ha = NLevelSpace(:atom, 2)
        P(g) = CustomOperator(PROJECTOR, ha, :P, g)
        @test isequal(P(1) * P(1), to_qadd(P(1)))
        @test iszero(P(1) * P(2))
        @test iszero(P(2) * P(1) * P(2))
        @test to_numeric(P(2), ha, 2; backend = QuantumOpticsBackend()).data ≈ [0 0; 0 1]
    end

    @testset "indexed operators" begin
        @variables N
        hp = FockSpace(:mode) ⊗ NLevelSpace(:atom, 2)
        i, j = Index(hp, :i, N, 1), Index(hp, :j, N, 1)
        K = CustomOperator(SU11, hp, :K, 1)
        @test isequal(K, CustomOperator(SU11, hp, :K, 1, 1))
        Ki, Kj = IndexedOperator(K, i), IndexedOperator(K, j)
        K0i = IndexedOperator(CustomOperator(SU11, hp, :K, 2), i)
        @test isequal(Ki' * Ki, Ki * Ki' + 2 * K0i)
        @test length(Ki' * Kj) == 1
        @test isequal(average(Km * Kp), average(Kp * Km) + 2 * average(K0))
    end

    @testset "canonical order" begin
        F = CustomOperator(FREE, h, :F, 1)
        a = Destroy(h, :a)
        @test canonical_isless(Kp, K0) && canonical_isless(K0, Km)
        @test canonical_isless(F, Kp) && !canonical_isless(Kp, F)
        @test canonical_isless(a, F)
        @test sort([Km, F, a, Kp]; lt = canonical_isless) == [a, F, Kp, Km]
    end

    @testset "printing" begin
        @test repr(Km * Kp) == "2 * K₀ + K₊ * K₋"
        @test latexify(Km * Kp) == raw"$2 {K}_{0} + {K}_{+}{K}_{-}$"
        @test repr(CustomOperator(FREE, h, :F, 2)) == "F₂"
        @test latexify(CustomOperator(FREE, h, :F, 2)) == raw"${F}_{2}$"
    end

    @testset "numeric conversion" begin
        n = 24
        block(m) = Matrix(m)[1:(n - 6), 1:(n - 6)]
        P = Km * Km * Kp * Kp
        qo(x) = to_numeric(x, h, n; backend = QuantumOpticsBackend())
        qt(x) = to_numeric(x, h, n; backend = QuantumToolboxBackend())
        @test block(qo(P).data) ≈ block((qo(Km) * qo(Km) * qo(Kp) * qo(Kp)).data)
        @test qo(P).data ≈ qt(P).data
        a = ladder(n + 1)
        @test Matrix(qo(Kp).data) ≈ a' * a' / 2
        C = K0 * K0 - (Kp * Km + Km * Kp) / 2
        @test block(qo(C).data) ≈ [k == l ? -3 / 16 : 0.0 for k in 1:(n - 6), l in 1:(n - 6)]
        for be in (QuantumOpticsBackend(), QuantumToolboxBackend())
            ctx = SQA.NumericContext(be, SQA.numeric_basis(be, h, n))
            sub = SQA.numeric_subbasis(be, ctx.basis, 1)
            @test Base.promote_op(SQA.numeric_leaf, Op, typeof(ctx)) ===
                Base.promote_op(SQA.numeric_operator, typeof(be), Op, typeof(sub))
        end
        @test_throws ArgumentError qo(CustomOperator(FREE, h, :F, 1))
    end

    @testset "invalid kinds" begin
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 2; adjoint = (2, 2))
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 2; adjoint = (1, 3))
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 0; adjoint = ())
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 1; adjoint = (1,), suffixes = ("a", "b"))
        @test_throws ArgumentError OperatorKind(
            :bad, FockSpace, 1; adjoint = (1,), relations = B -> ((B[1], B[1]) => Destroy(FockSpace(:mode), :a),),
        )
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 1; adjoint = (1,), relations = B -> ((B[1], B[1]) => "no",))
        @test_throws ArgumentError OperatorKind(
            :bad, FockSpace, 1; adjoint = (1,), relations = B -> ((B[1], B[1]) => B[1], (B[1], B[1]) => 0),
        )
        @variables θ
        @test_throws ArgumentError OperatorKind(
            :bad, FockSpace, 1; adjoint = (1,), realization = (θ * Destroy(FockSpace(:mode), :a),),
        )
        @test_throws ArgumentError OperatorKind(:bad, FockSpace, 1; adjoint = (1,), realization = (1,))
        hp = FockSpace(:mode) ⊗ FockSpace(:other)
        @test_throws ArgumentError CustomOperator(SU11, hp, :K, 1)
        @test_throws ArgumentError CustomOperator(SU11, hp, :K, 1, 3)
    end

    @testset "redefinition" begin
        OperatorKind(:redefined, FockSpace, 1; adjoint = (1,))
        R = OperatorKind(:redefined, FockSpace, 1; adjoint = (1,), relations = B -> ((B[1], B[1]) => 1,))
        r = CustomOperator(R, h, :r, 1)
        @test isequal(r * r, one(r))
    end
end
