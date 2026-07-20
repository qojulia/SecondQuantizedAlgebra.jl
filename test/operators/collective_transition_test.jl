using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: NO_INDEX, QAdd
using Latexify: latexify
using QuantumOpticsBase
using Test

@testset "collective transition" begin
    @testset "space and construction" begin
        h = CollectiveNLevelSpace(:atom, 3)
        S12 = CollectiveTransition(h, :S, 1, 2)
        @test is_collective_transition(S12)
        @test !is_transition(S12)
        @test S12.l1 == 1
        @test S12.l2 == 2
        @test S12.index == NO_INDEX
        @test operator_name(S12) == :S
        @test S12' == CollectiveTransition(h, :S, 2, 1)
        @test S12'' == S12
        @test S12 != CollectiveTransition(h, :S, 1, 3)
        @test S12 != CollectiveTransition(h, :T, 1, 2)
        @test hash(S12) == hash(CollectiveTransition(h, :S, 1, 2))

        @test CollectiveNLevelSpace(:atom, 3) == h
        @test CollectiveNLevelSpace(:atom, 2) != h
        @test_throws ArgumentError CollectiveNLevelSpace(:atom, 0)
        @test_throws ArgumentError CollectiveTransition(h, :S, 0, 1)
        @test_throws ArgumentError CollectiveTransition(h, :S, 1, 4)
        @test_throws MethodError CollectiveTransition(NLevelSpace(:atom, 3), :S, 1, 2)
        @test_throws MethodError Transition(h, :σ, 1, 2)
    end

    @testset "symbolic levels and product spaces" begin
        hc = CollectiveNLevelSpace(:atom, (:g, :e, :a))
        Sge = CollectiveTransition(hc, :S, :g, :e)
        @test Sge == CollectiveTransition(hc, :S, 1, 2)
        @test_throws ArgumentError CollectiveTransition(hc, :S, :missing, :e)

        h = FockSpace(:cavity) ⊗ hc
        @test CollectiveTransition(h, :S, 1, 2).space_index == 2
        @test CollectiveTransition(h, :S, :g, :e, 2) ==
            CollectiveTransition(h, :S, 1, 2, 2)
        @test_throws ArgumentError CollectiveTransition(h, :S, 1, 2, 1)

        duplicate = hc ⊗ CollectiveNLevelSpace(:other, 2)
        @test_throws ArgumentError CollectiveTransition(duplicate, :S, 1, 2)
        @test CollectiveTransition(duplicate, :T, 1, 2, 2).space_index == 2
    end

    @testset "fundamental operators and indexing" begin
        h = CollectiveNLevelSpace(:atom, 3)
        fund = fundamental_operators(h)
        @test length(fund) == 9
        @test all(is_collective_transition, fund)
        @test Set((Int(op.l1), Int(op.l2)) for op in fund) ==
            Set((i, j) for i in 1:3 for j in 1:3)

        S12 = CollectiveTransition(h, :S, 1, 2)
        i = Index(h, :i, 4, h)
        @test_throws ArgumentError IndexedOperator(S12, i)
        @test IndexedOperator(S12, NO_INDEX) == S12
    end

    @testset "su(N) algebra" begin
        h = CollectiveNLevelSpace(:atom, 3)
        S(i, j) = CollectiveTransition(h, :S, i, j)

        ordered = S(2, 1) * S(1, 3)
        @test ordered isa QAdd
        @test length(ordered) == 1
        @test only(keys(ordered.arguments)).ops == [S(2, 1), S(1, 3)]

        @test iszero(simplify(S(1, 2) * S(2, 3) - S(1, 3) - S(2, 3) * S(1, 2)))
        @test iszero(simplify(S(1, 2) * S(2, 1) - S(1, 1) + S(2, 2) - S(2, 1) * S(1, 2)))
        @test iszero(simplify(commutator(S(1, 2), S(2, 1)) - S(1, 1) + S(2, 2)))
        @test iszero(simplify(commutator(S(1, 2), S(2, 3)) - S(1, 3)))
        @test iszero(commutator(S(1, 2), S(1, 3)))
        @test expand_completeness(S(1, 1)) == simplify(S(1, 1))

        square = S(1, 2) * S(1, 2)
        @test length(square) == 1
        term = only(keys(square.arguments))
        @test term.ops == [S(1, 2), S(1, 2)]

        T12 = CollectiveTransition(h, :T, 1, 2)
        @test iszero(commutator(S(1, 2), T12))

        hp = h ⊗ CollectiveNLevelSpace(:other, 3)
        A12 = CollectiveTransition(hp, :A, 1, 2, 1)
        B21 = CollectiveTransition(hp, :B, 2, 1, 2)
        @test iszero(commutator(A12, B21))

        hm = FockSpace(:cavity) ⊗ h
        a = Destroy(hm, :a, 1)
        Sp = CollectiveTransition(hm, :S, 1, 2, 2)
        Sm = Sp'
        H = a * Sp + a' * Sm
        @test iszero(simplify(H' - H))
    end

    @testset "numeric conversion" begin
        N = 4
        h = CollectiveNLevelSpace(:atom, 3)
        S(i, j) = CollectiveTransition(h, :S, i, j)
        b1 = NLevelBasis(3)
        b = ManyBodyBasis(b1, bosonstates(b1, N))

        @test to_numeric(S(1, 2), b) == manybodyoperator(b, transition(b1, 1, 2))
        @test to_numeric(S(1, 2)', b) == QuantumOpticsBase.dagger(to_numeric(S(1, 2), b))
        @test sum(to_numeric(S(k, k), b) for k in 1:3) ≈ N * one(b)
        @test to_numeric(commutator(S(1, 2), S(2, 1)), b) ≈
            to_numeric(S(1, 1) - S(2, 2), b)

        wrong_onebody = FockBasis(3)
        wrong = ManyBodyBasis(wrong_onebody, bosonstates(wrong_onebody, 2))
        @test_throws ArgumentError to_numeric(S(1, 2), wrong)
        @test_throws ArgumentError to_numeric(Destroy(FockSpace(:f), :a), b)

        # A collective transition has no realisation on a plain (non-ManyBody) basis: the
        # closed operator ladder falls through to the throwing `Val` extension point.
        @test_throws ArgumentError to_numeric(S(1, 2), b1)

        # su(2) collective transitions reproduce spin-N/2 generators.
        h2 = CollectiveNLevelSpace(:twolevel, 2)
        C(i, j) = CollectiveTransition(h2, :S, i, j)
        b2 = NLevelBasis(2)
        bm = ManyBodyBasis(b2, bosonstates(b2, N))
        bs = SpinBasis(N // 2)
        @test to_numeric(C(1, 2), bm).data ≈ sigmap(bs).data
        @test to_numeric(C(2, 1), bm).data ≈ sigmam(bs).data
        @test to_numeric((C(1, 1) - C(2, 2)) / 2, bm).data ≈
            (0.5 * sigmaz(bs)).data
    end

    @testset "printing and type stability" begin
        h = CollectiveNLevelSpace(:atom, 3)
        S12 = CollectiveTransition(h, :S, 1, 2)
        S21 = S12'
        @test sprint(show, h) == "ℋ(atom)"
        @test sprint(show, S12) == "S₁₂"
        @test string(latexify(S12)) == raw"${S}^{{12}}$"
        SecondQuantizedAlgebra.transition_superscript(false)
        @test string(latexify(S12)) == raw"${S}_{{12}}$"
        SecondQuantizedAlgebra.transition_superscript(true)

        @inferred CollectiveTransition(h, :S, 1, 2)
        @inferred adjoint(S12)
        @inferred isequal(S12, S21)
        @inferred hash(S12, UInt(0))
        @inferred S12 * S21

        @static if VERSION >= v"1.12"
            @test @allocations(CollectiveTransition(h, :S, 1, 2)) == 0
            @test @allocations(adjoint(S12)) == 0
            @test @allocations(isequal(S12, S21)) == 0
            @test @allocations(hash(S12, UInt(0))) == 0
            S12 * S21 # warm the commute pipeline
            @test @allocations(S12 * S21) <= 20
        end
    end
end
