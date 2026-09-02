using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: OP_DESTROY, OP_PAULI, OP_TRANSITION, rename, set_acts_on

@testset "Operator discovery and identity" begin
    @testset "fundamental operators describe each space" begin
        hf = FockSpace(:cavity)
        fock = fundamental_operators(hf)
        @test fock == [Destroy(hf, :a)]
        @test is_destroy(only(fock))
        @test acts_on(only(fock)) == [1]

        hn = NLevelSpace(:atom, 2, 1)
        @test Set(fundamental_operators(hn)) ==
            Set([Transition(hn, :σ, 1, 2), Transition(hn, :σ, 2, 2)])

        hn3 = NLevelSpace(:atom3, 3, 1)
        @test length(fundamental_operators(hn3)) == 5

        hp = PauliSpace(:pauli)
        @test all(is_pauli, fundamental_operators(hp))
        @test length(fundamental_operators(hp)) == 3

        hs = SpinSpace(:spin)
        @test all(is_spin, fundamental_operators(hs))
        @test length(fundamental_operators(hs)) == 3

        hq = PhaseSpace(:phase)
        phase_ops = fundamental_operators(hq)
        @test phase_ops == [Position(hq, :x), Momentum(hq, :p)]
    end

    @testset "product-space discovery and names" begin
        h = FockSpace(:f) ⊗ NLevelSpace(:atom, 2, 1)
        @test length(fundamental_operators(h)) == 3
        @test acts_on.(fundamental_operators(h)) == [[1], [2], [2]]

        named = fundamental_operators(h; names = [:b, :τ])
        @test operator_name.(named) == [:b, :τ, :τ]

        duplicate = FockSpace(:left) ⊗ FockSpace(:right)
        discovered = find_operators(duplicate, 1)
        @test Set(operator_name.(discovered)) == Set([:a, :b])
    end

    @testset "fundamental set is complete across families" begin
        h = FockSpace(:s1) ⊗ NLevelSpace(:s2, 2, 1) ⊗ FockSpace(:s3) ⊗
            PauliSpace(:s4) ⊗ SpinSpace(:s5) ⊗ PhaseSpace(:s6)
        expected = [
            Destroy(h, :a, 1),
            Transition(h, :σ, 1, 2, 2), Transition(h, :σ, 2, 2, 2),
            Destroy(h, :b, 3),
            Pauli(h, :σP, 1, 4), Pauli(h, :σP, 2, 4), Pauli(h, :σP, 3, 4),
            Spin(h, :S, 1, 5), Spin(h, :S, 2, 5), Spin(h, :S, 3, 5),
            Position(h, :x, 6), Momentum(h, :p, 6),
        ]
        @test fundamental_operators(h; names = [:a, :σ, :b, :σP, :S, (:x, :p)]) == expected
    end

    @testset "adjoint-aware uniqueness" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test unique_ops([a, a']) == [a]
        @test unique_ops([a', a]) == [a']

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        @test unique_ops([σx, σx']) == [σx]
        @test unique_ops([σx, σy]) == [σx, σy]
    end

    @testset "operator products can be discovered" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test find_operators(h, 1) == [a]
        @test find_operators(h, 2) == [a, a * a, a' * a]
    end
end

@testset "Operator accessors and adjoints" begin
    @testset "symbolic coefficients follow qadjoint" begin
        @variables g::Number r::Real
        h = FockSpace(:cavity)
        a = Destroy(h, :a)

        @test isequal(qadjoint(3 + 2im), 3 - 2im)
        @test isequal(qadjoint(g), conj(g))
        @test isequal(qadjoint(3 * g * a), 3 * conj(g) * a')
        @test isequal(qadjoint(qadjoint(g)), g)
        @test isequal(qadjoint(r), r)
        @test qconj === qadjoint
        @test dagger(a) == a'
        @test dagger(2im * a') == qadjoint(2im * a')
    end

    @testset "average adjoints stay averages" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test is_average(inner_adjoint(average(a)))
        @test iszero(undo_average(inner_adjoint(average(a))) - a')

        expression = 2 * average(a) + average(a' * a)
        @test isequal(
            inner_adjoint(expression),
            2 * average(a') + average(a' * a),
        )
    end

    @testset "operator and index identity accessors" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        a = Destroy(h, :a, 1)
        σ = Transition(h, :σ, 1, 2, 2)
        i = Index(h, :i, 3, 2)

        @test !has_index(operator_index(a))
        ai = IndexedOperator(a, i)
        @test operator_index(ai) == i
        @test acts_on(i) == [2]
        @test acts_on(i) == acts_on(σ)

        moved = set_acts_on(a, 2)
        @test acts_on(moved) == [2]
        @test operator_name(moved) == :a
        @test operator_index(moved) == operator_index(a)
        @test set_acts_on(a, 1) == a

        renamed = rename(a, :b)
        @test operator_name(renamed) == :b
        @test acts_on(renamed) == acts_on(a)
        @test optype(renamed) == optype(a)
        @test rename(a, :a) == a

        renamed_index = rename(i, :j)
        @test index_name(renamed_index) == :j
        @test index_range(renamed_index) == index_range(i)
        @test acts_on(renamed_index) == acts_on(i)
        @test_throws ArgumentError rename(a, "b")
        @test_throws ArgumentError rename(i, "j")
    end

    @testset "one, zero, and inference" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test isone(one(a))
        @test isone(commutator(a, a'))
        @test !isone(a)
        @test !isone(zero(a + a'))

        @test @inferred operator_index(a') isa typeof(operator_index(a))
        @test @inferred acts_on(a) isa Vector{Int}
        @test @inferred set_acts_on(a, 1) isa typeof(a)
        @test @inferred rename(a, :b) isa typeof(a)
    end

    @testset "optype contract is observable" begin
        @test optype(Destroy(FockSpace(:k), :a)) == OP_DESTROY
        @test optype(Transition(NLevelSpace(:n, 2), :σ, 1, 2)) == OP_TRANSITION
        @test optype(Pauli(PauliSpace(:p), :σ, 1)) == OP_PAULI
    end

    @testset "ProductSpace constructors and discovery edge cases" begin
        h = FockSpace(:f) ⊗ NLevelSpace(:n, 2, 1)
        @test acts_on(Destroy(FockSpace(:solo), :a)) == [1]
        @test length(find_operators(h, 1)) == length(fundamental_operators(h))
        dup = FockSpace(:left) ⊗ FockSpace(:right)
        @test Set(operator_name.(find_operators(dup, 1))) == Set([:a, :b])
        hs = SpinSpace(:s)
        hp = PauliSpace(:p)
        ops = [Destroy(FockSpace(:f), :a), Transition(NLevelSpace(:n, 2), :σ, 1, 2), Pauli(hp, :σ, 1), Spin(hs, :S, 1)]
        @test length(unique_ops(ops)) == length(ops)
        @test length(unique_ops([ops[1], ops[1]'])) == 1

        # Public auto-detect: single subspace of each type infers space_index
        h_mixed = FockSpace(:c) ⊗ NLevelSpace(:atom, 2) ⊗ PauliSpace(:p) ⊗
            SpinSpace(:s) ⊗ PhaseSpace(:q)
        @test Destroy(h_mixed, :a) == Destroy(h_mixed, :a, 1)
        @test Create(h_mixed, :a) == Create(h_mixed, :a, 1)
        @test Transition(h_mixed, :σ, 1, 2) == Transition(h_mixed, :σ, 1, 2, 2)
        @test Pauli(h_mixed, :p, 1) == Pauli(h_mixed, :p, 1, 3)
        @test Spin(h_mixed, :S, 2) == Spin(h_mixed, :S, 2, 4)
        @test Position(h_mixed, :x) == Position(h_mixed, :x, 5)
        @test Momentum(h_mixed, :p) == Momentum(h_mixed, :p, 5)
        h_sym = FockSpace(:c) ⊗ NLevelSpace(:atom, (:g, :e))
        @test Transition(h_sym, :σ, :g, :e) == Transition(h_sym, :σ, :g, :e, 2)
        h_jc = FockSpace(:c) ⊗ NLevelSpace(:atom, 2)
        @qnumbers a::Destroy(h_jc) σ::Transition(h_jc, 1, 2)
        @test a == Destroy(h_jc, :a, 1)
        @test σ == Transition(h_jc, :σ, 1, 2, 2)
        h_two_fock = FockSpace(:a) ⊗ FockSpace(:b)
        @test_throws ArgumentError Destroy(h_two_fock, :a)
        @test_throws ArgumentError Create(h_two_fock, :a)
        h_no_fock = NLevelSpace(:a, 2) ⊗ NLevelSpace(:b, 2)
        @test_throws ArgumentError Destroy(h_no_fock, :a)
        @test_throws ArgumentError Pauli(h_no_fock, :p, 1)
    end

    @testset "qadjoint reaches symbolic expressions" begin
        @variables g::Number
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        @test isequal(qadjoint(g * a + 2im * a'), conj(g) * a' - 2im * a)
        @test isequal(inner_adjoint(average(a)), average(a'))
    end
end
