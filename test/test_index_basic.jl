using Test
using SecondQuantizedAlgebra
using QuantumOpticsBase
using SymbolicUtils
using Symbolics

const sqa = SecondQuantizedAlgebra

@testset "index_basic" begin
    @testset "Index Creation and Comparison" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha) # transition index
        indF(i) = Index(h, i, N, hf) # fock index
        i_ind = indT(:i)
        j_ind = indT(:j)

        @test !isequal(indT(:i), indT(:j))
        @test !isequal(indT(:i), indF(:j))
        @test !isequal(indT(:i), indF(:i))
        @test isequal(indT(:i), Index(h, :i, 10, ha))
    end

    @testset "Indexed Variables" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)

        g(k) = IndexedVariable(:g, k)
        @test !isequal(g(indT(:i)), g(indT(:j)))
        @test isequal(g(indT(:i)), g(Index(h, :i, 10, ha)))

        Γij = DoubleIndexedVariable(:Γ, i_ind, j_ind)
        k_ind = indT(:k)
        @test isequal(
            change_index(Γij, j_ind, k_ind), DoubleIndexedVariable(:Γ, i_ind, k_ind)
        )
        @test isequal(change_index(g(k_ind), k_ind, j_ind), g(j_ind))
    end

    @testset "Indexed Operators and Basic Operations" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
        σ12i = σ(1, 2, indT(:i))

        @test isequal(σ12i, σ(1, 2, i_ind))
        @test !isequal(σ12i, σ(2, 2, i_ind))
        @test !isequal(σ12i, σ(1, 2, j_ind))

        @test isequal(0, σ12i*σ(1, 2, i_ind))
        @test isequal(σ(2, 2, i_ind), σ(2, 1, i_ind)*σ12i)
        @test isequal(adjoint(σ(1, 2, i_ind)), σ(2, 1, i_ind))
        @test isequal(acts_on(σ12i), 2)
        @test i_ind < j_ind
    end

    @testset "Single Sums" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        ind(a) = indT(a)
        i_ind = indT(:i)
        j_ind = indT(:j)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
        g(k) = IndexedVariable(:g, k)
        a = Destroy(h, :a)

        sum1 = SingleSum(σ(1, 2, i_ind)*a', i_ind)
        sum2 = SingleSum(σ(2, 1, i_ind)*a, i_ind)
        sum3 = SingleSum(a'*σ(1, 2, i_ind) + a*σ(2, 1, i_ind), i_ind)

        @test isequal(adjoint(sum1), sum2)
        @test isequal(sum3, (sum1+sum2))
        @test isequal(0, Σ(0, i_ind))
        @test isequal(0, Σ(σ(2, 1, i_ind)*σ(2, 1, i_ind), i_ind))

        # Sum equivalences
        Γij = DoubleIndexedVariable(:Γ, i_ind, j_ind)
        @test isequal(N*g(ind(:j)), Σ(g(ind(:j)), ind(:i)))
        @test Σ(g(ind(:j)), ind(:j)) isa sqa.SingleSum
        @test isequal(N*Γij, Σ(Γij, ind(:k)))
        @test Σ(Γij, ind(:i)) isa sqa.SingleSum

        # Sum notation equivalence
        @test isequal(∑(σ(1, 2, i_ind), i_ind), Σ(σ(1, 2, i_ind), i_ind))
        @test isequal(change_index(∑(2g(i_ind), i_ind), i_ind, j_ind), ∑(2g(j_ind), j_ind))
    end

    @testset "Double Sums" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

        innerSum = SingleSum(σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind)
        @test isequal(
            DoubleSum(innerSum, j_ind),
            DoubleSum(SingleSum(σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind, [j_ind]), j_ind) +
            SingleSum(σ(2, 2, j_ind), j_ind),
        )

        @test isequal(
            ∑(σ(1, 2, i_ind)*σ(2, 1, j_ind), i_ind, j_ind),
            Σ(σ(1, 2, i_ind)*σ(2, 1, j_ind), i_ind, j_ind),
        )
    end

    @testset "Index Reordering and Special Terms" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)
        k_ind = indT(:k)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

        @test isequal(change_index(σ(1, 2, j_ind)*σ(1, 2, i_ind), j_ind, i_ind), 0)
        @test isequal(
            order_by_index(σ(1, 2, k_ind)*σ(1, 2, j_ind)*σ(1, 2, i_ind), [i_ind]),
            σ(1, 2, i_ind)*σ(1, 2, k_ind)*σ(1, 2, j_ind),
        )

        @test isequal(
            reorder(σ(1, 2, k_ind)*σ(1, 2, j_ind)*σ(1, 2, i_ind), [(i_ind, j_ind)]),
            SpecialIndexedTerm(
                σ(1, 2, k_ind)*σ(1, 2, i_ind)*σ(1, 2, j_ind), [(i_ind, j_ind)]
            ),
        )

        # Test reordering edge cases
        Ωij = DoubleIndexedVariable(:Ω, i_ind, j_ind; identical=false)
        @test change_index(Ωij, i_ind, j_ind) == 0
        @test reorder(sqa.QAdd([]), [(i_ind, j_ind)]) == 0
        @test reorder(sqa.QAdd([0]), [(i_ind, j_ind)]) == 0
        @test reorder(sqa.QAdd([σ(1, 2, i_ind), σ(2, 1, j_ind)]), [(i_ind, j_ind)]) isa
            sqa.QAdd
        @test reorder(average(sqa.QAdd([0])), [(i_ind, j_ind)]) == 0
    end

    @testset "Addition Operations (QAdd)" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
        g(k) = IndexedVariable(:g, k)
        a = Destroy(h, :a)

        sum1 = SingleSum(σ(1, 2, i_ind)*a', i_ind)
        qadd = a + a'
        qmul = a'*a

        # Basic addition tests
        @test (sum1 + a') isa sqa.QAdd
        @test (sum1 + σ(1, 2, i_ind)) isa sqa.QAdd
        @test (σ(2, 1, j_ind) + σ(1, 2, i_ind)) isa sqa.QAdd
        @test (sum1 + g(i_ind)) isa sqa.QAdd
        @test isequal(sum1 + g(i_ind), g(i_ind) + sum1)
        @test (a + σ(1, 2, i_ind)) isa sqa.QAdd
        @test (σ(1, 2, i_ind)+a) isa sqa.QAdd

        # Complex addition tests
        @test length((qadd + sum1).arguments) == 3
        @test isequal((sum1+qadd), (qadd + sum1))
        @test length((qadd + σ(1, 2, i_ind)).arguments) == 3
        @test isequal((σ(1, 2, i_ind)+qadd), (qadd + σ(1, 2, i_ind)))
        @test length((qadd + g(i_ind)).arguments) == 3
        @test isequal((g(i_ind)+qadd), (qadd + g(i_ind)))

        @test sum1+qmul isa sqa.QAdd
        @test isequal((g(i_ind)+qmul), (qmul + g(i_ind)))
        @test isequal(g(i_ind) + σ(1, 2, j_ind), σ(1, 2, j_ind) + g(i_ind))

        # Special indexed terms
        specTerm = sqa.SpecialIndexedTerm(σ(1, 2, i_ind)*σ(1, 2, j_ind), [(i_ind, j_ind)])
        @test isequal((g(i_ind)+specTerm), (specTerm + g(i_ind)))
        @test isequal((specTerm+qadd), (qadd + specTerm))
        @test isequal((specTerm+2), (2 + specTerm))

        # Negation
        @test isequal(-σ(1, 2, i_ind), -1*σ(1, 2, i_ind))
        @test isequal(-g(i_ind), -1*g(i_ind))
        @test isequal(g(i_ind) + a, a + g(i_ind))
        @test isequal(qadd+g(i_ind), g(i_ind)+qadd)
    end

    @testset "Multiplication Operations (QMul)" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)
        k_ind = indT(:k)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
        g(k) = IndexedVariable(:g, k)
        a = Destroy(h, :a)
        sum1 = SingleSum(σ(1, 2, i_ind)*a', i_ind)
        qadd = a + a'
        qmul = a'*a

        # Basic multiplication tests
        @test g(i_ind)*a isa sqa.QMul
        @test (g(i_ind)*a).args_nc == [a]
        @test g(i_ind)*a' isa sqa.QMul
        @test (g(i_ind)*a').args_nc == [a']

        @test a*g(i_ind) isa sqa.QMul
        @test (a*g(i_ind)).args_nc == [a]
        @test a'*g(i_ind) isa sqa.QMul
        @test (a'*g(i_ind)).args_nc == [a']

        @test σ(1, 2, i_ind)*g(i_ind) isa sqa.QMul
        @test (σ(1, 2, i_ind)*g(i_ind)).args_nc == [σ(1, 2, i_ind)]
        @test g(i_ind)*σ(1, 2, i_ind) isa sqa.QMul
        @test (g(i_ind)*σ(1, 2, i_ind)).args_nc == [σ(1, 2, i_ind)]

        @test length((qmul*g(i_ind)).args_nc) == 2
        @test isequal((qmul*g(i_ind)).arg_c, g(i_ind))
        @test length((g(i_ind)*qmul).args_nc) == 2
        @test isequal((g(i_ind)*qmul).arg_c, g(i_ind))

        # Sum multiplication
        @test isequal(
            σ(1, 2, k_ind) * sum1,
            simplify(SingleSum(σ(1, 2, k_ind)*σ(1, 2, i_ind)*a', i_ind)),
        )

        @test isequal(
            simplify(σ(2, 1, k_ind) * sum1),
            simplify(
                SingleSum(σ(2, 1, k_ind)*σ(1, 2, i_ind)*a', i_ind, [k_ind]) +
                a'*σ(2, 2, k_ind),
            ),
        )

        # Special indexed term multiplication
        specTerm = sqa.SpecialIndexedTerm(σ(1, 2, i_ind)*σ(1, 2, j_ind), [(i_ind, j_ind)])
        asdf = specTerm*σ(1, 2, k_ind)
        asdf2 = σ(1, 2, k_ind)*specTerm
        @test isequal(asdf, asdf2)
        @test isequal(specTerm*qmul, qmul*specTerm)
        @test isequal(qadd*specTerm, specTerm*qadd)
        @test isequal(2*specTerm, specTerm*2)
    end

    @testset "Commutators" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
        a = Destroy(h, :a)
        qadd = a + a'
        qmul = a'*a

        @test isequal(
            commutator(σ(1, 2, i_ind), σ(2, 1, i_ind)),
            σ(1, 2, i_ind)*σ(2, 1, i_ind) - σ(2, 1, i_ind)*σ(1, 2, i_ind),
        )
        @test isequal(simplify(commutator(σ(1, 2, i_ind), qadd)), 0)
        @test isequal(simplify(commutator(σ(1, 2, i_ind), qmul)), 0)
    end

    @testset "Fock Operators and Commutation" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indF(i) = Index(h, i, N, hf)
        ai(k) = IndexedOperator(Destroy(h, :a), k)

        @test isequal((ai(indF(:m))*ai(indF(:m))'), ai(indF(:m))'*ai(indF(:m)) + 1)
    end

    @testset "Utility Functions and Equivalences" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        indT(i) = Index(h, i, N, ha)
        i_ind = indT(:i)
        j_ind = indT(:j)

        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

        @test isequal(NumberedOperator(Transition(h, :σ, 1, 2), 1), σ(1, 2, 1))
        @test isequal([i_ind, j_ind], sqa.get_indices(σ(1, 2, i_ind) + σ(2, 1, j_ind)))
        @test isequal(
            [i_ind, j_ind],
            sort(sqa.get_indices(average(σ(1, 2, i_ind)) + 3 + average(σ(2, 1, j_ind)))),
        )

        @test isequal(IndexedVariable(:Ω, 1, 2), sqa.DoubleNumberedVariable(:Ω, 1, 2))
        @test isequal(IndexedVariable(:Ω, 2), sqa.SingleNumberedVariable(:Ω, 2))
    end

    @testset "Multi-Space Operations" begin
        N = 10
        hc = FockSpace(:cavity)
        hf = FockSpace(:filter)
        h = hc ⊗ hf

        i = Index(h, :i, N, hf)
        j = Index(h, :j, N, hf)
        k = Index(h, :k, N, hf)

        b(k) = IndexedOperator(Destroy(h, :b, 2), k)

        @test reorder(b(i)*b(k)*b(i)', [(i, k)]) isa sqa.QAdd
        @test reorder(b(i)'*b(i)*b(k), [(i, k)]) isa sqa.SpecialIndexedTerm
        @test isequal(
            reorder(b(i)*b(k)*b(i)', [(i, k)]),
            reorder(b(i)'*b(i)*b(k), [(i, k)]) + reorder(b(k), [(i, k)]),
        )
    end

    @testset "Basis Conversion" begin
        hfock = FockSpace(:fock)
        bfock = FockBasis(3)
        hnlevel = NLevelSpace(:nlevel, 2)
        bnlevel = NLevelBasis(2)
        h_ = hnlevel ⊗ hfock

        N_n = 4
        N_f = 2
        ranges = [N_n, N_f]

        b_1 = ⊗([bnlevel for i in 1:N_n]...)
        b_2 = ⊗([bfock for i in 1:N_f]...)
        b_ = b_1 ⊗ b_2

        ai(i) = IndexedOperator(Destroy(h_, :a), i)
        σi(i, j, k) = IndexedOperator(Transition(h_, :σ, i, j), k)

        @test to_numeric(ai(1), b_; ranges=ranges) == LazyTensor(b_, [5], (destroy(bfock),))
        @test to_numeric(ai(2), b_; ranges=ranges) == LazyTensor(b_, [6], (destroy(bfock),))
        @test to_numeric(σi(1, 2, 4), b_; ranges=ranges) ==
            LazyTensor(b_, [4], (QuantumOpticsBase.transition(bnlevel, 1, 2),))
        @test_throws MethodError to_numeric(σi(1, 2, 5), b_; ranges=ranges)

        ai2(i) = IndexedOperator(Destroy(hfock, :a), i)
        @test to_numeric(ai2(1), b_2; ranges=[2]) isa LazyTensor
        @test to_numeric(ai2(2), b_2; ranges=[2]) isa LazyTensor
        @test_throws BoundsError to_numeric(ai2(3), b_2; ranges=[2])
    end

    @testset "Single Hilbert Space Operations" begin
        N = 10
        h = NLevelSpace(:atom, 2)
        i = Index(h, :i, N, h)
        σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)

        @test isa(i.hilb, ProductSpace) == false
        @test σ(1, 2, i) == (σ(1, 2, i)')'
    end

    @testset "Mixed Operator Types" begin
        hc = NLevelSpace(:cavity, 3)
        ha = NLevelSpace(:atom, 2)
        h = hc ⊗ ha

        @cnumbers N α
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, hc)

        S(x, y) = Transition(h, :S, x, y, 1)
        σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)

        @test S(2, 1)*σ(1, 2, i) isa SecondQuantizedAlgebra.QMul
        @test σ(1, 2, i)*S(2, 1) isa SecondQuantizedAlgebra.QMul
        @test σ(1, 2, 2)*S(2, 1) isa SecondQuantizedAlgebra.QMul
        @test S(2, 1)*σ(1, 2, 3) isa SecondQuantizedAlgebra.QMul
        @test σ(1, 2, 2) isa NumberedOperator
        @test isequal(S(2, 1)*σ(1, 2, i), σ(1, 2, i)*S(2, 1))
    end

    @testset "Index Array Creation and Utilities" begin
        hc = NLevelSpace(:cavity, 3)
        ha = NLevelSpace(:atom, 2)
        h = hc ⊗ ha

        @cnumbers N
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, hc)

        ranges1 = [1:10, 1:5]
        arr = sqa.create_index_arrays([i, j], ranges1)
        @test isequal(vec(collect(collect(Iterators.product(ranges1...)))), arr)

        arr = sqa.create_index_arrays([i], [1:10])
        @test isequal(1:10, arr)
    end

    @testset "Internal Utility Functions" begin
        ha = NLevelSpace(:atom, 2)
        σ(i, j, k) = IndexedOperator(Transition(ha, :σ, i, j), k)

        @test isequal(sqa.inorder!(σ(2, 1, 1)*σ(2, 2, 2)*σ(1, 2, 1)), σ(2, 2, 1)*σ(2, 2, 2))
        @test isequal(
            sqa.inadjoint(σ(2, 1, 1)*σ(2, 2, 2)*σ(1, 2, 1)), σ(2, 2, 1)*σ(2, 2, 2)
        )
        @test isequal(
            sqa._inconj(average(σ(2, 1, 1)*σ(2, 2, 2)*σ(1, 2, 1))),
            (average(σ(2, 2, 1)*σ(2, 2, 2))),
        )
        @test sqa.ismergeable(σ(2, 1, 5), σ(1, 2, 5))
    end

    @testset "QuantumCumulants Issue 188" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        i = Index(h, :i, N, ha)
        @cnumbers α
        gi = IndexedVariable(:g, i)

        @test isa(∑(5gi, i), SingleSum)
        @test isa(∑(gi*α, i), SingleSum)
        @test isequal(∑(α, i), N*α)
        @test isequal(∑(5α, i), 5*N*α)
    end
end
