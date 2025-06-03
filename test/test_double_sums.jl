using Test
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics

const sqa = SecondQuantizedAlgebra

@testset "double_sums" begin
    N = 10
    ha = NLevelSpace(Symbol(:atom), 2)
    hf = FockSpace(:cavity)
    h = hf⊗ha

    ind(i) = Index(h, i, N, ha)

    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    g(k) = IndexedVariable(:g, k)

    innerSum = SingleSum(σ(2, 1, ind(:i))*σ(1, 2, ind(:j)), ind(:i))
    Dsum = DoubleSum(innerSum, ind(:j), [ind(:i)])
    @test(
        isequal(
            DoubleSum(innerSum, ind(:j)),
            DoubleSum(
                SingleSum(σ(2, 1, ind(:i))*σ(1, 2, ind(:j)), ind(:i), [ind(:j)]), ind(:j)
            ) + SingleSum(σ(2, 2, ind(:j)), ind(:j)),
        )
    )

    @test(
        isequal(
            DoubleSum(innerSum, ind(:j)),
            DoubleSum(SingleSum(σ(2, 1, ind(:i)), ind(:i))*σ(1, 2, ind(:j)), ind(:j)),
        )
    )

    N_atoms = 4
    N_modes = 2

    hf = FockSpace(:cavity)
    ha = NLevelSpace(Symbol(:atom), 2)
    h = hf ⊗ ha

    i_ind = Index(h, :i, N_atoms, ha)
    j_ind = Index(h, :j, N_atoms, ha)
    k_ind = Index(h, :k, N_modes, hf)
    l_ind = Index(h, :l, N_modes, hf)

    g_ik = IndexedVariable(:g, i_ind, k_ind)

    a(k) = IndexedOperator(Destroy(h, :a), k)
    σ(i, j, k) = IndexedOperator(Transition(h, Symbol("σ"), i, j), k)

    Ssum1 = Σ(g_ik*a(k_ind)*σ(2, 1, i_ind), i_ind)
    Ssum2 = Σ(g_ik*a(k_ind)'*σ(1, 2, i_ind), i_ind)

    ### issue 221 (DoubleSum)
    @cnumbers c1 N1
    i_ind2 = Index(h, :i, N1, ha)
    j_ind2 = Index(h, :j, N1, ha)
    @test isequal(
        simplify(Σ(-σ(2, 2, i_ind), i_ind, j_ind)), simplify(Σ(-σ(2, 2, i_ind), i_ind)*4)
    )
    @test isequal(
        simplify(Σ(3*σ(2, 2, i_ind), i_ind, j_ind)), simplify(Σ(3*σ(2, 2, i_ind), i_ind)*4)
    )
    @test isequal(
        simplify(Σ(c1*σ(2, 2, i_ind), i_ind, j_ind)),
        simplify(Σ(c1*σ(2, 2, i_ind), i_ind)*4),
    )
    @test isequal(
        simplify(Σ(-σ(2, 2, i_ind2), i_ind2, j_ind2)),
        simplify(Σ((1-N1)*σ(2, 2, i_ind2), i_ind2)) - Σ(σ(2, 2, i_ind2), i_ind2),
    )
    @test isequal(
        simplify(Σ(3*σ(2, 2, i_ind2), i_ind2, j_ind2)),
        simplify(Σ(3*(N1-1)*σ(2, 2, i_ind2), i_ind2)) + 3*Σ(σ(2, 2, i_ind2), i_ind2),
    )
    @test isequal(
        simplify(Σ(c1*σ(2, 2, i_ind2), i_ind2, j_ind2)),
        simplify(c1*Σ(σ(2, 2, i_ind2), i_ind2) + Σ(c1*(N1-1)*σ(2, 2, i_ind2), i_ind2)),
    )
    ###

    @test isequal(Σ(conj(g_ik)*a(k_ind)'*σ(1, 2, i_ind), i_ind), Ssum1')

    @test isequal(
        Σ(Σ(g_ik*(a(k_ind)*σ(2, 1, i_ind) + a(k_ind)'*σ(1, 2, i_ind)), i_ind), k_ind),
        Σ(Ssum1+Ssum2, k_ind),
    )

    @test Ssum1 isa SingleSum
    @test Σ(Ssum1, k_ind) isa DoubleSum

    H = Σ(Σ(g_ik*(a(k_ind)*σ(2, 1, i_ind) + a(k_ind)'*σ(1, 2, i_ind)), i_ind), k_ind)

    for arg in H.arguments
        @test arg isa DoubleSum
    end

    H2 = H*a(l_ind)

    @test Ssum1*a(l_ind) isa SingleSum

    DSum1 = H.arguments[1]
    Dsum2 = H.arguments[2]

    @test isequal(DSum1*a(l_ind), Σ(Ssum1*a(l_ind), k_ind))

    @test isequal(
        Σ(Σ(σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind, [j_ind]), j_ind, [i_ind]),
        Σ(Σ(σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind, [j_ind]), j_ind),
    )
    @test isequal(
        Σ(Σ(σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind, [j_ind]), j_ind, [i_ind]),
        Σ(Σ(σ(2, 1, i_ind), i_ind, [j_ind])*σ(1, 2, j_ind), j_ind),
    )

    innerSum = Σ(σ(1, 2, i_ind)*σ(2, 1, j_ind), i_ind)
    DSum = Σ(innerSum, j_ind)
    @test innerSum isa sqa.QAdd
    @test Dsum isa sqa.QAdd

    @test DSum isa sqa.QAdd

    @test Dsum.arguments[1] isa sqa.DoubleSum
    @test Dsum.arguments[2] isa sqa.SingleSum

    @test isequal(
        Σ(Σ(σ(1, 2, i_ind)*σ(2, 1, j_ind), i_ind, [j_ind]), j_ind) +
        Σ(σ(1, 2, j_ind)*σ(2, 1, j_ind), j_ind),
        DSum,
    )

    sum1 = Σ(σ(1, 2, i_ind), i_ind)
    sum2 = Σ(σ(2, 1, i_ind), i_ind)

    mul1 = *(sum1, sum2; ind=j_ind)

    sum2_ = Σ(σ(2, 1, j_ind), j_ind)
    mul2 = sum1*sum2_

    @test mul2 isa sqa.QAdd
    @test mul1 isa sqa.QAdd

    @test isequal(mul1, mul2)

    # Double indexed variable
    @cnumbers N
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    Γ(i, j) = IndexedVariable(:Γ, i, j)
    Ω(i, j) = IndexedVariable(:Ω, i, j; identical=false)
    @test iszero(Ω(i, i))
    @test iszero(Ω(2, 2))
    @test !isequal(Ω(2, 3), 0)
    @test !isequal(Γ(i, i), 0)
    @test !isequal(Γ(2, 2), 0)

    ### issue 223
    ha = NLevelSpace(:atom, 2)
    σ(α, β, i) = IndexedOperator(Transition(ha, :σ, α, β), i)
    @cnumbers N g
    i = Index(ha, :i, N, 1)
    j = Index(ha, :j, N, 1)
    k = Index(ha, :k, N, 1)
    #
    H = Σ(2*σ(2, 2, j), i, j)
    H_ji = Σ(2*σ(2, 2, j), j, i)
    H_s = simplify(H)
    H_g = Σ(g*σ(2, 2, j), i, j)
    H_ji_g = Σ(g*σ(2, 2, j), j, i)
    H_ji_g_s = simplify(Σ(g*σ(2, 2, j), j, i))

    dict_N = Dict(N => 10)
    sub_dict(x) = simplify(substitute(x, dict_N))
    #
    @test isequal(sub_dict(simplify(commutator(H, σ(2, 1, k)))), 20*σ(2, 1, k))
    @test isequal(sub_dict(simplify(commutator(H_s, σ(2, 1, k)))), 20*σ(2, 1, k))
    #
    @test isequal(sub_dict(simplify(commutator(H, σ(1, 2, k)))), -20*σ(1, 2, k))
    @test isequal(sub_dict(simplify(commutator(H_s, σ(1, 2, k)))), -20*σ(1, 2, k))

    @test isequal(sub_dict(simplify(commutator(H_g, σ(2, 1, k)))), 10*g*σ(2, 1, k))
    @test isequal(sub_dict(simplify(commutator(H_ji_g, σ(2, 1, k)))), 10*g*σ(2, 1, k))
    @test isequal(sub_dict(simplify(commutator(H_ji_g_s, σ(2, 1, k)))), 10*g*σ(2, 1, k))
    #
    @test isequal(sub_dict(simplify(commutator(H_g, σ(1, 2, k)))), -10g*σ(1, 2, k))
    @test isequal(sub_dict(simplify(commutator(H_ji_g, σ(1, 2, k)))), -10g*σ(1, 2, k))
    @test isequal(sub_dict(simplify(commutator(H_ji_g_s, σ(1, 2, k)))), -10g*σ(1, 2, k))
end
