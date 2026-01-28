using Test
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics

const sqa = SecondQuantizedAlgebra

@testset "double_sums" begin
    @testset "Basic DoubleSum Construction and Equivalence" begin
        N = 10
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf⊗ha

        ind(i) = Index(h, i, N, ha)
        σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

        innerSum = SingleSum(σ(2, 1, ind(:i))*σ(1, 2, ind(:j)), ind(:i))
        Dsum = DoubleSum(innerSum, ind(:j), [ind(:i)])

        @test isequal(
            DoubleSum(innerSum, ind(:j)),
            DoubleSum(
                SingleSum(σ(2, 1, ind(:i))*σ(1, 2, ind(:j)), ind(:i), [ind(:j)]), ind(:j)
            ) + SingleSum(σ(2, 2, ind(:j)), ind(:j)),
        )

        @test isequal(
            DoubleSum(innerSum, ind(:j)),
            DoubleSum(SingleSum(σ(2, 1, ind(:i)), ind(:i))*σ(1, 2, ind(:j)), ind(:j)),
        )
    end

    @testset "Multi-Mode Multi-Atom Systems" begin
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
        expected = Σ(Ssum1*a(l_ind), k_ind)
        @test isequal(DSum1*a(l_ind), expected) || isequal(Dsum2*a(l_ind), expected)
    end

    @testset "Nested Sum Equivalences" begin
        N_atoms = 4
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

        i_ind = Index(h, :i, N_atoms, ha)
        j_ind = Index(h, :j, N_atoms, ha)
        σ(i, j, k) = IndexedOperator(Transition(h, Symbol("σ"), i, j), k)

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
        @test DSum isa sqa.QAdd

        @test isequal(
            Σ(Σ(σ(1, 2, i_ind)*σ(2, 1, j_ind), i_ind, [j_ind]), j_ind) +
            Σ(σ(1, 2, j_ind)*σ(2, 1, j_ind), j_ind),
            DSum,
        )
    end

    @testset "Sum Multiplication" begin
        N_atoms = 4
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

        i_ind = Index(h, :i, N_atoms, ha)
        j_ind = Index(h, :j, N_atoms, ha)
        σ(i, j, k) = IndexedOperator(Transition(h, Symbol("σ"), i, j), k)

        sum1 = Σ(σ(1, 2, i_ind), i_ind)
        sum2 = Σ(σ(2, 1, i_ind), i_ind)
        sum2_ = Σ(σ(2, 1, j_ind), j_ind)

        mul1 = *(sum1, sum2; ind=j_ind)
        mul2 = sum1*sum2_

        @test mul2 isa sqa.QAdd
        @test mul1 isa sqa.QAdd
        @test isequal(mul1, mul2)
    end

    @testset "Indexed Variables and Symmetries" begin
        ha = NLevelSpace(Symbol(:atom), 2)
        hf = FockSpace(:cavity)
        h = hf ⊗ ha

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
    end

    @testset "QuantumCumulants Issues" begin
        @testset "Issue 221 - DoubleSum Simplification" begin
            N = 10
            ha = NLevelSpace(Symbol(:atom), 2)
            hf = FockSpace(:cavity)
            h = hf⊗ha

            i_ind = Index(h, :i, N, ha)
            j_ind = Index(h, :j, N, ha)
            σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

            @cnumbers c1 N1
            i_ind2 = Index(h, :i, N1, ha)
            j_ind2 = Index(h, :j, N1, ha)

            # Tests with fixed N - corrected expected values
            @test isequal(
                simplify(Σ(-σ(2, 2, i_ind), i_ind, j_ind)),
                simplify(Σ(-σ(2, 2, i_ind), i_ind)*N),
            )
            @test isequal(
                simplify(Σ(3*σ(2, 2, i_ind), i_ind, j_ind)),
                simplify(Σ(3*σ(2, 2, i_ind), i_ind)*N),
            )
            @test isequal(
                simplify(Σ(c1*σ(2, 2, i_ind), i_ind, j_ind)),
                simplify(Σ(c1*σ(2, 2, i_ind), i_ind)*N),
            )

            # Tests with symbolic N1
            @test isequal(
                simplify(Σ(-σ(2, 2, i_ind2), i_ind2, j_ind2)),
                simplify(Σ((1-N1)*σ(2, 2, i_ind2), i_ind2)) - Σ(σ(2, 2, i_ind2), i_ind2),
            )
            @test isequal(
                simplify(Σ(3*σ(2, 2, i_ind2), i_ind2, j_ind2)),
                simplify(Σ(3*(N1-1)*σ(2, 2, i_ind2), i_ind2)) +
                3*Σ(σ(2, 2, i_ind2), i_ind2),
            )
            @test isequal(
                simplify(Σ(c1*σ(2, 2, i_ind2), i_ind2, j_ind2)),
                simplify(
                    c1*Σ(σ(2, 2, i_ind2), i_ind2) + Σ(c1*(N1-1)*σ(2, 2, i_ind2), i_ind2)
                ),
            )
        end

        @testset "Issue 223 - Commutator Relations" begin
            ha = NLevelSpace(:atom, 2)
            σ(α, β, i) = IndexedOperator(Transition(ha, :σ, α, β), i)

            @cnumbers N g
            i = Index(ha, :i, N, 1)
            j = Index(ha, :j, N, 1)
            k = Index(ha, :k, N, 1)

            H = Σ(2*σ(2, 2, j), i, j)
            H_ji = Σ(2*σ(2, 2, j), j, i)
            H_s = simplify(H)
            H_g = Σ(g*σ(2, 2, j), i, j)
            H_ji_g = Σ(g*σ(2, 2, j), j, i)
            H_ji_g_s = simplify(Σ(g*σ(2, 2, j), j, i))

            dict_N = Dict(N => 10)
            sub_dict(x) = simplify(substitute(x, dict_N))

            # Commutator tests with numerical substitution
            @test isequal(sub_dict(simplify(commutator(H, σ(2, 1, k)))), 20*σ(2, 1, k))
            @test isequal(sub_dict(simplify(commutator(H_s, σ(2, 1, k)))), 20*σ(2, 1, k))
            @test isequal(sub_dict(simplify(commutator(H, σ(1, 2, k)))), -20*σ(1, 2, k))
            @test isequal(sub_dict(simplify(commutator(H_s, σ(1, 2, k)))), -20*σ(1, 2, k))

            @test isequal(sub_dict(simplify(commutator(H_g, σ(2, 1, k)))), 10*g*σ(2, 1, k))
            @test isequal(
                sub_dict(simplify(commutator(H_ji_g, σ(2, 1, k)))), 10*g*σ(2, 1, k)
            )
            @test isequal(
                sub_dict(simplify(commutator(H_ji_g_s, σ(2, 1, k)))), 10*g*σ(2, 1, k)
            )
            @test isequal(sub_dict(simplify(commutator(H_g, σ(1, 2, k)))), -10g*σ(1, 2, k))
            @test isequal(
                sub_dict(simplify(commutator(H_ji_g, σ(1, 2, k)))), -10g*σ(1, 2, k)
            )
            @test isequal(
                sub_dict(simplify(commutator(H_ji_g_s, σ(1, 2, k)))), -10g*σ(1, 2, k)
            )
        end
    end
    @testset "DoubleSum Multiplication Fix (Issue #256 QC)" begin
        # This is the core test case from the GitHub issue
        @cnumbers N M γ Ω Δ V
        hA = NLevelSpace(:atomA, (:g, :r))
        hB = NLevelSpace(:atomB, (:g, :r))
        h = hA ⊗ hB

        σA(α, β, i) = IndexedOperator(Transition(h, Symbol("σ_A"), α, β, 1), i)
        σB(α, β, i) = IndexedOperator(Transition(h, Symbol("σ_B"), α, β, 2), i)

        i = Index(h, Symbol("i"), N, 1)  # index for subspace 1 (hA)
        j = Index(h, Symbol("j"), M, 2)  # index for subspace 2 (hB)
        k = Index(h, Symbol("k"), M, 2)  # index for subspace 2 (hB)

        op = σA(:g, :r, i)  # operator in subspace 1
        @testset "Basic IndexedOperator * DoubleSum with Different Subspaces" begin
            # This should work: IndexedOperator * SingleSum
            result_single = op * Σ(σB(:r, :r, j), j)
            @test result_single isa SingleSum

            # This was the failing case before the fix: IndexedOperator * DoubleSum
            # where the operator is in a different subspace than the sum
            doublesum = Σ(σB(:r, :r, j) * σB(:r, :r, k), j, k)
            @test doublesum.arguments[1] isa DoubleSum

            # This should NOT throw an AssertionError anymore
            result_double = op * doublesum
            @test result_double.arguments[1] isa DoubleSum

            # The result should be equivalent to multiplying the operator with the inner sum
            # and wrapping it in a DoubleSum with the same outer sum index
            doublesum = doublesum.arguments[1]
            result_double = result_double.arguments[1]
            expected = DoubleSum(
                op * doublesum.innerSum, doublesum.sum_index, doublesum.NEI
            )
            @test isequal(result_double, expected)
        end

        @testset "Symmetric Case: DoubleSum * IndexedOperator" begin
            doublesum = Σ(σB(:r, :r, j) * σB(:r, :r, k), j, k)
            doublesum = doublesum.arguments[1]

            # Test DoubleSum * IndexedOperator
            result = doublesum * op
            @test result isa DoubleSum

            # The result should be equivalent to multiplying the inner sum with the operator
            # and wrapping it in a DoubleSum with the same outer sum index
            expected = DoubleSum(
                doublesum.innerSum * op, doublesum.sum_index, doublesum.NEI
            )
            @test isequal(result, expected)
        end

        @testset "Complex Multi-Subspace Scenario" begin
            # Test a more complex scenario similar to the full example in the issue
            @cnumbers N1 N2 Va Vb V
            i1 = Index(h, Symbol("i"), N1, 1)
            j1 = Index(h, Symbol("j"), N1, 1)
            k1 = Index(h, Symbol("k"), N1, 1)
            i2 = Index(h, Symbol("i"), N2, 2)
            j2 = Index(h, Symbol("j"), N2, 2)
            k2 = Index(h, Symbol("k"), N2, 2)

            # Create DoubleSum terms as in the original issue
            VA = Va / 2 * Σ(σA(:r, :r, i1) * Σ(σA(:r, :r, j1), j1, [i1]), i1)
            VB = Vb / 2 * Σ(σB(:r, :r, i2) * Σ(σB(:r, :r, j2), j2, [i2]), i2)
            HAB = V / 2 * Σ(σA(:r, :r, i1) * σB(:r, :r, j2), i1, j2)

            # Test multiplying operators from different subspaces with these terms
            op_A = σA(:g, :r, k1)
            op_B = σB(:g, :r, k2)

            # These operations should not throw errors
            result_A_VB = op_A * VB  # operator from A times DoubleSum from B
            result_B_VA = op_B * VA  # operator from B times DoubleSum from A

            @test !iszero(result_A_VB)
            @test !iszero(result_B_VA)
        end

        @testset "Regression Test for Original Commutator Problem" begin
            # This tests the scenario from the issue where meanfield equations couldn't be generated
            @cnumbers N1 N2
            i1 = Index(h, Symbol("i"), N1, 1)
            j1 = Index(h, Symbol("j"), N1, 1)
            k1 = Index(h, Symbol("k"), N1, 1)
            i2 = Index(h, Symbol("i"), N2, 2)
            j2 = Index(h, Symbol("j"), N2, 2)

            # Create the problematic terms from the issue
            VA = Σ(σA(:r, :r, i1) * Σ(σA(:r, :r, j1), j1, [i1]), i1)
            VB = Σ(σB(:r, :r, i2) * Σ(σB(:r, :r, j2), j2, [i2]), i2)

            # This should work now: operators from different subspaces
            op_A = σA(:g, :r, k1)
            op_B = σB(:g, :r, j2)

            # These operations were causing AssertionErrors before the fix
            comm_A_VB = commutator(op_A, VB)  # [σA, VB] should be zero (different subspaces)
            comm_B_VA = commutator(op_B, VA)  # [σB, VA] should be zero (different subspaces)

            @test iszero(simplify(comm_A_VB))
            @test iszero(simplify(comm_B_VA))

            # Self-commutators should not be zero
            comm_A_VA = commutator(op_A, VA)
            @test !iszero(simplify(comm_A_VA))
        end

        @testset "LaTeX Rendering Fix" begin
            # Test that the LaTeX rendering doesn't produce parse errors
            doublesum = Σ(σB(:r, :r, j) * σB(:r, :r, k), j, k)

            # The LaTeX representation should not contain problematic characters
            latex_str = repr(MIME"text/latex"(), doublesum)
            @test !contains(latex_str, "≠")  # Should not contain the problematic ≠ character
            @test contains(latex_str, "neq") || contains(latex_str, "\\neq")  # Should contain proper LaTeX
        end
    end
end
