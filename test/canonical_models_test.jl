using SecondQuantizedAlgebra
using Test
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: QAdd, Index, simplify

# ============================================================================
# Hand-derived expected forms for the canonical quantum-optics models, asserted
# by exact `==` against the package output. Each expected form comes from a
# pen-and-paper derivation written in the comment above the test; the goal is
# an independent reference that catches math errors in operations on indexed
# sums (Σ-construction, *, +, commutator, range-factor handling).
# ============================================================================

@testset "Canonical model commutators" begin

    # Shared Tavis-Cummings setup: a single bosonic cavity coupled to N
    # two-level atoms via dipole coupling g_i.
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    @variables N Δ
    g(k) = IndexedVariable(:g, k)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    H_TC = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    @testset "Tavis-Cummings: [H, a] = Δa - Σ_i g(i) σ_i^{12}" begin
        # Derivation: [a'a, a] = -a. [Σg(i)(a'σ12_i + aσ21_i), a]:
        # [a'σ12_i, a] = [a',a] σ12_i = -σ12_i (cross-space commutators trivial).
        # [aσ21_i, a] = 0. Sum: -Σ_i g(i) σ12_i.
        expected = Δ * a - Σ(g(i) * σ(1, 2, i), i)
        @test commutator(H_TC, a) == expected
    end

    @testset "Tavis-Cummings: [H, a'] = -Δa' + Σ_i g(i) σ_i^{21}" begin
        # Adjoint of the [H, a] result with H Hermitian.
        expected = -Δ * a' + Σ(g(i) * σ(2, 1, i), i)
        @test commutator(H_TC, a') == expected
    end

    @testset "Tavis-Cummings: [H, a'a] = Σ_i g(i)(a σ_i^{21} - a' σ_i^{12})" begin
        # [a'a, a'a] = 0. For each i:
        # [g(i) a' σ12_i, a'a] = g(i) [a', a'a] σ12_i = -g(i) a' σ12_i.
        # [g(i) a σ21_i, a'a] = g(i) [a, a'a] σ21_i = g(i) a σ21_i.
        expected = Σ(g(i) * a * σ(2, 1, i) - g(i) * a' * σ(1, 2, i), i)
        @test commutator(H_TC, a' * a) == expected
    end

    @testset "Tavis-Cummings: [H, σ_j^{22}] = g(j) a' σ_j^{12} - g(j) a σ_j^{21}" begin
        # Atom-j population EOM. Only the i=j term in the sum survives:
        # [σ12_j, σ22_j] = σ12_j (since σ12 σ22 = σ12 and σ22 σ12 = 0).
        # [σ21_j, σ22_j] = -σ21_j.
        expected = g(j) * a' * σ(1, 2, j) - g(j) * a * σ(2, 1, j)
        @test commutator(H_TC, σ(2, 2, j)) == expected
    end

    @testset "Tavis-Cummings: [H, σ_j^{12}] = g(j) a (σ_j^{22} - σ_j^{11})" begin
        # The classic optical-Bloch result for the atomic coherence EOM.
        expected = g(j) * a * σ(2, 2, j) - g(j) * a * σ(1, 1, j)
        @test commutator(H_TC, σ(1, 2, j)) == expected
    end

    @testset "Tavis-Cummings: [H, σ_j^{12} σ_k^{21}] merges multi-σ_j and σ_k" begin
        # Two-atom-correlator EOM with j ≠ k. The Σ_i collapses on i = j and
        # i = k diagonals; on each diagonal two same-site Transitions (e.g.
        # σ_j^{21} from the substitution adjacent to σ_j^{12} from X) must
        # compose. Regression for the bug where `σ_j^{12} σ_k^{21} σ_j^{21}`
        # was left un-merged (the two σ_j operators were separated by σ_k
        # and the canonical sort treated j, k as Undetermined).
        # Pen-and-paper: [a' σ_i^{12}, σ_j^{12} σ_k^{21}] → i=k gives
        # σ_j^{12}(σ_k^{11} - σ_k^{22}); [a σ_i^{21}, σ_j^{12} σ_k^{21}] →
        # i=j gives (σ_j^{22} - σ_j^{11}) σ_k^{21}. Off-diagonal vanishes
        # because i ≠ j, k commutes with both factors of X.
        k = Index(h, :k, N, ha)
        expected = g(j) * a * (σ(2, 2, j) - σ(1, 1, j)) * σ(2, 1, k) +
            g(k) * a' * σ(1, 2, j) * (σ(1, 1, k) - σ(2, 2, k))
        # The package result carries an implicit j ≠ k constraint from the
        # diagonal partition; assume_distinct_index aligns the expected
        # term-key with that ne so == compares like-for-like.
        @test assume_distinct_index(commutator(H_TC, σ(1, 2, j) * σ(2, 1, k)), [(j, k)]) ==
            assume_distinct_index(expected, [(j, k)])

        # Structural check independent of completeness rewriting: the raw
        # commutator must reduce to four 3-operator terms (no residual
        # 4-operator terms with un-merged same-index Transitions).
        raw = commutator(H_TC, σ(1, 2, j) * σ(2, 1, k))
        @test length(raw.arguments) == 4
        for term in keys(raw.arguments)
            @test length(term.ops) == 3
        end
    end

end

@testset "Dicke model commutators" begin

    # Single cavity + single collective spin (Pauli realization, 2-level atoms).
    hcD = FockSpace(:cavity)
    hsD = SpinSpace(:S)
    hD = hcD ⊗ hsD
    @qnumbers a::Destroy(hD, 1)
    Sx = Spin(hD, :S, 1, 2)
    Sy = Spin(hD, :S, 2, 2)
    Sz = Spin(hD, :S, 3, 2)
    @variables ω₀ ωₐ λ

    H_dicke = ω₀ * a' * a + ωₐ * Sz + λ * (a + a') * Sx

    @testset "Dicke: [H, a] = -ω₀ a - λ Sx" begin
        # [a'a, a] = -a. [Sz, a] = 0 (different subspaces). [(a+a')Sx, a]:
        # (a+a')[Sx,a] + [a+a',a]Sx = 0 + (-1)Sx = -Sx.
        expected = -ω₀ * a - λ * Sx
        @test commutator(H_dicke, a) == expected
    end

    @testset "Dicke: [H, Sx] = i ωₐ Sy" begin
        # [Sz, Sx] = i Sy. The cavity and coupling terms commute with Sx.
        expected = im * ωₐ * Sy
        @test commutator(H_dicke, Sx) == expected
    end

    @testset "Dicke: [H, Sy] = -i ωₐ Sx + i λ (a + a') Sz" begin
        # [Sz, Sy] = -i Sx. [(a+a')Sx, Sy] = (a+a')[Sx, Sy] = i (a+a') Sz.
        expected = -im * ωₐ * Sx + im * λ * (a + a') * Sz
        @test commutator(H_dicke, Sy) == expected
    end

    @testset "Dicke: [H, Sz] = i λ (a + a') Sy" begin
        # [(a+a')Sx, Sz] = (a+a')[Sx, Sz] = -i (a+a') Sy.
        # Note: [Sx, Sz] = -i Sy, so [Sz, Sx] = i Sy; here we want [Sx, Sz].
        expected = -im * λ * (a + a') * Sy
        @test commutator(H_dicke, Sz) == expected
    end

end

@testset "Sum identities" begin

    hf = FockSpace(:f)
    @qnumbers a::Destroy(hf)
    @variables N
    i = Index(hf, :i, N, hf)
    j = Index(hf, :j, N, hf)
    k = Index(hf, :k, N, hf)
    ai = IndexedOperator(a, i)
    aj = IndexedOperator(a, j)
    ai_dag = IndexedOperator(a', i)
    aj_dag = IndexedOperator(a', j)

    @testset "[Σ_i a_i, Σ_j a'_j] = N (canonical commutator)" begin
        c = commutator(Σ(ai, i), Σ(aj_dag, j))
        @test iszero(simplify(c - Complex(N)))
    end

    @testset "Σ_i (a_i a'_i) = N + Σ_i a'_i a_i" begin
        # The +1 from the bosonic commutator, summed over i, gives N.
        @test Σ(ai * ai_dag, i) == Σ(ai_dag * ai, i) + Complex(N)
    end

    @testset "Σ_i (a_i + c)^2 = Σ_i a_i^2 + 2c Σ_i a_i + N c^2 (Fock, c-number)" begin
        @variables c_sym
        expected = Σ(ai * ai, i) + 2 * c_sym * Σ(ai, i) + Complex(N * c_sym^2)
        @test Σ((ai + c_sym)^2, i) == expected
    end

    @testset "Pauli: [Σ_i σx_i, Σ_j σy_j] = 2i Σ_k σz_k" begin
        hP = PauliSpace(:p)
        iP = Index(hP, :ip, N, hP)
        jP = Index(hP, :jp, N, hP)
        σxi = IndexedOperator(Pauli(hP, :σ, 1), iP)
        σyj = IndexedOperator(Pauli(hP, :σ, 2), jP)
        σzi = IndexedOperator(Pauli(hP, :σ, 3), iP)
        # Off-diagonal vanishes (Pauli on different sites commute);
        # diagonal contributes 2i σz_i for each i.
        # Bound-index name will be `ip` (the canonically-earlier name survives the pin).
        expected = Σ(2im * σzi, iP)
        @test commutator(Σ(σxi, iP), Σ(σyj, jP)) == expected
    end

    @testset "Σ_i Σ_j (i≠j) g(i,j) σ_i σ_j builds a clean double sum" begin
        # Double-indexed prefactor and two-site product, no diagonal collapse.
        ha = NLevelSpace(:atom, 2)
        iN = Index(ha, :ia, N, ha)
        jN = Index(ha, :ja, N, ha)
        @variables g_dd
        Γ(p, q) = DoubleIndexedVariable(:Γ, p, q; identical = false)
        σ12_i = IndexedOperator(Transition(ha, :σ, 1, 2), iN)
        σ21_j = IndexedOperator(Transition(ha, :σ, 2, 1), jN)
        r = Σ(Σ(Γ(iN, jN) * σ12_i * σ21_j, jN), iN)
        # Two bound indices preserved; one off-diagonal term plus one pinned diag.
        @test iN in r.indices
        @test jN in r.indices
    end

end

@testset "expand_completeness × indexed sums" begin

    @variables N
    ha2 = NLevelSpace(:atom2, 2)
    ha3 = NLevelSpace(:atom3, 3)
    i2 = Index(ha2, :i2, N, ha2)
    i3 = Index(ha3, :i3, N, ha3)
    σ11_2(k) = IndexedOperator(Transition(ha2, :σ, 1, 1), k)
    σ22_2(k) = IndexedOperator(Transition(ha2, :σ, 2, 2), k)
    σ11_3(k) = IndexedOperator(Transition(ha3, :σ, 1, 1), k)
    σ22_3(k) = IndexedOperator(Transition(ha3, :σ, 2, 2), k)
    σ33_3(k) = IndexedOperator(Transition(ha3, :σ, 3, 3), k)

    @testset "Σ_i σ_i^gg = N - Σ_i (sum over excited)" begin
        # 2-level: σ^11 = 1 - σ^22, so Σ_i σ_i^11 = N - Σ_i σ_i^22.
        # The +1 inside Σ_i picks up the range factor N.
        @test expand_completeness(Σ(σ11_2(i2), i2)) ==
            Complex(N) - Σ(σ22_2(i2), i2)

        # 3-level: σ^11 = 1 - σ^22 - σ^33.
        @test expand_completeness(Σ(σ11_3(i3), i3)) ==
            Complex(N) - Σ(σ22_3(i3), i3) - Σ(σ33_3(i3), i3)
    end

    @testset "Σ_i σ_i^ee (non-GS) passes through" begin
        # σ^22 is not the GS in NLevelSpace(:atom2, 2) (GS defaults to 1).
        @test expand_completeness(Σ(σ22_2(i2), i2)) == Σ(σ22_2(i2), i2)
        @test expand_completeness(Σ(σ22_3(i3), i3)) == Σ(σ22_3(i3), i3)
        @test expand_completeness(Σ(σ33_3(i3), i3)) == Σ(σ33_3(i3), i3)
    end

    @testset "expand_completeness on bare projector keeps single-site form" begin
        # No indexed sum scope: σ^gg goes to (1 - σ^ee), no N factor.
        @test expand_completeness(σ11_2(i2)) == 1 - σ22_2(i2)
    end

end

@testset "average ∘ Σ roundtrips" begin

    @variables N
    hf = FockSpace(:f)
    @qnumbers a::Destroy(hf)
    i = Index(hf, :i, N, hf)
    j = Index(hf, :j, N, hf)
    ai = IndexedOperator(a, i)
    ai_dag = IndexedOperator(a', i)
    aj_dag = IndexedOperator(a', j)

    @testset "undo_average(average(x)) == x for indexed sums" begin
        s1 = Σ(ai_dag * ai, i)
        @test undo_average(average(s1)) == s1

        # Sum with constant residual baked in (from a*a' = a'a + 1, scaled by N).
        s2 = Σ(ai * ai_dag, i)
        @test undo_average(average(s2)) == s2

        # Mixed: indexed sum plus standalone constant.
        s3 = s1 + 3
        @test undo_average(average(s3)) == s3

        # Distinct-index sum (two separate scopes in one QAdd).
        s4 = Σ(ai, i) + Σ(aj_dag, j)
        @test undo_average(average(s4)) == s4

        # Double-indexed product with diagonal pin and N residual.
        s5 = Σ(ai, i) * Σ(aj_dag, j)
        @test undo_average(average(s5)) == s5
    end

    @testset "average of indexed sum is linear in the sum body" begin
        lhs = average(2 * Σ(ai_dag * ai, i))
        rhs = 2 * average(Σ(ai_dag * ai, i))
        @test undo_average(lhs) == undo_average(rhs)
    end

end

@testset "Lindblad dissipator with indexed jump operators" begin

    # D[J]ρ = JρJ' - ½ {J'J, ρ}.
    # Heisenberg adjoint: D'[J](x) = J' x J - ½ {J' J, x}.

    @variables N γ
    ha = NLevelSpace(:atom, 2)
    σ(α, β, k) = IndexedOperator(Transition(ha, :σ, α, β), k)
    j = Index(ha, :j, N, ha)

    function dissipator_adjoint(J, J_dag, x)
        return J_dag * x * J - (1 // 2) * (J_dag * J * x + x * J_dag * J)
    end

    @testset "Single-site jump: D'[σ^{12}](σ^{22}) = -σ^{22}" begin
        J = σ(1, 2, j)
        Jd = σ(2, 1, j)
        x = σ(2, 2, j)
        # σ21 σ22 σ12 = σ21 · 0 = 0 (since σ22 σ12 = 0).
        # σ21 σ12 σ22 = σ22 σ22 = σ22; σ22 σ21 σ12 = σ22 σ22 = σ22.
        # D'(σ22) = 0 - (σ22 + σ22)/2 = -σ22.
        @test dissipator_adjoint(J, Jd, x) == -σ(2, 2, j)
    end

    @testset "Single-site jump: D'[σ^{12}](σ^{11}) = +σ^{22}" begin
        J = σ(1, 2, j)
        Jd = σ(2, 1, j)
        x = σ(1, 1, j)
        # σ21 σ11 σ12 = σ21 σ12 = σ22.
        # σ21 σ12 σ11 = σ22 σ11 = 0; σ11 σ21 σ12 = 0.
        # D'(σ11) = σ22.
        @test dissipator_adjoint(J, Jd, x) == 1 * σ(2, 2, j)
    end

    @testset "Collective decay: well-formedness and j-dependence" begin
        # J = Σ_i σ_i^{12}, J' = Σ_k σ_k^{21}, x = σ_j^{22}.
        i = Index(ha, :i, N, ha)
        k = Index(ha, :k, N, ha)
        Jsum = Σ(σ(1, 2, i), i)
        Jsum_dag = Σ(σ(2, 1, k), k)
        x = σ(2, 2, j)
        res = dissipator_adjoint(Jsum, Jsum_dag, x)
        @test !iszero(res)
        @test res isa SecondQuantizedAlgebra.QAdd
        has_j = false
        for (term, _) in res.arguments
            for op in term.ops
                op.index == j && (has_j = true; break)
            end
            has_j && break
        end
        @test has_j
    end

end
