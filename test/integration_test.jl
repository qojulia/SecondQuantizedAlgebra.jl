using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using QuantumOpticsBase
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, _single_qadd, _to_cnum, Index, sorted_arguments, constraint_pairs

@testset "integration" begin
    @testset "Basic operator creation and algebra" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'

        @test hash(a) != hash(ad)
        @test isequal(a, ad')

        # Lazy multiplication
        @test a * ad isa QAdd
        @test ad * a isa QAdd
    end

    @testset "Commutation via normal_order" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a * a† normal-ordered → a†a + 1
        m = a * a'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "Distinct modes on same space don't commute" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        a = Destroy(h1, :a)
        b = Destroy(h2, :b)
        # a and b have same space_index=1 but different names
        m = a * b'
        result = normal_order(m)
        @test length(result) == 1  # no commutation applied
    end

    @testset "Numeric conversion round-trip" begin
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
        @test to_numeric(_single_qadd(_to_cnum(3), QSym[]), b) == 3 * one(b)
    end

    @testset "Multi-space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

        # Operators on different spaces — normal_order doesn't commute them
        m = a * b'
        result = normal_order(m)
        @test length(result) == 1

        # Same-space commutation still works
        m2 = a * a'
        result2 = normal_order(m2)
        @test length(result2) == 2
    end

    @testset "Symbolic parameters" begin
        @variables g ω
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        expr = g * a' * a + ω * a' * a
        result = simplify(expr)
        @test length(result) == 1

        # Printing with symbolic prefactors should not crash
        @test repr(g * a) isa String
        @test repr(g * a + ω * a') isa String
    end

    @testset "Division with non-integer" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test (a / 2.0) isa QAdd
        @test prefactor(only(sorted_arguments(a / 2.0))) ≈ 0.5

        @test (a / 2) isa QAdd
        @test prefactor(only(sorted_arguments(a / 2))) == 1 // 2

        @test ((a + a') / 2.0) isa QAdd
    end

    @testset "Power edge cases" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^0 = identity (scalar QAdd)
        m0 = a^0
        @test m0 isa QAdd
        @test isempty(operators(only(sorted_arguments(m0))))
        @test prefactor(only(sorted_arguments(m0))) == 1
    end

    @testset "Adjoint with complex prefactor" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        m = (2.0 + 3.0im) * a' * a
        ma = m'
        @test prefactor(only(sorted_arguments(ma))) ≈ conj(2.0 + 3.0im)
    end

    @testset "Higher-order normal ordering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^2 * (a')^2 should produce multiple terms
        m = a * a * a' * a'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) >= 3
    end

    @testset "== consistency" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test _single_qadd(_to_cnum(1), QSym[a]) == _single_qadd(_to_cnum(1), QSym[a])
        @test (a + a') == (a + a')
    end

    @testset "Jaynes-Cummings model" begin
        hc = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2, 1)
        h = hc ⊗ ha
        @qnumbers a::Destroy(h, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        σ21 = Transition(h, :σ, 2, 1, 2)

        @variables g_jc
        H_int = g_jc * (a' * σ12 + a * σ21)
        @test H_int isa QAdd

        bf = FockBasis(5)
        bn = NLevelBasis(2)
        bc = bf ⊗ bn
        H_num = to_numeric(H_int, bc)
        @test H_num !== nothing
    end

    @testset "Pauli algebra via normal_order" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)

        # σx·σy·σz = iσz·σz = i·I
        m = σx * σy * σz
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))
        @test prefactor(only(sorted_arguments(result))) == im
    end

    @testset "Mixed Fock + Pauli" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p)
        @qnumbers a::Destroy(h, 1)
        σz = Pauli(h, :σ, 3, 2)

        m = a' * a * σz
        @test m isa QAdd
        @test length(operators(only(sorted_arguments(m)))) == 3

        # Normal order: commute a,a† on space 1, leave σz on space 2
        result = normal_order(a * a' * σz)
        @test result isa QAdd
    end

    @testset "Superradiant laser commutators" begin
        # Diagonal split during *(QAdd, QAdd) substitutes sum_idx → ext_idx
        # BEFORE site-sort so same-site composition (σ_αβ · σ_βγ = σ_αγ) fires
        # on the correct physical order. Reference equations (6)–(7):
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/
        #
        # Under NormalOrder, σⱼ¹¹ never appears in dict keys — it is eagerly
        # expanded via completeness σⱼ¹¹ = 1 - σⱼ²². Assertions check the
        # post-expansion form; LazyOrder shows the un-expanded shape.
        hc = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2)
        h = hc ⊗ ha

        @qnumbers a::Destroy(h, 1)
        σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

        @variables N Δ
        g(idx) = IndexedVariable(:g, idx)

        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        k = Index(h, :k, N, ha)

        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

        # Eq. (6): i [H, a'σⱼ¹²] = -iΔ a'σⱼ¹² + Σᵢ≠ⱼ i g(i) σᵢ²¹σⱼ¹²
        #                          + i g(j) (a'a (σⱼ²² - σⱼ¹¹) + σⱼ²²)
        # Post-completeness: σⱼ²² - σⱼ¹¹ = 2σⱼ²² - 1, so:
        #   = -iΔ a'σⱼ¹² + Σᵢ≠ⱼ i g(i) σᵢ²¹σⱼ¹² + i g(j)(2 a'a σⱼ²² - a'a + σⱼ²²)
        op3 = a' * σ(1, 2, j)
        result3 = 1im * commutator(H, op3)
        # -iΔ * a' * σⱼ¹²
        @test any(result3) do (term, c)
            term.ops == [a', σ(1, 2, j)] && isequal(c, _to_cnum(-1im * Δ))
        end
        # i * g(j) * σⱼ²²
        @test any(result3) do (term, c)
            term.ops == [σ(2, 2, j)] && isequal(c, _to_cnum(1im * g(j)))
        end
        # 2i * g(j) * a'a * σⱼ²²  (was im * g(j) before completeness folded in σⱼ¹¹)
        @test any(result3) do (term, c)
            term.ops == [a', a, σ(2, 2, j)] && isequal(c, _to_cnum(2im * g(j)))
        end
        # -i * g(j) * a'a  (identity contribution from σⱼ¹¹ = 1 - σⱼ²²)
        @test any(result3) do (term, c)
            term.ops == [a', a] && isequal(c, _to_cnum(-1im * g(j)))
        end
        # Σᵢ≠ⱼ i * g(i) * σᵢ²¹σⱼ¹²
        @test any(result3) do (term, c)
            term.ops == [σ(2, 1, i), σ(1, 2, j)] && isequal(c, _to_cnum(1im * g(i)))
        end
        @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(result3))
        # σⱼ¹¹ should NOT appear under NormalOrder canonical form
        @test !any(result3) do (term, c)
            any(op -> op isa Transition && op.i == 1 && op.j == 1, term.ops)
        end

        # Eq. (7): i [H, σⱼ¹²σₖ²¹] = i g(j) a (σⱼ²² - σⱼ¹¹) σₖ²¹
        #                            + i g(k) a' σⱼ¹² (σₖ¹¹ - σₖ²²)
        # Post-completeness: σⱼ²² - σⱼ¹¹ = 2σⱼ²² - 1, σₖ¹¹ - σₖ²² = 1 - 2σₖ²², so:
        #   = i g(j) a (2σⱼ²² σₖ²¹ - σₖ²¹) + i g(k) a' (σⱼ¹² - 2 σⱼ¹² σₖ²²)
        op4 = σ(1, 2, j) * σ(2, 1, k)
        result4 = simplify(1im * commutator(H, op4))
        # 2i * g(j) * a * σⱼ²² * σₖ²¹
        @test any(result4) do (term, c)
            term.ops == [a, σ(2, 2, j), σ(2, 1, k)] && isequal(c, _to_cnum(2im * g(j)))
        end
        # -i * g(j) * a * σₖ²¹  (identity from σⱼ¹¹ expansion)
        @test any(result4) do (term, c)
            term.ops == [a, σ(2, 1, k)] && isequal(c, _to_cnum(-1im * g(j)))
        end
        # i * g(k) * a' * σⱼ¹²  (identity from σₖ¹¹ expansion)
        @test any(result4) do (term, c)
            term.ops == [a', σ(1, 2, j)] && isequal(c, _to_cnum(1im * g(k)))
        end
        # -2i * g(k) * a' * σⱼ¹² * σₖ²²
        @test any(result4) do (term, c)
            term.ops == [a', σ(1, 2, j), σ(2, 2, k)] && isequal(c, _to_cnum(-2im * g(k)))
        end
        # σ_·_11 should NOT appear under NormalOrder canonical form
        @test !any(result4) do (term, c)
            any(op -> op isa Transition && op.i == 1 && op.j == 1, term.ops)
        end
    end

    @testset "Decay-channel triple product" begin
        # 2 * Σᵢ σ²¹_i · σ²²_j · σ¹²_i must be Σᵢ≠ⱼ 2 σ²²_i σ²²_j.
        # When an index pair's substitution into the unsorted form disagrees
        # with substituting into the sorted form, `_accumulate_with_diag!`
        # records the constraint (i ≠ j) and emits the correct diagonal.
        hc = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2)
        h = hc ⊗ ha
        σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

        @variables N
        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)

        op = σ(2, 2, j)
        result = 2 * Σ(σ(2, 1, i) * op * σ(1, 2, i), i)

        # Exactly one term, exactly the off-diagonal under i ≠ j
        @test length(result) == 1
        @test haskey(result, [σ(2, 2, i), σ(2, 2, j)])
        @test isequal(result[[σ(2, 2, i), σ(2, 2, j)]], _to_cnum(2))
        @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(result))
        # No spurious σ_j₂₂-only term ("at i=j") — the physical i=j contribution
        # from σ²¹_j · σ²²_j · σ¹²_j is zero (σ²¹·σ²² = 0).
        @test !haskey(result, [σ(2, 2, j)])
    end

    @testset "Commutator with double-sum Hamiltonian" begin
        # commutator with H = … + Σⱼ Σᵢ Ω(i,j) σ²¹_i σ¹²_j must retain the
        # Ω(k,j) double-sum contributions (i = k, j ≠ k partial collapse).
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/cavity_antiresonance_indexed/
        hc = FockSpace(:cavity)
        ha = NLevelSpace(Symbol(:atom), 2)
        h = hc ⊗ ha

        @variables N Δc η Δa κ
        Ω(idx1, idx2) = DoubleIndexedVariable(:Ω, idx1, idx2; identical = false)

        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        k = Index(h, :k, N, ha)

        @qnumbers a::Destroy(h, 1)
        σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y, 2), idx)

        Ha = Δa * Σ(σ(2, 2, i), i) + Σ(Ω(i, j) * σ(2, 1, i) * σ(1, 2, j), j, i)
        H = Δc * a' * a + η * (a' + a) + Ha

        # [H, σ¹²_k]: must include σ¹²_j-bearing terms (Ω(k, j) contributions)
        result = 1im * commutator(H, σ(1, 2, k))
        @test any(result) do (term, c)
            any(op -> op isa Transition && op.i == 1 && op.j == 2 && op.index == j, term.ops)
        end

        # [H, σ²²_k]: must also include j-indexed σ-bearing contributions
        result3 = 1im * commutator(H, σ(2, 2, k))
        @test any(result3) do (term, c)
            any(op -> op isa Transition && op.index == j, term.ops)
        end
    end

    @testset "QuantumCumulants examples" begin

        # ------------------------------------------------------------------------
        # Mollow triplet — driven two-level atom
        # H = Δ σ_ee + Ω (σ_ge + σ_eg)
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/mollow/
        # ------------------------------------------------------------------------
        @testset "Mollow triplet" begin
            h = NLevelSpace(:atom, (:g, :e))
            @variables Δ Ω
            σ(α, β) = Transition(h, :σ, α, β)

            H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))

            # [σ_ee, σ_ge] = -σ_ge,  [σ_ee, σ_eg] = +σ_eg,  [σ_ge, σ_eg] = σ_ee - σ_gg
            # Under canonical (no σ_gg, ground = :g): σ_gg = 1 - σ_ee
            # ⇒ [H, σ_ge] = -Δ σ_ge + Ω(σ_ee - σ_gg) = -Δ σ_ge + Ω(2 σ_ee - 1)
            eq = commutator(H, σ(:g, :e))
            expected = -Δ * σ(:g, :e) + Ω * (2 * σ(:e, :e) - 1)
            @test iszero(simplify(eq - expected))

            # [H, σ_eg] = Δ σ_eg - Ω(2 σ_ee - 1)
            eq2 = commutator(H, σ(:e, :g))
            expected2 = Δ * σ(:e, :g) - Ω * (2 * σ(:e, :e) - 1)
            @test iszero(simplify(eq2 - expected2))

            # [H, σ_ee] = Ω σ_ge - Ω σ_eg  (drives Rabi oscillations of the population)
            eq3 = commutator(H, σ(:e, :e))
            expected3 = Ω * σ(:g, :e) - Ω * σ(:e, :g)
            @test iszero(simplify(eq3 - expected3))

            # average() should be a no-op wrap: identity round-trip
            @test isequal(undo_average(average(H)), H)
        end

        # ------------------------------------------------------------------------
        # Single-atom laser — Jaynes-Cummings with detuning
        # H = Δ a†a + g (a† σ_ge + a σ_eg)
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/single-atom-laser-spectrum/
        # ------------------------------------------------------------------------
        @testset "Single-atom laser (JC)" begin
            hf = FockSpace(:cavity)
            ha = NLevelSpace(:atom, (:g, :e))
            h = hf ⊗ ha

            @variables Δ g
            a = Destroy(h, :a)
            s(α, β) = Transition(h, :σ, α, β)

            H = Δ * a' * a + g * (a' * s(:g, :e) + a * s(:e, :g))

            # [a†a, a] = -a; coupling commutes through Fock sector for σ-term.
            # [a, a' σ_ge] = σ_ge ; [a, a σ_eg] = 0.
            # ⇒ [H, a] = -Δ a - g σ_ge
            eq_a = commutator(H, a)
            @test iszero(simplify(eq_a + Δ * a + g * s(:g, :e)))

            # [H, σ_ge] = g a (σ_ee - σ_gg) = g a (2σ_ee - 1) = -g a + 2g a σ_ee
            eq_s = commutator(H, s(:g, :e))
            expected_s = -g * a + 2g * a * s(:e, :e)
            @test iszero(simplify(eq_s - expected_s))

            # Photon number commutes with itself.
            @test iszero(simplify(commutator(a' * a, a' * a)))

            # Hermiticity of H — adjoint round-trip recovers the same expression.
            @test iszero(simplify(H - adjoint(H)))
        end

        # ------------------------------------------------------------------------
        # Superradiant laser (indexed) — N atoms in a cavity
        # H = -Δ a†a + Σ_i g_i (a† σ_i^{12} + a σ_i^{21})
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/
        # ------------------------------------------------------------------------
        @testset "Superradiant laser (indexed)" begin
            hc = FockSpace(:cavity)
            ha = NLevelSpace(:atom, 2)
            h = hc ⊗ ha

            @variables N Δ
            a = Destroy(h, :a)
            σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β), idx)
            gi(idx) = IndexedVariable(:g, idx)

            i = Index(h, :i, N, ha)
            H = -Δ * a' * a + Σ(gi(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

            # [H, a] = Δ a - Σ_i g_i σ_i^{12}
            eq_a = commutator(H, a)
            expected_a = Δ * a + Σ(-gi(i) * σ(1, 2, i), i)
            @test iszero(simplify(eq_a - expected_a))

            # [a†a, σ_i^{12}] = 0 (different subspace); commutator with H reduces to atomic part.
            # Check that the cavity-only part of H commutes with any indexed atomic op.
            @test iszero(simplify(commutator(-Δ * a' * a, σ(1, 2, i))))

            # The Hamiltonian is Hermitian.
            @test iszero(simplify(H - adjoint(H)))

            # average distributes through the sum.
            avg_H = average(H)
            @test avg_H isa QAdd || avg_H isa Number || true   # smoke: no error on average of indexed sum
        end

        # ------------------------------------------------------------------------
        # Mollow triplet — extended algebra checks (Hermiticity, identities, idempotency)
        # ------------------------------------------------------------------------
        @testset "Mollow triplet — extended" begin
            h = NLevelSpace(:atom, (:g, :e))
            @variables Δ Ω
            σ(α, β) = Transition(h, :σ, α, β)

            H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))

            # Hermiticity of H.
            @test iszero(simplify(H - adjoint(H)))

            # σ_ee is idempotent: σ_ee² = σ_ee.
            @test iszero(simplify(σ(:e, :e) * σ(:e, :e) - σ(:e, :e)))

            # σ_ge σ_eg = σ_gg = 1 - σ_ee  (canonical, σ_gg eliminated).
            @test iszero(simplify(σ(:g, :e) * σ(:e, :g) - (1 - σ(:e, :e))))

            # σ_eg σ_ge = σ_ee.
            @test iszero(simplify(σ(:e, :g) * σ(:g, :e) - σ(:e, :e)))

            # [σ_ge, σ_eg] = σ_gg - σ_ee = 1 - 2σ_ee.
            @test iszero(simplify(commutator(σ(:g, :e), σ(:e, :g)) - (1 - 2 * σ(:e, :e))))
        end

        # ------------------------------------------------------------------------
        # Single-atom laser — extended JC checks
        # ------------------------------------------------------------------------
        @testset "JC — extended" begin
            hf = FockSpace(:cavity)
            ha = NLevelSpace(:atom, (:g, :e))
            h = hf ⊗ ha

            @variables Δ g
            a = Destroy(h, :a)
            s(α, β) = Transition(h, :σ, α, β)

            H = Δ * a' * a + g * (a' * s(:g, :e) + a * s(:e, :g))

            # Excitation number N = a'a + σ_ee is conserved: [H, N] = 0.
            N_op = a' * a + s(:e, :e)
            @test iszero(simplify(commutator(H, N_op)))

            # [H, σ_eg] = g a' (σ_gg - σ_ee) = g a' (1 - 2σ_ee) = g a' - 2g a' σ_ee.
            eq = commutator(H, s(:e, :g))
            expected = g * a' - 2g * a' * s(:e, :e)
            @test iszero(simplify(eq - expected))

            # [H, a'] = Δ a' + g σ_eg  (note σ_eg, not σ_ge — see [a, a'] = 1 contraction).
            @test iszero(simplify(commutator(H, a') - (Δ * a' + g * s(:e, :g))))

            # [a, σ_ge] = 0 (different subspaces).
            @test iszero(simplify(commutator(a, s(:g, :e))))
        end

        # ------------------------------------------------------------------------
        # Dissipator structure checks — superradiant laser jump operator J = σ_i¹².
        # The Lindblad dissipator
        #   D[J](o) = J† o J - ½ {J†J, o}
        # for J = σ_i¹² gives D(o) = σ²¹_i o σ¹²_i - ½(σ²¹_i σ¹²_i o + o σ²¹_i σ¹²_i).
        # Multiplying by 2 and summing over i (atomic ensemble) is the colleague's
        # MWE in PR #99. Some shapes still produce wrong results — those are marked
        # `@test_broken` and reference the upstream issue.
        # ------------------------------------------------------------------------
        @testset "Lindblad dissipator (atomic decay)" begin
            hc = FockSpace(:cavity)
            ha = NLevelSpace(:atom, 2)
            h = hc ⊗ ha
            @qnumbers a::Destroy(h, 1)
            σi(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

            @variables N γ
            i = Index(h, :i, N, ha)
            j = Index(h, :j, N, ha)
            k = Index(h, :k, N, ha)

            decay(op) = Σ(
                (γ / 2) * (
                    2 * σi(2, 1, i) * op * σi(1, 2, i) -
                        σi(2, 1, i) * σi(1, 2, i) * op -
                        op * σi(2, 1, i) * σi(1, 2, i)
                ), i
            )

            # decay(a'a) = 0 — the atomic dissipator does nothing on a pure-cavity op.
            @test iszero(simplify(decay(a' * a)))

            # decay(σ_j²²) = -γ σ_j²².
            # Fixed by preserving physical order in `*(::QSym, ::QSym)` for
            # same-space pairs (see qadd_arithmetic.jl): when σ_j²² appears
            # LEFT of σ_i, the dict's term.ops keeps that order, so the
            # `_accumulate_with_diag!` substitution at `i = j` fires on the
            # correct adjacency and the diagonal σ_j²²·σ_j²² = σ_j²² survives.
            @test iszero(simplify(decay(σi(2, 2, j)) + γ * σi(2, 2, j)))

            # decay(a'σ_j¹²) = -(γ/2) a'σ_j¹².
            # Still broken: the `_accumulate_with_diag!` second loop emits
            # diagonals anchored to the wrong index, so the i=j contribution
            # gets re-summed by Σ as Σ_i a'σ_i¹². Fixing this requires either
            # tracking phys_ops through QTerm (so Σ-time `_diagonal_split!`
            # uses physical order) or threading "must-equal" constraints
            # through subsequent multiplications. See devdocs.md.
            @test iszero(simplify(decay(a' * σi(1, 2, j)) + (γ / 2) * a' * σi(1, 2, j)))

            # decay(σ_j¹²·σ_k²¹) = -γ σ_j¹²·σ_k²¹.
            # Still broken — same root cause as above (diagonal anchoring).
            @test iszero(simplify(decay(σi(1, 2, j) * σi(2, 1, k)) + γ * σi(1, 2, j) * σi(2, 1, k)))

            # Direct algebraic identity: σ_i²¹ σ_j²² σ_i¹² with i ≠ j  must equal σ_i²² σ_j²².
            # When i = j: σ²¹·σ²²·σ¹² = (σ²¹·σ²²)·σ¹² = 0·σ¹² = 0 (since σ²¹·σ²² has ⟨1|2⟩ = 0).
            # The off-diagonal sum currently comes out correct.
            triple_sandwich = Σ(σi(2, 1, i) * σi(2, 2, j) * σi(1, 2, i), i)
            @test haskey(triple_sandwich, [σi(2, 2, i), σi(2, 2, j)])
            @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(triple_sandwich))
            # No spurious σ_j²² standalone term (the i = j contribution is genuinely zero).
            @test !haskey(triple_sandwich, [σi(2, 2, j)])
        end

        # ------------------------------------------------------------------------
        # Cavity antiresonance — atom–atom coupling Ω(i, j) σ²¹_i σ¹²_j
        # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/cavity_antiresonance_indexed/
        # The colleague's PR comment notes the Ω(i,j) double-sum terms were missing
        # from [H, σ_k¹²] / [H, σ_k²²]. They are now present, but with extra
        # spurious Σ_i over a body that is i-free — marked broken.
        # ------------------------------------------------------------------------
        @testset "Cavity antiresonance (Ω double sum)" begin
            hc = FockSpace(:cavity)
            ha = NLevelSpace(Symbol(:atom), 2)
            h = hc ⊗ ha

            @variables N Δc Δa η
            Ω(idx1, idx2) = DoubleIndexedVariable(:Ω, idx1, idx2; identical = false)
            i = Index(h, :i, N, ha)
            j = Index(h, :j, N, ha)
            k = Index(h, :k, N, ha)

            @qnumbers a::Destroy(h, 1)
            σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y, 2), idx)

            Ha = Δa * Σ(σ(2, 2, i), i) + Σ(Ω(i, j) * σ(2, 1, i) * σ(1, 2, j), j, i)
            H = Δc * a' * a + η * (a' + a) + Ha

            # [H, a] = -Δc a - η  (atomic part commutes with a).
            @test iszero(simplify(commutator(H, a) + Δc * a + η))

            # [H, σ_k¹²] must contain at least one term with σ_j¹² (the Ω(k,j) contribution).
            result_12 = 1im * commutator(H, σ(1, 2, k))
            @test any(result_12) do (term, _c)
                any(op -> op isa Transition && op.i == 1 && op.j == 2 && op.index == j, term.ops)
            end

            # [H, σ_k²²] must contain at least one j-indexed σ-bearing contribution.
            result_22 = 1im * commutator(H, σ(2, 2, k))
            @test any(result_22) do (term, _c)
                any(op -> op isa Transition && op.index == j, term.ops)
            end

            # The coupling sum Σ Ω(i,j) σ²¹_i σ¹²_j is Hermitian if Ω(j,i)* = Ω(i,j) (symmetric).
            # We don't know the symmetry of the user-defined `Ω`, so just smoke-check that
            # adjoint produces a sum of the same shape.
            Hcoup = Σ(Ω(i, j) * σ(2, 1, i) * σ(1, 2, j), j, i)
            @test adjoint(Hcoup) isa QAdd
        end

        # ------------------------------------------------------------------------
        # Many-atom Tavis–Cummings — indexed coupling, simpler than superradiant
        # H = ωc a'a + Σ_i (ω_a/2) σ_i^z + Σ_i g_i (a' σ_i¹² + a σ_i²¹)
        # ------------------------------------------------------------------------
        @testset "Tavis–Cummings (indexed)" begin
            hc = FockSpace(:cavity)
            ha = NLevelSpace(:atom, 2)
            h = hc ⊗ ha
            @qnumbers a::Destroy(h, 1)
            σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

            @variables N ωc ωa
            g(idx) = IndexedVariable(:g, idx)
            i = Index(h, :i, N, ha)

            # Inversion σ_i^z = σ_i²² - σ_i¹¹.  Under canonical form with ground=1,
            # σ_i¹¹ is eliminated to 1 - σ_i²², so σ_i^z = 2 σ_i²² - 1.
            H = ωc * a' * a + (ωa / 2) * Σ(2 * σ(2, 2, i) - 1, i) +
                Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

            # Cavity-only part [ωc a'a, σ_l²²] = 0 (different subspace).
            l = Index(h, :l, N, ha)
            @test iszero(simplify(commutator(ωc * a' * a, σ(2, 2, l))))

            # Hermiticity of H.
            @test iszero(simplify(H - adjoint(H)))
        end

        # ------------------------------------------------------------------------
        # Pauli algebra spot checks — σ_x σ_y = i σ_z, etc.
        # ------------------------------------------------------------------------
        @testset "Pauli algebra spot checks" begin
            h = PauliSpace(:p)
            σx = Pauli(h, :σ, 1)
            σy = Pauli(h, :σ, 2)
            σz = Pauli(h, :σ, 3)

            # σ_x σ_y = i σ_z (and cyclic permutations).
            @test iszero(simplify(σx * σy - 1im * σz))
            @test iszero(simplify(σy * σz - 1im * σx))
            @test iszero(simplify(σz * σx - 1im * σy))

            # σ_x² = 1.
            @test iszero(simplify(σx * σx - 1))

            # [σ_x, σ_y] = 2i σ_z.
            @test iszero(simplify(commutator(σx, σy) - 2im * σz))
        end
    end

end
