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
        result4 = 1im * commutator(H, op4)
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
end
