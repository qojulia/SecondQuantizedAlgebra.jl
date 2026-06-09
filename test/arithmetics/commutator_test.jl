using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, _single_qadd, _to_cnum, Index,
    sorted_arguments, constraint_pairs

@testset "commutator" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Scalar short-circuits" begin
        @test iszero(commutator(1, 2))
        @test iszero(commutator(3, a))
        @test iszero(commutator(a, 5))
        @test iszero(commutator(2.0, ad))
    end

    @testset "Fock: [a, a†] = 1" begin
        result = commutator(a, ad)
        @test result isa QAdd
        # simplify(a*a' - a'*a) = 1 (scalar identity)
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))
        @test prefactor(only(sorted_arguments(result))) == 1
    end

    @testset "Fock: [a†, a] = -1" begin
        result = commutator(ad, a)
        @test result isa QAdd
        @test length(result) == 1
        @test prefactor(only(sorted_arguments(result))) == -1
    end

    @testset "Identical operators: [a, a] = 0" begin
        @test iszero(commutator(a, a))
    end

    @testset "Different spaces: [a, b] = 0" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1, a2))
        @test iszero(commutator(a1, a2'))
    end

    @testset "QAdd with QSym — shared space" begin
        # Number-operator identities: [a'a, a] = -a, [a'a, a'] = a'.
        # RHS wrapped via `1 *` to promote QSym to QAdd for equality.
        @test commutator(ad * a, a) == -a
        @test commutator(ad * a, ad) == 1 * ad
    end

    @testset "QAdd with QSym — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1' * a1, a2))
        @test iszero(commutator(a2, a1' * a1))
    end

    @testset "QAdd with QAdd — no shared space" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(commutator(a1' * a1, a2' * a2))
    end

    @testset "QAdd bilinearity: [a + a†, a] = [a,a] + [a†,a]" begin
        expr = a + ad
        result = commutator(expr, a)
        @test result isa QAdd
        # [a,a] = 0, [a†,a] = -1 → total = -1
        simplified = simplify(result)
        @test length(simplified) == 1
        @test prefactor(only(sorted_arguments(simplified))) == -1
    end

    @testset "QAdd, QAdd bilinearity" begin
        result = commutator(a + ad, a + ad)
        @test result isa QAdd
        # [a+a†, a+a†] = [a,a]+[a,a†]+[a†,a]+[a†,a†] = 0+1-1+0 = 0
        simplified = simplify(result)
        @test iszero(simplified)
    end

    @testset "Sum-sum commutator drops dead indices" begin
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        ai = IndexedOperator(a, i)
        adj = IndexedOperator(ad, j)

        s1 = Σ(ai, i)
        s2 = Σ(adj, j)
        result = commutator(s1, s2)
        @test result isa QAdd
        # All off-diagonal (i≠j) terms cancel between s1*s2 and s2*s1, and the
        # surviving diagonal contributions are c-numbers, so no surviving term
        # depends on i or j and neither index should remain in `.indices`.
        @test isempty(result.indices)
    end

    @testset "Nested commutator" begin
        # [a, [a, a†]] = [a, 1] = 0
        inner = commutator(a, ad)
        result = commutator(a, inner)
        @test result isa QAdd
        @test iszero(simplify(result))
    end

    @testset "Spin: [Sx, Sy] = iSz (and cyclic)" begin
        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test commutator(Sx, Sy) == im * Sz
        @test commutator(Sy, Sz) == im * Sx
        @test commutator(Sz, Sx) == im * Sy
        # Antisymmetry
        @test commutator(Sy, Sx) == -im * Sz
    end

    @testset "PhaseSpace: [X, P] = i" begin
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        result = commutator(x, p)
        # Canonical commutator is the c-number i (empty-ops QAdd with prefactor im).
        @test length(result) == 1
        (term, c) = only(collect(result))
        @test isempty(term.ops)
        @test c == im
        @test commutator(p, x) == -commutator(x, p)
    end

    @testset "Return type is always QAdd" begin
        @test commutator(a, ad) isa QAdd
        @test commutator(a, a) isa QAdd
        @test commutator(1, a) isa QAdd
        @test commutator(a + ad, a) isa QAdd
        @test commutator(a + ad, a + ad) isa QAdd
        @test commutator(ad * a, a) isa QAdd
    end

    @testset "Type stability" begin
        @inferred commutator(a, ad)
        @inferred commutator(ad, a)
        @inferred commutator(a, a)
        @inferred commutator(1, a)
        @inferred commutator(a, 1)
        @inferred commutator(ad * a, a)
        @inferred commutator(a, ad * a)
        @inferred commutator(ad * a, a * ad)
        @inferred commutator(a + ad, a)
        @inferred commutator(a, a + ad)
        @inferred commutator(a + ad, a + ad)

        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @inferred commutator(a1, a2)
        @inferred commutator(a1, a2')
    end

    @testset "Allocations" begin
        # Warmup
        commutator(a, ad)
        commutator(a, a)
        commutator(ad * a, a)
        commutator(a + ad, a)

        # Short-circuit (zero) should be very cheap
        @test @allocations(commutator(a, a)) < 100
        @test @allocations(commutator(1, a)) < 100

        # Non-trivial but simple
        @test @allocations(commutator(a, ad)) < 20000
        @test @allocations(commutator(ad * a, a)) < 30000

        # Bilinearity over QAdd
        @test @allocations(commutator(a + ad, a)) < 40000
    end

    @testset "anticommutator" begin
        # Fock: {a, a†} = a a† + a† a = (a†a + 1) + a†a = 1 + 2 a†a
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'
        r = anticommutator(a, ad)
        @test r isa QAdd
        @test isequal(r, 1 + 2 * ad * a)

        # {a, a} = 2 a²
        @test isequal(anticommutator(a, a), 2 * a * a)

        # Pauli: {σⱼ, σₖ} = 2 δⱼₖ
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test iszero(anticommutator(σx, σy))
        @test iszero(anticommutator(σx, σz))
        @test iszero(anticommutator(σy, σz))
        # σx² = I, so {σx, σx} = 2
        @test isequal(simplify(anticommutator(σx, σx)), _single_qadd(_to_cnum(2), QSym[]))

        # Different sites: {A, B} = 2 A B
        h2 = FockSpace(:c1) ⊗ FockSpace(:c2)
        @qnumbers a1::Destroy(h2, 1) a2::Destroy(h2, 2)
        @test isequal(anticommutator(a1, a2), 2 * a1 * a2)

        # Scalars
        @test anticommutator(2, 3) == 12
        @test isequal(anticommutator(2, a), 4 * a)
        @test isequal(anticommutator(a, 3), 6 * a)
    end
end

@testset "Scenario: Superradiant laser commutators" begin
    # Diagonal split during *(QAdd, QAdd) substitutes sum_idx → ext_idx
    # BEFORE site-sort so same-site composition (σ_αβ · σ_βγ = σ_αγ) fires
    # on the correct physical order. Reference equations (6)–(7):
    # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/
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

    # Eq. (6): i [H, a'σⱼ¹²]
    op3 = a' * σ(1, 2, j)
    result3 = expand_completeness(1im * commutator(H, op3))
    @test any(result3) do (term, c)
        term.ops == [a', σ(1, 2, j)] && isequal(c, _to_cnum(-1im * Δ))
    end
    @test any(result3) do (term, c)
        term.ops == [σ(2, 2, j)] && isequal(c, _to_cnum(1im * g(j)))
    end
    @test any(result3) do (term, c)
        term.ops == [a', a, σ(2, 2, j)] && isequal(c, _to_cnum(2im * g(j)))
    end
    @test any(result3) do (term, c)
        term.ops == [a', a] && isequal(c, _to_cnum(-1im * g(j)))
    end
    @test any(result3) do (term, c)
        term.ops == [σ(2, 1, i), σ(1, 2, j)] && isequal(c, _to_cnum(1im * g(i)))
    end
    @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(result3))
    # σⱼ¹¹ should NOT appear under canonical form
    @test !any(result3) do (term, _c)
        any(op -> op isa Transition && op.i == 1 && op.j == 1, term.ops)
    end

    # Eq. (7): i [H, σⱼ¹²σₖ²¹]
    op4 = σ(1, 2, j) * σ(2, 1, k)
    result4 = simplify(
        assume_distinct_index(
            expand_completeness(1im * commutator(H, op4)), [(j, k)],
        ),
    )
    @test any(result4) do (term, c)
        term.ops == [a, σ(2, 2, j), σ(2, 1, k)] && isequal(c, _to_cnum(2im * g(j)))
    end
    @test any(result4) do (term, c)
        term.ops == [a, σ(2, 1, k)] && isequal(c, _to_cnum(-1im * g(j)))
    end
    @test any(result4) do (term, c)
        term.ops == [a', σ(1, 2, j)] && isequal(c, _to_cnum(1im * g(k)))
    end
    @test any(result4) do (term, c)
        term.ops == [a', σ(1, 2, j), σ(2, 2, k)] && isequal(c, _to_cnum(-2im * g(k)))
    end
    @test !any(result4) do (term, _c)
        any(op -> op isa Transition && op.i == 1 && op.j == 1, term.ops)
    end
end

@testset "Scenario: Decay-channel triple product" begin
    # 2 * Σᵢ σ²¹_i · σ²²_j · σ¹²_i must be Σᵢ≠ⱼ 2 σ²²_i σ²²_j.
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

    @variables N
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    op = σ(2, 2, j)
    result = expand_completeness(2 * Σ(σ(2, 1, i) * op * σ(1, 2, i), i))

    @test length(result) == 1
    @test haskey(result, [σ(2, 2, i), σ(2, 2, j)])
    @test isequal(result[[σ(2, 2, i), σ(2, 2, j)]], _to_cnum(2))
    @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(result))
    @test !haskey(result, [σ(2, 2, j)])
end

@testset "Scenario: Commutator with double-sum Hamiltonian" begin
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

    result = 1im * commutator(H, σ(1, 2, k))
    @test any(result) do (term, _c)
        any(op -> op isa Transition && op.i == 1 && op.j == 2 && op.index == j, term.ops)
    end

    result3 = 1im * commutator(H, σ(2, 2, k))
    @test any(result3) do (term, _c)
        any(op -> op isa Transition && op.index == j, term.ops)
    end
end

@testset "Scenario: Lindblad dissipator (atomic decay)" begin
    # D[J](o) = J† o J - ½ {J†J, o} with J = σ_i¹².
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
        ), i,
    )

    @test iszero(simplify(decay(a' * a)))
    @test iszero(simplify(decay(σi(2, 2, j)) + γ * σi(2, 2, j)))
    @test iszero(simplify(decay(a' * σi(1, 2, j)) + (γ / 2) * a' * σi(1, 2, j)))
    @test iszero(
        simplify(
            assume_distinct_index(
                decay(σi(1, 2, j) * σi(2, 1, k)) + γ * σi(1, 2, j) * σi(2, 1, k),
                [(j, k)],
            ),
        ),
    )

    # Direct algebraic identity: σ_i²¹ σ_j²² σ_i¹² with i ≠ j must equal σ_i²² σ_j²².
    triple_sandwich = Σ(σi(2, 1, i) * σi(2, 2, j) * σi(1, 2, i), i)
    @test haskey(triple_sandwich, [σi(2, 2, i), σi(2, 2, j)])
    @test any(p -> p == (i, j) || p == (j, i), constraint_pairs(triple_sandwich))
    @test !haskey(triple_sandwich, [σi(2, 2, j)])
end

@testset "Regression: dissipator (PR #99)" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    a = Destroy(h, :a, 1)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

    @variables N Δ κ γ R ν
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    Jd = adjoint.(J)
    rates = [κ, γ, R, ν]

    decay(op) = simplify(
        expand_completeness(
            0.5 * rates[1] * (2 * Jd[1] * op * J[1] - Jd[1] * J[1] * op - op * Jd[1] * J[1]) +
                sum(
                Σ(
                        0.5 * rates[x] *
                        (2 * Jd[x] * op * J[x] - Jd[x] * J[x] * op - op * Jd[x] * J[x]),
                        i,
                    )
                    for x in 2:length(J)
            ),
        ),
    )

    @test isequal(decay(a' * a), simplify(-κ * (a' * a)))
    @test isequal(decay(σ(2, 2, j)), simplify(R - (R + γ) * σ(2, 2, j)))
end

@testset "Regression: rational coefficients cancel exactly" begin
    # Symbolics leaves `J0/D + (-J0)/D` un-combined, so a commuting `[Hₖ,op]=0` with a
    # rational coupling J0/|xᵢ-xⱼ|³ must still cancel exactly, not survive and bloat RHSs.
    h = NLevelSpace(:a1, (:g, :e)) ⊗ NLevelSpace(:a2, (:g, :e)) ⊗ NLevelSpace(:a3, (:g, :e))
    σ(x, y, k) = Transition(h, Symbol(:σ_, k), x, y, k)
    @variables J0 x1 x2 x3
    D(p, q) = abs(p - q)^3
    Hk = (J0 / D(x2, x3)) * (σ(:e, :g, 2) * σ(:g, :e, 3) + σ(:g, :e, 2) * σ(:e, :g, 3))
    @test iszero(commutator(Hk, σ(:g, :e, 1)))                       # disjoint atoms
    @test iszero((J0 / D(x1, x2)) * σ(:g, :e, 1) - (J0 / D(x1, x2)) * σ(:g, :e, 1))
    # A genuine (shared-space) rational-coefficient commutator is preserved:
    hf = FockSpace(:c)
    a = Destroy(hf, :a)
    @variables Δ
    @test iszero(commutator((J0 / Δ) * (a' * a), a) + (J0 / Δ) * a)  # [J/Δ·a†a, a] = -J/Δ·a
end
