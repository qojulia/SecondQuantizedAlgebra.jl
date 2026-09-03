using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra:
    QAdd,
    QField,
    constraint_pairs,
    expim,
    get_sum_body,
    get_sum_indices,
    get_sum_non_equal,
    has_sum_metadata,
    indexed_sum

@testset "Public algebraic properties" begin
    @testset "canonicalization and average round trips" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        a = Destroy(h, :a, 1)
        σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
        @variables N Δ g
        i = Index(h, :i, N, 2)
        j = Index(h, :j, N, 2)
        k = Index(h, :k, N, 2)

        @variables κS
        corpus = [
            ("a a'", a * a'),
            ("a' a + a a'", a' * a + a * a'),
            ("Σ_gσ", Σ(g * σ(1, 2, i), i)),
            ("Σ_gσ * σ", Σ(g * σ(1, 2, i), i) * σ(2, 1, j)),
            (
                "commutator H σ_j",
                commutator(-Δ * a' * a + Σ(g * a' * σ(1, 2, i), i), σ(1, 2, j)),
            ),
            (
                "commutator H σ_j σ_k",
                commutator(
                    -Δ * a' * a + Σ(g * a' * σ(1, 2, i), i),
                    σ(1, 2, j) * σ(2, 1, k),
                ),
            ),
            ("expand_completeness Σ", expand_completeness(Σ(σ(1, 1, i), i))),
            ("substitute Σ", substitute(Σ(g * σ(1, 2, i), i), Dict(g => 2))),
            ("change_index Σ", change_index(Σ(g * σ(1, 2, i), i), i, j)),
            ("Σ minus itself", Σ(g * σ(1, 2, i), i) - Σ(g * σ(1, 2, i), i)),
            ("Σ product diag split", Σ(g * σ(1, 2, i), i) * σ(2, 1, j)),
            (
                "simplify H-H",
                simplify(
                    (-Δ * a' * a + Σ(g * a' * σ(1, 2, i), i)) -
                        (-Δ * a' * a + Σ(g * a' * σ(1, 2, i), i)),
                ),
            ),
            ("adjoint Σ-product", (Σ(g * σ(1, 2, i), i) * σ(2, 1, j))'),
            ("normal_order sum prod", normal_order(Σ(g * σ(2, 2, i), i) * a' * a)),
            (
                "partial cancel i dies j lives",
                (Σ(g * σ(1, 2, i), i) + Σ(g * σ(2, 1, j), j)) - Σ(g * σ(1, 2, i), i),
            ),
            (
                "double pin i twice",
                Σ(g * σ(1, 2, i) * σ(2, 1, i), i) * σ(2, 1, j) * σ(1, 2, k),
            ),
            ("conjugate Displace", conjugate(a' * a, Displace(a, Δ))),
            ("conjugate Squeeze", conjugate(a' * a, Squeeze(a, κS))),
            (
                "sum-sum commutator Bug A",
                let hf = FockSpace(:fb), bi = IndexedOperator(Destroy(hf, :b), Index(hf, :iF, N, hf)), jF = Index(hf, :jF, N, hf), bj_dag = IndexedOperator(Create(hf, :b), jF)
                    commutator(Σ(bi, Index(hf, :iF, N, hf)), Σ(bj_dag, jF))
                end,
            ),
            (
                "anticommutator on sums",
                let hf = FockSpace(:fb2), iF = Index(hf, :iF, N, hf), jF = Index(hf, :jF, N, hf), bi = IndexedOperator(Destroy(hf, :b), iF), bj_dag = IndexedOperator(Create(hf, :b), jF)
                    anticommutator(Σ(bi, iF), Σ(bj_dag, jF))
                end,
            ),
            (
                "Pauli sum-sum",
                let hP = PauliSpace(:p_inv), iP = Index(hP, :ip, N, hP), jP = Index(hP, :jp, N, hP), σxi = IndexedOperator(Pauli(hP, :σ, 1), iP), σyj = IndexedOperator(Pauli(hP, :σ, 2), jP)
                    commutator(Σ(σxi, iP), Σ(σyj, jP))
                end,
            ),
            (
                "Dicke [H, σy_j]",
                let hcD = FockSpace(:cD), hpD = PauliSpace(:aD), hD = hcD ⊗ hpD, aD = Destroy(hD, :a, 1), iD = Index(hD, :id, N, hpD), jD = Index(hD, :jd, N, hpD), σxD(kk) = IndexedOperator(Pauli(hD, :σ, 1, 2), kk), σyD(kk) = IndexedOperator(Pauli(hD, :σ, 2, 2), kk), σzD(kk) = IndexedOperator(Pauli(hD, :σ, 3, 2), kk), ω₀ = 1.0, ωₐ = 1.0, λ = 0.5
                    H_dicke = ω₀ * aD' * aD + ωₐ * Σ(σzD(iD), iD) + λ * (aD + aD') * Σ(σxD(iD), iD)
                    commutator(H_dicke, σyD(jD))
                end,
            ),
            (
                "coefficient-index leak u(j,i)",
                let ha = NLevelSpace(:atom_leak, 2), h = FockSpace(:cav_leak) ⊗ ha, iL = Index(h, :iL, N, ha), jL = Index(h, :jL, N, ha), kL = Index(h, :kL, N, ha), u(k, l) = DoubleIndexedVariable(:u, k, l), σL(kk) = IndexedOperator(Transition(h, :σ, 1, 2, 2), kk)
                    inner = Σ(u(jL, kL) * (σL(jL) + σL(kL)), jL)
                    outer = Σ(inner, kL)
                    commutator(σL(iL)' * σL(iL), outer)
                end,
            ),
            ("expand_completeness on Σ", expand_completeness(Σ(σ(1, 1, i) * a, i))),
            (
                "substitute on Σ",
                substitute(Σ(g * σ(1, 2, i), i), Dict(σ(1, 2, j) => 2 * σ(1, 2, j))),
            ),
        ]

        for (label, expression) in corpus
            @test iszero(simplify(expression - normal_order(expression))) ||
                error("canonical invariant failed: $label")
            @test iszero(simplify(undo_average(average(expression)) - expression)) ||
                error("average roundtrip failed: $label")
        end
        # Rotation unitary
        @test iszero(simplify(conjugate(a, Rotation(a, Δ)) - expim(-Δ) * a))
        hF = FockSpace(:f)
        @qnumbers b::Destroy(hF)
        @variables κd tt Δt ttt
        @test iszero(
            simplify(
                conjugate(conjugate(b' * b, Rotation(b, Δ) * Squeeze(b, κd)), inv(Rotation(b, Δ) * Squeeze(b, κd))) - b' * b,
            ),
        )
        gt = gauge_term(Squeeze(b, κd * ttt, Δt * ttt, ttt))
        @test gt isa QAdd || gt isa QField
    end

    @testset "cancellation never leaves observable scope" begin
        h = FockSpace(:sites)
        a = Destroy(h, :a)
        @variables N
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        left = Σ(IndexedOperator(a, i), i)
        right = Σ(IndexedOperator(a, j), j)

        @test isempty(get_indices(left - left))
        @test Set(get_indices(left - right)) == Set([i, j])
        @test get_indices((left + right) - left) == [j]
        @test isempty(get_indices(commutator(left, Σ(IndexedOperator(a', j), j))))
    end

    @testset "sum-sum commutator and indexed cancellations" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        a = Destroy(h, :a, 1)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
        @variables N
        i = Index(h, :i, N, 2)
        j = Index(h, :j, N, 2)
        hf = FockSpace(:fmodes)
        b = Destroy(hf, :b)
        iF = Index(hf, :iF, N, hf)
        jF = Index(hf, :jF, N, hf)
        # sum-sum commutator yields scalar, no indices
        c = commutator(Σ(IndexedOperator(b, iF), iF), Σ(IndexedOperator(b', jF), jF))
        @test iszero(simplify(c - N))
        @test isempty(get_indices(c))
        # constants inside Σ pick up range factor
        @test iszero(simplify(Σ(IndexedOperator(b, iF) + 1, iF) - (Σ(IndexedOperator(b, iF), iF) + N)))
        # Pauli sum-sum
        hP = PauliSpace(:pauli_inv)
        iP = Index(hP, :ip, N, hP)
        jP = Index(hP, :jp, N, hP)
        σxi = IndexedOperator(Pauli(hP, :σ, 1), iP)
        σyj = IndexedOperator(Pauli(hP, :σ, 2), jP)
        p = Σ(σxi, iP) * Σ(σyj, jP)
        @test length(p) == 2
        @test any(pair -> pair == (iP, jP) || pair == (jP, iP), constraint_pairs(p))
    end

    @testset "public algebraic laws" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        b = Destroy(h, :b)
        ha = NLevelSpace(:atom, 2)
        hmix = h ⊗ ha
        af = Destroy(hmix, :a, 1)
        σ(α, β, k) = IndexedOperator(Transition(hmix, :σ, α, β, 2), k)
        @variables N gval
        i = Index(hmix, :i, N, ha)
        j = Index(hmix, :j, N, ha)
        base = [
            ("a", a),
            ("a'", a'),
            ("a'*a", a' * a),
            ("a+a'", a + a'),
            ("σ_j^{12}", σ(1, 2, j)),
            ("σ_j^{11}", σ(1, 1, j)),
            ("a*σ_j", af * σ(1, 2, j)),
            ("Σ g σ", Σ(gval * σ(1, 2, i), i)),
        ]
        for (la, A) in base, (lb, B) in base
            (A isa QAdd && B isa QAdd && !isempty(get_indices(A)) && !isempty(get_indices(B)) && any(idx -> idx in get_indices(B), get_indices(A))) && continue
            @test iszero(simplify(commutator(A, B) + commutator(B, A))) ||
                error("antisymmetry $la,$lb")
        end
        for (A, B, C) in ((a, a', b), (a + a', a, b'), (a' * a, a + b, a' - b'))
            @test iszero(simplify(commutator(A, B + C) - commutator(A, B) - commutator(A, C)))
        end
        # average linearity over base
        for (la, A) in base, (lb, B) in base
            (A isa QAdd && B isa QAdd && !isempty(get_indices(A)) && !isempty(get_indices(B)) && any(idx -> idx in get_indices(B), get_indices(A))) && continue
            d = simplify(undo_average(average(A + B)) - undo_average(average(A) + average(B)))
            @test iszero(d) || error("average linearity $la $lb")
        end
        for (label, A) in base
            qa = A isa QAdd ? A : (1 * A)
            @test iszero(simplify(undo_average(average(qa)) - qa)) || error("roundtrip $label")
            @test iszero(commutator(A, A))
            @test isequal(adjoint(adjoint(qa)), qa)
        end

        hs = SpinSpace(:spin)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        jacobi = commutator(commutator(Sx, Sy), Sz) +
            commutator(commutator(Sy, Sz), Sx) +
            commutator(commutator(Sz, Sx), Sy)
        @test iszero(simplify(jacobi))
    end

    @testset "indexed average metadata describes its scope" begin
        h = NLevelSpace(:atom, 2)
        @variables N
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        σ(k) = IndexedOperator(Transition(h, :σ, 1, 2), k)

        averaged = average(Σ(σ(i), i, [j]))
        @test is_indexed_sum(averaged)
        @test has_sum_metadata(averaged)
        @test get_sum_indices(averaged) == [i]
        @test get_sum_non_equal(averaged) == [(i, j)]
        @test iszero(undo_average(averaged) - Σ(σ(i), i, [j]))

        plain = average(Σ(σ(i), i))
        @test get_sum_indices(plain) == [i]
        @test isempty(get_sum_non_equal(plain))

        # Sum-independent term must not carry metadata
        hf = FockSpace(:f)
        @qnumbers bf::Destroy(hf)
        iF = Index(hf, :iF, N, hf)
        jF = Index(hf, :jF, N, hf)
        c = commutator(Σ(IndexedOperator(bf, iF), iF), Σ(IndexedOperator(bf', jF), jF))
        @test !has_sum_metadata(average(c))

        # Split metadata: dep vs indep
        hmix = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        a2 = Destroy(hmix, :a, 1)
        σ2(k) = IndexedOperator(Transition(hmix, :σ, 1, 2, 2), k)
        g2(k) = IndexedVariable(:g, k)
        i2 = Index(hmix, :i2, N, 2)
        j2 = Index(hmix, :j2, N, 2)
        prod = Σ(g2(i2) * σ2(i2), i2) * σ2(j2)
        avg = average(prod)
        @test !isequal(avg, 0)
        # off-diagonal constraint should be visible via constraint_pairs on undo
        restored = undo_average(avg)
        @test !iszero(restored)
        @test iszero(simplify(undo_average(average(Σ(g2(i2) * σ2(i2), i2))) - Σ(g2(i2) * σ2(i2), i2)))
    end
end
