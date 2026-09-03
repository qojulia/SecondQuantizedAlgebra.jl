using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: constraint_pairs

@testset "Canonical model workflows" begin
    @testset "Tavis-Cummings equations" begin
        cavity = FockSpace(:cavity)
        atom = NLevelSpace(:atom, 2)
        h = cavity ⊗ atom
        a = Destroy(h, :a, 1)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
        @variables N Δ
        g(k) = IndexedVariable(:g, k)
        i = Index(h, :i, N, atom)
        j = Index(h, :j, N, atom)
        k = Index(h, :k, N, atom)

        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
        @test iszero(
            simplify(
                commutator(H, a) - (Δ * a - Σ(g(i) * σ(1, 2, i), i)),
            )
        )
        @test iszero(
            simplify(
                commutator(H, a') - (-Δ * a' + Σ(g(i) * σ(2, 1, i), i)),
            )
        )
        @test iszero(
            simplify(
                commutator(H, a' * a) -
                    Σ(g(i) * (a * σ(2, 1, i) - a' * σ(1, 2, i)), i),
            )
        )
        @test iszero(
            simplify(
                commutator(H, σ(2, 2, j)) -
                    (g(j) * a' * σ(1, 2, j) - g(j) * a * σ(2, 1, j)),
            )
        )
        @test iszero(
            simplify(
                commutator(H, σ(1, 2, j)) -
                    g(j) * a * (σ(2, 2, j) - σ(1, 1, j)),
            )
        )

        expected = g(j) * a * (σ(2, 2, j) - σ(1, 1, j)) * σ(2, 1, k) +
            g(k) * a' * σ(1, 2, j) * (σ(1, 1, k) - σ(2, 2, k))
        actual = assume_distinct_index(
            commutator(H, σ(1, 2, j) * σ(2, 1, k)), [(j, k)],
        )
        @test iszero(simplify(actual - assume_distinct_index(expected, [(j, k)])))
        @test iszero(simplify(H - adjoint(H)))
    end

    @testset "Dicke model equations" begin
        h = FockSpace(:cavity) ⊗ SpinSpace(:collective_spin)
        a = Destroy(h, :a, 1)
        Sx = Spin(h, :S, 1, 2)
        Sy = Spin(h, :S, 2, 2)
        Sz = Spin(h, :S, 3, 2)
        @variables ω₀ ωₐ λ
        H = ω₀ * a' * a + ωₐ * Sz + λ * (a + a') * Sx

        @test iszero(simplify(commutator(H, a) - (-ω₀ * a - λ * Sx)))
        @test iszero(simplify(commutator(H, Sx) - im * ωₐ * Sy))
        @test iszero(
            simplify(
                commutator(H, Sy) - (-im * ωₐ * Sx + im * λ * (a + a') * Sz),
            )
        )
        @test iszero(simplify(commutator(H, Sz) + im * λ * (a + a') * Sy))
    end

    @testset "sum and completeness identities" begin
        h = FockSpace(:modes)
        a = Destroy(h, :a)
        @variables N c
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)

        @test iszero(simplify(commutator(Σ(ai, i), Σ(aj', j)) - N))
        @test iszero(simplify(Σ(ai * ai', i) - (Σ(ai' * ai, i) + N)))
        @test iszero(
            simplify(
                Σ((ai + c)^2, i) -
                    (Σ(ai * ai, i) + 2c * Σ(ai, i) + N * c^2),
            )
        )

        hp = PauliSpace(:pauli_sites)
        ip = Index(hp, :ip, N, hp)
        jp = Index(hp, :jp, N, hp)
        σx = IndexedOperator(Pauli(hp, :σ, 1), ip)
        σy = IndexedOperator(Pauli(hp, :σ, 2), jp)
        σz = IndexedOperator(Pauli(hp, :σ, 3), ip)
        @test iszero(
            simplify(
                commutator(Σ(σx, ip), Σ(σy, jp)) - Σ(2im * σz, ip),
            )
        )

        hn2 = NLevelSpace(:two_level, 2)
        n2 = Index(hn2, :n2, N, hn2)
        σ11(k) = IndexedOperator(Transition(hn2, :σ, 1, 1), k)
        σ22(k) = IndexedOperator(Transition(hn2, :σ, 2, 2), k)
        @test iszero(
            simplify(
                expand_completeness(Σ(σ11(n2), n2)) -
                    (N - Σ(σ22(n2), n2)),
            )
        )
        @test isequal(expand_completeness(Σ(σ22(n2), n2)), Σ(σ22(n2), n2))
    end

    @testset "averages and dissipators" begin
        h = NLevelSpace(:atom, 2)
        @variables N γ
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)

        s = Σ(σ(1, 2, i), i)
        @test is_indexed_sum(average(s))
        @test iszero(undo_average(average(s)) - s)

        dissipator(J, Jd, x) = Jd * x * J - (1 // 2) * (Jd * J * x + x * Jd * J)
        @test iszero(
            simplify(
                dissipator(σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)) + σ(2, 2, j),
            )
        )
        @test iszero(
            simplify(
                dissipator(σ(1, 2, j), σ(2, 1, j), σ(1, 1, j)) - σ(2, 2, j),
            )
        )

        collective = dissipator(
            Σ(σ(1, 2, i), i), Σ(σ(2, 1, j), j), σ(2, 2, j),
        )
        @test !iszero(collective)
        @test j in get_indices(collective)

        # The collective dissipator is linear in its observable argument. This
        # checks the distributed sum/product result through the public algebra,
        # rather than accepting a nonzero expression as a smoke test.
        observable_sum = dissipator(
            Σ(σ(1, 2, i), i), Σ(σ(2, 1, j), j),
            σ(2, 2, j) + σ(1, 1, j),
        )
        @test iszero(
            simplify(
                observable_sum -
                    (
                    collective + dissipator(
                        Σ(σ(1, 2, i), i), Σ(σ(2, 1, j), j), σ(1, 1, j),
                    )
                ),
            ),
        )
    end

    @testset "structural regression: two-atom correlator" begin
        cavity = FockSpace(:cavity)
        atom = NLevelSpace(:atom, 2)
        h = cavity ⊗ atom
        a = Destroy(h, :a, 1)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
        @variables N Δ
        g(k) = IndexedVariable(:g, k)
        i = Index(h, :i, N, atom)
        j = Index(h, :j, N, atom)
        k = Index(h, :k, N, atom)
        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
        raw = commutator(H, σ(1, 2, j) * σ(2, 1, k))
        # diagonal partition must reduce to four 3-operator terms
        @test length(raw) == 4
        for term in keys(raw.arguments)
            @test length(term.ops) == 3
        end
        # coefficient-only drive regression: no diagonal leak
        l = Index(h, :l, N, atom)
        u(α, β) = DoubleIndexedVariable(:u, α, β)
        Hdrive = Σ(Σ(u(j, k) * (σ(2, 1, j) + σ(1, 2, j)), j), k)
        comm = commutator(σ(2, 2, l), Hdrive)
        expected = u(l, l) * σ(2, 1, l) - u(l, l) * σ(1, 2, l) +
            Σ(u(l, k) * σ(2, 1, l) - u(l, k) * σ(1, 2, l), k, Index[l])
        @test iszero(simplify(commutator(σ(2, 2, l), Hdrive) - expected))
        @test any(pair -> pair == (k, l) || pair == (l, k), constraint_pairs(comm))
    end

    @testset "hand-derived exact forms are reproducible" begin
        # Tavis-Cummings exact commutators with public iszero(simplify) as independent reference
        cavity = FockSpace(:cavity)
        atom = NLevelSpace(:atom2, 2)
        h = cavity ⊗ atom
        a = Destroy(h, :a, 1)
        σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
        @variables N Δ2
        g(k) = IndexedVariable(:g, k)
        i = Index(h, :i, N, atom)
        j = Index(h, :j, N, atom)
        H = -Δ2 * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
        @test iszero(simplify(commutator(H, a) - (Δ2 * a - Σ(g(i) * σ(1, 2, i), i))))
        @test iszero(simplify(commutator(H, σ(2, 2, j)) - (g(j) * a' * σ(1, 2, j) - g(j) * a * σ(2, 1, j))))
        # Dicke exact form without simplify tolerance
        hD = FockSpace(:cavD) ⊗ SpinSpace(:spinD)
        aD = Destroy(hD, :a, 1)
        Sx = Spin(hD, :S, 1, 2)
        Sy = Spin(hD, :S, 2, 2)
        Sz = Spin(hD, :S, 3, 2)
        @variables ω₀ ωₐ λ
        HD = ω₀ * aD' * aD + ωₐ * Sz + λ * (aD + aD') * Sx
        @test iszero(simplify(commutator(HD, Sx) - im * ωₐ * Sy))
        @test iszero(simplify(commutator(HD, Sy) - (-im * ωₐ * Sx + im * λ * (aD + aD') * Sz)))
    end
end
