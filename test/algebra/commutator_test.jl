using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: QAdd

@testset "Commutators and anticommutators" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = adjoint(a)

    @testset "bilinearity and scalar boundaries" begin
        @test commutator(1, 2) isa QAdd
        @test commutator(a, ad) isa QAdd
        @test iszero(commutator(1, 2))
        @test iszero(commutator(3, a))
        @test iszero(commutator(a, 5))
        @test iszero(commutator(2.0, ad))
        @test iszero(commutator(a, a))
        @test iszero(commutator(a, ad) - 1)
        @test iszero(commutator(ad, a) + 1)
        @test iszero(commutator(ad * a, a) + a)
        @test iszero(commutator(ad * a, ad) - ad)
        @test iszero(commutator(a + ad, a) + 1)
        @test iszero(commutator(a + ad, a + ad))

        left = a + 2 * ad
        right = a - ad
        @test iszero(
            commutator(left, right) -
                (
                commutator(a, a) - commutator(a, ad) +
                    2 * commutator(ad, a) - 2 * commutator(ad, ad)
            ),
        )
        @test iszero(commutator(a, commutator(a, ad)))
    end

    @testset "products on disjoint spaces remain independent" begin
        h2 = FockSpace(:left) ⊗ FockSpace(:right)
        left = Destroy(h2, :left, 1)
        right = Destroy(h2, :right, 2)
        @test iszero(commutator(left' * left, right))
        @test iszero(commutator(right, left' * left))
        @test iszero(commutator(left' * left, right' * right))
    end

    @testset "independent sites" begin
        h2 = FockSpace(:left) ⊗ FockSpace(:right)
        a1 = Destroy(h2, :left, 1)
        a2 = Destroy(h2, :right, 2)
        @test iszero(commutator(a1, a2))
        @test iszero(commutator(a1, a2'))
        @test iszero(commutator(a1' * a1, a2' * a2))

        # Different named modes in one Fock space are independent sites too.
        hs = FockSpace(:same_space)
        mode_a = Destroy(hs, :a)
        mode_b = Destroy(hs, :b)
        @test iszero(commutator(mode_a, mode_b'))
        @test iszero(commutator(mode_a, mode_b))
        @test iszero(commutator(mode_a', mode_b'))
    end

    @testset "spin, Pauli, and phase-space relations" begin
        hs = SpinSpace(:spin)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test iszero(commutator(Sx, Sy) - im * Sz)
        @test iszero(commutator(Sy, Sz) - im * Sx)
        @test iszero(commutator(Sz, Sx) - im * Sy)
        @test iszero(commutator(Sy, Sx) + im * Sz)

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test iszero(commutator(σx, σy) - 2im * σz)
        @test iszero(commutator(σy, σz) - 2im * σx)
        @test iszero(commutator(σz, σx) - 2im * σy)

        hq = PhaseSpace(:phase)
        x = Position(hq, :x)
        p = Momentum(hq, :p)
        @test iszero(commutator(x, p) - im)
        @test iszero(commutator(p, x) + im)
    end

    @testset "anticommutator" begin
        @test iszero(anticommutator(a, ad) - (1 + 2 * ad * a))
        @test iszero(anticommutator(a, a) - 2 * a * a)
        @test anticommutator(2, 3) == 12
        @test iszero(anticommutator(2, a) - 4 * a)
        @test iszero(anticommutator(a, 3) - 6 * a)

        hp = PauliSpace(:pauli)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        @test iszero(anticommutator(σx, σy))
        @test iszero(anticommutator(σx, σx) - 2)

        h2 = FockSpace(:one) ⊗ FockSpace(:two)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        @test iszero(anticommutator(a1, a2) - 2 * a1 * a2)
    end

    @testset "indexed sums" begin
        i = Index(h, :i, 10, h)
        j = Index(h, :j, 10, h)
        ai = IndexedOperator(a, i)
        ajd = IndexedOperator(ad, j)

        result = commutator(Σ(ai, i), Σ(ajd, j))
        @test iszero(result - 10)
        @test isempty(get_indices(result))

        @test iszero(commutator(IndexedOperator(a, i), IndexedOperator(a, j)'))
        @test iszero(
            commutator(IndexedOperator(a, i), IndexedOperator(a, i)') - 1
        )
    end

    @testset "public entry points remain inferable" begin
        @inferred commutator(a, ad)
        @inferred commutator(ad * a, a)
        @inferred commutator(a + ad, a + ad)
        @inferred anticommutator(a, ad)
    end
end

@testset "Scenario: indexed superradiant equations" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

    @variables N Δ
    g(idx) = IndexedVariable(:g, idx)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    expected3 = -1im * Δ * a' * σ(1, 2, j) +
        1im * g(j) * σ(2, 2, j) +
        2im * g(j) * a' * a * σ(2, 2, j) -
        1im * g(j) * a' * a +
        Σ(1im * g(i) * σ(2, 1, i) * σ(1, 2, j), i, [j])
    actual3 = expand_completeness(1im * commutator(H, a' * σ(1, 2, j)))
    @test iszero(simplify(actual3 - expected3))

    expected4 = 2im * g(j) * a * σ(2, 2, j) * σ(2, 1, k) -
        1im * g(j) * a * σ(2, 1, k) +
        1im * g(k) * a' * σ(1, 2, j) -
        2im * g(k) * a' * σ(1, 2, j) * σ(2, 2, k)
    actual4 = simplify(
        assume_distinct_index(
            expand_completeness(1im * commutator(H, σ(1, 2, j) * σ(2, 1, k))),
            [(j, k)],
        ),
    )
    @test iszero(
        simplify(
            actual4 - assume_distinct_index(expected4, [(j, k)]),
        )
    )
end

@testset "Scenario: indexed dissipator" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    @variables N γ
    i = Index(h, :i, N, 2)
    j = Index(h, :j, N, 2)
    k = Index(h, :k, N, 2)

    decay(op) = Σ(
        (γ / 2) * (
            2 * σ(2, 1, i) * op * σ(1, 2, i) -
                σ(2, 1, i) * σ(1, 2, i) * op -
                op * σ(2, 1, i) * σ(1, 2, i)
        ), i,
    )

    @test iszero(simplify(decay(a' * a)))
    @test iszero(simplify(decay(σ(2, 2, j)) + γ * σ(2, 2, j)))
    @test iszero(simplify(decay(a' * σ(1, 2, j)) + (γ / 2) * a' * σ(1, 2, j)))
    @test iszero(
        simplify(
            assume_distinct_index(
                decay(σ(1, 2, j) * σ(2, 1, k)) +
                    γ * σ(1, 2, j) * σ(2, 1, k),
                [(j, k)],
            ),
        ),
    )

    expected = 2 * Σ(σ(2, 2, i) * σ(2, 2, j), i, [j])
    actual = expand_completeness(
        2 * Σ(σ(2, 1, i) * σ(2, 2, j) * σ(1, 2, i), i),
    )
    @test iszero(simplify(actual - expected))
end

@testset "Scenario: rational coefficients and site identity" begin
    h = NLevelSpace(:a1, (:g, :e)) ⊗ NLevelSpace(:a2, (:g, :e)) ⊗
        NLevelSpace(:a3, (:g, :e))
    σ(x, y, slot) = Transition(h, Symbol(:σ_, slot), x, y, slot)
    @variables J0 x1 x2 x3
    distance(p, q) = abs(p - q)^3
    H = (J0 / distance(x2, x3)) *
        (σ(:e, :g, 2) * σ(:g, :e, 3) + σ(:g, :e, 2) * σ(:e, :g, 3))
    @test iszero(commutator(H, σ(:g, :e, 1)))
    @test iszero(
        (J0 / distance(x1, x2)) * σ(:g, :e, 1) -
            (J0 / distance(x1, x2)) * σ(:g, :e, 1)
    )

    hf = FockSpace(:cavity)
    a = Destroy(hf, :a)
    @variables Δ
    @test iszero(commutator((J0 / Δ) * a' * a, a) + (J0 / Δ) * a)
end

@testset "Scenario: multi-channel dissipator regression" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
    a = Destroy(h, :a, 1)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
    @variables N κ γ R ν
    i = Index(h, :i, N, 2)
    j = Index(h, :j, N, 2)

    jumps = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    adjoints = adjoint.(jumps)
    rates = [κ, γ, R, ν]
    decay(op) = simplify(
        expand_completeness(
            0.5 * rates[1] *
                (
                2 * adjoints[1] * op * jumps[1] -
                    adjoints[1] * jumps[1] * op - op * adjoints[1] * jumps[1]
            ) +
                sum(
                Σ(
                    0.5 * rates[channel] *
                        (
                        2 * adjoints[channel] * op * jumps[channel] -
                            adjoints[channel] * jumps[channel] * op -
                            op * adjoints[channel] * jumps[channel]
                    ),
                    i,
                ) for channel in 2:length(jumps)
            ),
        ),
    )

    @test isequal(decay(a' * a), simplify(-κ * (a' * a)))
    @test isequal(
        decay(σ(2, 2, j)),
        simplify(R - (R + γ) * σ(2, 2, j)),
    )
end
