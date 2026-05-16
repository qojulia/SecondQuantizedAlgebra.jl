using SecondQuantizedAlgebra
using QuantumOpticsBase
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, CNum, sorted_arguments

@testset "normal_order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Already normal-ordered" begin
        m = ad * a
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
        @test isequal(result, m)
    end

    @testset "Single commutation: a*a† → a†a + 1" begin
        m = a * ad
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "QSym passthrough" begin
        result = normal_order(a)
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "QAdd — normal-orders each term" begin
        expr = a * ad + ad * a
        result = normal_order(expr)
        @test result isa QAdd
    end

    @testset "Different spaces don't commute" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "Return type is always QAdd" begin
        @test normal_order(a) isa QAdd
        @test normal_order(2 * a) isa QAdd
        @test normal_order(a + ad) isa QAdd
    end
end

@testset "normal_order — Transition" begin
    h = NLevelSpace(:atom, 3, 1)

    @testset "Composition: |1⟩⟨2| · |2⟩⟨3| = |1⟩⟨3|" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ23 = Transition(h, :σ, 2, 3)
        result = σ12 * σ23  # eager composition
        @test result isa QAdd
        @test length(result) == 1
        ops = first(keys(result.arguments)).ops
        @test ops[1] isa Transition
        @test ops[1].i == 1 && ops[1].j == 3
    end

    @testset "Orthogonality: |1⟩⟨2| · |3⟩⟨1| = 0" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        result = σ12 * σ31
        @test iszero(result)
    end

    @testset "Orthogonality in longer product" begin
        σ12 = Transition(h, :σ, 1, 2)
        σ31 = Transition(h, :σ, 3, 1)
        σ23 = Transition(h, :σ, 2, 3)
        result = σ12 * σ31 * σ23
        @test iszero(result)
    end

    @testset "Ground state rewriting — single space" begin
        σ11 = Transition(h, :σ, 1, 1)
        result = expand_completeness(normal_order(σ11))
        @test result isa QAdd
        @test length(result) == 3
        simplified = simplify(result)
        @test length(simplified) == 3
    end

    @testset "Ground state rewriting — ProductSpace" begin
        hc = FockSpace(:c)
        ha = NLevelSpace(:atom, 2, 1)
        hp = hc ⊗ ha
        σ11 = Transition(hp, :σ, 1, 1, 2)
        result = expand_completeness(normal_order(σ11))
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "Ground state rewriting — recursive (double projection)" begin
        h2 = NLevelSpace(:atom, 2, 1)
        σ11 = Transition(h2, :σ, 1, 1)
        m = σ11 * σ11
        result = expand_completeness(normal_order(m))
        @test result isa QAdd
    end
end

@testset "normal_order — Pauli" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)

    @testset "σⱼ·σⱼ = 1 (idempotent)" begin
        for σ in (σx, σy, σz)
            result = σ * σ
            @test result isa QAdd
            @test length(result) == 1
            (term, c) = only(collect(result))
            @test isempty(term.ops)
            @test c == 1
        end
    end

    @testset "Pauli cycle σₐ·σᵦ = i·σ_c, σᵦ·σₐ = -i·σ_c" begin
        # (a, b, c) with a × b = c under Levi-Civita; a < b ⇒ +i, b > a ⇒ -i.
        cycle = (
            (σx, σy, σz, 3),
            (σy, σz, σx, 1),
            (σz, σx, σy, 2),
        )
        for (a, b, _c, axis) in cycle
            result = a * b
            @test result isa QAdd
            @test length(result) == 1
            (term, c) = only(collect(result))
            @test c == im
            @test length(term.ops) == 1
            @test term.ops[1] isa Pauli
            @test term.ops[1].axis == axis

            (term2, c2) = only(collect(b * a))
            @test c2 == -im
            @test term2.ops[1].axis == axis
        end
    end
end

@testset "normal_order — Spin" begin
    h = SpinSpace(:s)
    Sx = Spin(h, :S, 1)
    Sy = Spin(h, :S, 2)
    Sz = Spin(h, :S, 3)

    @testset "Sy·Sx gets reordered (swap + commutator)" begin
        result = Sy * Sx  # eager ordering
        @test result isa QAdd
        @test length(result) >= 2
    end

    @testset "Same axis — already ordered" begin
        result = Sx * Sx
        @test result isa QAdd
        @test length(result) == 1
    end

    @testset "Sx·Sy already in order" begin
        result = Sx * Sy
        @test result isa QAdd
        @test length(result) == 1
    end
end

@testset "normal_order — Type stability" begin
    hf = FockSpace(:c)
    a = Destroy(hf, :a)

    @inferred normal_order(a)
    @inferred normal_order(a * a')
    @inferred normal_order(a' * a)
    @inferred normal_order(a * a' + a' * a)

    hn = NLevelSpace(:atom, 3, 1)
    σ12 = Transition(hn, :σ, 1, 2)
    σ23 = Transition(hn, :σ, 2, 3)
    @inferred normal_order(σ12 * σ23)
    @inferred expand_completeness(normal_order(σ12))

    hp = PauliSpace(:p)
    σx = Pauli(hp, :σ, 1)
    σy = Pauli(hp, :σ, 2)
    @inferred normal_order(σx * σy)

    hs = SpinSpace(:s)
    Sy = Spin(hs, :S, 2)
    Sx = Spin(hs, :S, 1)
    @inferred normal_order(Sy * Sx)

    hps = PhaseSpace(:q)
    x = Position(hps, :x)
    p = Momentum(hps, :p)
    @inferred normal_order(p * x)
end

@testset "normal_order — Allocations" begin
    hf = FockSpace(:c)
    a = Destroy(hf, :a)

    normal_order(a * a')
    normal_order(a' * a)
    normal_order(a * a' + a' * a)

    @test @allocations(normal_order(a * a')) < 1000
    @test @allocations(normal_order(a' * a)) < 1000
    @test @allocations(normal_order(a * a' + a' * a)) < 2000

    hn = NLevelSpace(:atom, 3, 1)
    σ11 = Transition(hn, :σ, 1, 1)
    expand_completeness(normal_order(σ11))
    @test @allocations(expand_completeness(normal_order(σ11))) < 5000

    hps = PhaseSpace(:q)
    x = Position(hps, :x)
    p = Momentum(hps, :p)
    normal_order(p * x)
    @test @allocations(normal_order(p * x)) < 1000
end

@testset "Cross-cutting smoke" begin
    @testset "Commutation via normal_order" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        # a * a† normal-ordered → a†a + 1
        result = normal_order(a * a')
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "Distinct modes on same space don't commute" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        a = Destroy(h1, :a)
        b = Destroy(h2, :b)
        result = normal_order(a * b')
        @test length(result) == 1
    end

    @testset "Multi-space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

        @test length(normal_order(a * b')) == 1   # different subspaces
        @test length(normal_order(a * a')) == 2   # same-subspace commutator
    end

    @testset "Higher-order normal ordering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        result = normal_order(a * a * a' * a')
        @test result isa QAdd
        @test length(result) >= 3
    end
end

@testset "Scenario: Jaynes-Cummings model" begin
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

@testset "Scenario: Pauli algebra via normal_order" begin
    h = PauliSpace(:p)
    σx = Pauli(h, :σ, 1)
    σy = Pauli(h, :σ, 2)
    σz = Pauli(h, :σ, 3)

    # σx·σy·σz = iσz·σz = i·I
    result = normal_order(σx * σy * σz)
    @test result isa QAdd
    @test length(result) == 1
    @test isempty(operators(only(sorted_arguments(result))))
    @test prefactor(only(sorted_arguments(result))) == im
end

@testset "Scenario: Mixed Fock + Pauli" begin
    h = FockSpace(:c) ⊗ PauliSpace(:p)
    @qnumbers a::Destroy(h, 1)
    σz = Pauli(h, :σ, 3, 2)

    m = a' * a * σz
    @test m isa QAdd
    @test length(operators(only(sorted_arguments(m)))) == 3
    @test normal_order(a * a' * σz) isa QAdd
end

@testset "Scenario: Mollow triplet" begin
    # H = Δ σ_ee + Ω (σ_ge + σ_eg) — driven two-level atom
    # https://qojulia.github.io/QuantumCumulants.jl/stable/examples/mollow/
    h = NLevelSpace(:atom, (:g, :e))
    @variables Δ Ω
    σ(α, β) = Transition(h, :σ, α, β)

    H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))

    @test iszero(
        simplify(
            expand_completeness(commutator(H, σ(:g, :e))) -
                (-Δ * σ(:g, :e) + Ω * (2 * σ(:e, :e) - 1))
        )
    )
    @test iszero(
        simplify(
            expand_completeness(commutator(H, σ(:e, :g))) -
                (Δ * σ(:e, :g) - Ω * (2 * σ(:e, :e) - 1))
        )
    )
    @test iszero(
        simplify(
            commutator(H, σ(:e, :e)) -
                (Ω * σ(:g, :e) - Ω * σ(:e, :g))
        )
    )

    # average() is a no-op wrap: identity round-trip
    @test isequal(undo_average(average(H)), H)
end

@testset "Scenario: Single-atom laser (JC)" begin
    # H = Δ a†a + g (a† σ_ge + a σ_eg)
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha

    @variables Δ g
    a = Destroy(h, :a)
    s(α, β) = Transition(h, :σ, α, β)

    H = Δ * a' * a + g * (a' * s(:g, :e) + a * s(:e, :g))

    @test iszero(simplify(commutator(H, a) + Δ * a + g * s(:g, :e)))
    @test iszero(
        simplify(
            expand_completeness(commutator(H, s(:g, :e))) -
                (-g * a + 2g * a * s(:e, :e))
        )
    )
    @test iszero(simplify(commutator(a' * a, a' * a)))
    @test iszero(simplify(H - adjoint(H)))
end

@testset "Scenario: Mollow triplet — extended" begin
    h = NLevelSpace(:atom, (:g, :e))
    @variables Δ Ω
    σ(α, β) = Transition(h, :σ, α, β)

    H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))

    @test iszero(simplify(H - adjoint(H)))
    @test iszero(simplify(σ(:e, :e) * σ(:e, :e) - σ(:e, :e)))
    @test iszero(simplify(expand_completeness(σ(:g, :e) * σ(:e, :g)) - (1 - σ(:e, :e))))
    @test iszero(simplify(σ(:e, :g) * σ(:g, :e) - σ(:e, :e)))
    @test iszero(
        simplify(
            expand_completeness(commutator(σ(:g, :e), σ(:e, :g))) -
                (1 - 2 * σ(:e, :e))
        )
    )
end

@testset "Scenario: JC — extended" begin
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
    @test iszero(
        simplify(
            expand_completeness(commutator(H, s(:e, :g))) -
                (g * a' - 2g * a' * s(:e, :e))
        )
    )
    @test iszero(simplify(commutator(H, a') - (Δ * a' + g * s(:e, :g))))
    @test iszero(simplify(commutator(a, s(:g, :e))))
end
