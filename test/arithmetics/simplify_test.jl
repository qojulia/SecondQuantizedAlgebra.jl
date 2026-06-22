using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QField, sorted_arguments, _single_qadd, _to_cnum

@testset "simplify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Collect like terms" begin
        s = _single_qadd(_to_cnum(5), Op[ad, a])
        result = simplify(s)
        @test result isa QAdd
        @test length(result) == 1
        @test prefactor(only(sorted_arguments(result))) == 5
    end

    @testset "Remove zero terms" begin
        s = _single_qadd(_to_cnum(3), Op[ad, a])
        result = simplify(s)
        @test length(result) == 1
        @test prefactor(only(sorted_arguments(result))) == 3
    end

    @testset "a + a = 2a" begin
        s = a + a
        result = simplify(s)
        @test length(result) == 1
        @test prefactor(only(sorted_arguments(result))) == 2
    end

    @testset "Symbolic prefactors" begin
        @variables g h_sym
        s = _single_qadd(_to_cnum(g + h_sym), Op[ad, a])
        result = simplify(s)
        @test length(result) == 1
    end

    @testset "Fock: simplify preserves eager normal ordering" begin
        m = a * ad  # eager `*` already produces a†a + 1
        result = simplify(m)
        @test result isa QAdd
        @test length(result) == 2  # a†a + 1 (from eager ordering)
    end

    @testset "Transition: simplify applies composition" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        result = simplify(σ12 * σ23)
        @test length(result) == 1
        @test operators(only(sorted_arguments(result)))[1].l1 == 1
        @test operators(only(sorted_arguments(result)))[1].l2 == 3
    end

    @testset "Transition: simplify applies orthogonality" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ31 = Transition(hn, :σ, 3, 1)
        result = simplify(σ12 * σ31)
        @test iszero(result)
    end

    @testset "Pauli: simplify applies σⱼ²=I and product rule" begin
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)

        # σx² = I
        result = simplify(σx * σx)
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))
        @test prefactor(only(sorted_arguments(result))) == 1

        # σx·σy = iσz
        result = simplify(σx * σy)
        @test prefactor(only(sorted_arguments(result))) == im
        @test operators(only(sorted_arguments(result)))[1].l1 == 3
    end

    @testset "simplify on eagerly ordered expression" begin
        result = simplify(a * ad)
        @test result isa QAdd
        @test length(result) == 2  # eager `*` already applied
    end

    @testset "SymbolicUtils.simplify on QField" begin
        @variables g
        s = _single_qadd(_to_cnum(2g), Op[ad, a])
        result = SymbolicUtils.simplify(s)
        @test result isa QAdd
    end

    @testset "Symbolics.expand on QField" begin
        expr = (a + ad) * (a + ad)
        result = Symbolics.expand(expr)
        @test result isa QAdd
        @test length(result) == 4

        # Prefactor-level distribution: polynomial coefficients are stored in
        # canonical expanded form, so a symbolic prefactor like (g+h)^2 is
        # already g^2 + 2gh + h^2 at construction (the package-wide 'always
        # expand' invariant), and `Symbolics.expand` on it is a no-op.
        @variables g h
        coeff_expr = ((g + h)^2) * a
        @test length(coeff_expr) == 1
        @test isequal(coeff_expr, (g^2 + 2 * g * h + h^2) * a)

        expanded = Symbolics.expand(coeff_expr)
        @test length(expanded) == 1
        @test isequal(expanded, (g^2 + 2 * g * h + h^2) * a)
    end

    @testset "Idempotency: simplify(simplify(x)) == simplify(x)" begin
        # Fock — basic
        @test isequal(simplify(simplify(a * ad)), simplify(a * ad))
        @test isequal(simplify(simplify(a * ad * a * ad)), simplify(a * ad * a * ad))
        @test isequal(simplify(simplify(a * ad + ad * a)), simplify(a * ad + ad * a))

        # Fock — higher power
        @test isequal(simplify(simplify((a * ad)^4)), simplify((a * ad)^4))

        # Transition — composition
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @test isequal(simplify(simplify(σ12 * σ23)), simplify(σ12 * σ23))

        # Transition — ground state rewriting
        σ11 = Transition(hn, :σ, 1, 1)
        σ21 = Transition(hn, :σ, 2, 1)
        @test isequal(
            expand_completeness(simplify(expand_completeness(simplify(σ11 * σ12 * σ21)))),
            expand_completeness(simplify(σ11 * σ12 * σ21))
        )

        # Pauli — long chain
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)
        @test isequal(simplify(simplify(σx * σy * σz)), simplify(σx * σy * σz))
        @test isequal(
            simplify(simplify(σx * σy * σz * σx * σy)),
            simplify(σx * σy * σz * σx * σy)
        )

        # Pauli — sum-of-products
        pauli_expr = (σx + σy) * (σy + σz)
        @test isequal(simplify(simplify(pauli_expr)), simplify(pauli_expr))

        # Spin — Casimir S²
        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @test isequal(simplify(simplify(Sy * Sx)), simplify(Sy * Sx))
        S2 = Sx * Sx + Sy * Sy + Sz * Sz
        @test isequal(simplify(simplify(S2)), simplify(S2))

        # PhaseSpace
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @test isequal(simplify(simplify(p * x)), simplify(p * x))
        @test isequal(simplify(simplify(x * p * x * p)), simplify(x * p * x * p))

        # Multi-mode Fock
        h3 = FockSpace(:a) ⊗ FockSpace(:b) ⊗ FockSpace(:c)
        a1 = Destroy(h3, :a, 1)
        a2 = Destroy(h3, :b, 2)
        a3 = Destroy(h3, :c, 3)
        H_bs = a1' * a2 + a2' * a1 + a2' * a3 + a3' * a2
        @test isequal(simplify(simplify(H_bs)), simplify(H_bs))
        cross = a1 * a2' * a1' * a2 * a3 * a3'
        @test isequal(simplify(simplify(cross)), simplify(cross))

        # Mixed Hilbert space: Fock + NLevel
        hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a_jc = Destroy(hjc, :a, 1)
        σ_jc = Transition(hjc, :σ, 1, 2, 2)
        H_jc = a_jc' * σ_jc + a_jc * σ_jc'
        @test isequal(simplify(simplify(H_jc)), simplify(H_jc))
        c_jc = commutator(H_jc, a_jc' * a_jc)
        @test isequal(simplify(simplify(c_jc)), simplify(c_jc))

        # Indexed operators
        hf = FockSpace(:f)
        a_idx = Destroy(hf, :a)
        i_idx = Index(hf, :i, 10, hf)
        j_idx = Index(hf, :j, 10, hf)
        ai = IndexedOperator(a_idx, i_idx)
        aj = IndexedOperator(a_idx, j_idx)
        N_sum = Σ(ai' * ai, i_idx)
        @test isequal(simplify(simplify(N_sum)), simplify(N_sum))
        c_idx = commutator(N_sum, aj)
        @test isequal(simplify(simplify(c_idx)), simplify(c_idx))
    end

    @testset "Return type is QAdd" begin
        @test simplify(a + ad) isa QAdd
        @test simplify(a) isa QAdd
        @test simplify(2 * a) isa QAdd
    end

    @testset "Type stability" begin
        @inferred simplify(a)
        @inferred simplify(a * ad)
        @inferred simplify(a + ad)
        @inferred simplify(a * ad)
        @inferred simplify(ad * a)

        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @inferred simplify(σ12 * σ23)
        @inferred simplify(σ12 + σ23)
        @inferred expand_completeness(simplify(σ12))

        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        @inferred simplify(σx * σy)
        @inferred simplify(σx * σx)
        @inferred simplify(σx + σy)

        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        @inferred simplify(Sz * Sy)
        @inferred simplify(Sx + Sy)

        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @inferred simplify(p * x)
        @inferred simplify(x * p)
        @inferred simplify(x + p)
    end

    @testset "Allocations — Fock" begin
        # Warmup
        simplify(a * ad)
        simplify(ad * a)
        simplify(a + ad)

        allocs_swap = @allocations simplify(a * ad)
        allocs_ordered = @allocations simplify(ad * a)
        allocs_add = @allocations simplify(a + ad)
        allocs_sym = @allocations simplify(a)

        @test allocs_swap < 1000
        @test allocs_ordered < 1000
        @test allocs_add < 1000
        @test allocs_sym < 200

        # Scaling: (a*a')^n allocations should grow, but not explosively
        simplify((a * ad)^2)
        simplify((a * ad)^3)
        allocs_2 = @allocations simplify((a * ad)^2)
        allocs_3 = @allocations simplify((a * ad)^3)
        @test allocs_2 < 3000
        @test allocs_3 < 10000
    end

    @testset "Allocations — Transition" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        σ31 = Transition(hn, :σ, 3, 1)

        # Warmup
        simplify(σ12 * σ23)
        simplify(σ12 * σ31)
        simplify(σ12 + σ23)

        @test @allocations(simplify(σ12 * σ23)) < 1000  # composition
        @test @allocations(simplify(σ12 * σ31)) < 1000  # orthogonality → 0
        @test @allocations(simplify(σ12 + σ23)) < 1000
    end

    @testset "Allocations — Pauli" begin
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        σz = Pauli(hp, :σ, 3)

        # Warmup
        simplify(σx * σx)
        simplify(σx * σy)
        simplify(σx * σy * σz)

        @test @allocations(simplify(σx * σx)) < 1000   # σ² = I
        @test @allocations(simplify(σx * σy)) < 1000   # product rule
        @test @allocations(simplify(σx * σy * σz)) < 1000
    end

    @testset "Allocations — Spin" begin
        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)

        # Warmup
        simplify(Sz * Sy)
        simplify(Sz * Sy * Sx)

        @test @allocations(simplify(Sz * Sy)) < 3000
        @test @allocations(simplify(Sz * Sy * Sx)) < 3000
    end

    @testset "Allocations — PhaseSpace" begin
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)

        # Warmup
        simplify(p * x)
        simplify(x * p)

        @test @allocations(simplify(p * x)) < 1000   # P·X → X·P - i
        @test @allocations(simplify(x * p)) < 1000   # already ordered
    end

    @testset "Allocations — collect like terms" begin
        # Many duplicate terms: Dict-based collection should handle efficiently
        many = sum(a' * a for _ in 1:20) + sum(a * a' for _ in 1:20)
        simplify(many)  # warmup
        @test @allocations(simplify(many)) < 50000
    end

    # ====================================================================
    #  Eager-ordering invariants
    # ====================================================================

    @testset "Eager `*` produces canonical form" begin
        # Fock
        @test length(a * ad) == 2          # a†a + 1
        @test length(ad * a) == 1
        @test length((a * ad)^2) >= 3
        @test isequal(simplify(a * ad), a * ad)
        @test isequal(normal_order(a * ad), a * ad)

        # Transition: |i⟩⟨j|·|j⟩⟨k| → |i⟩⟨k|
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @test length(σ12 * σ23) == 1
        @test isequal(normal_order(σ12 * σ23), σ12 * σ23)

        # Pauli: σx·σy = iσz
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        result = σx * σy
        @test prefactor(only(sorted_arguments(result))) == im
        @test operators(only(sorted_arguments(result)))[1].l1 == 3

        # Spin: out-of-order axes swap eagerly
        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        @test length(Sy * Sx) >= 2
        @test length(Sx * Sy) == 1

        # PhaseSpace: p·x = x·p - i
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @test length(p * x) == 2
        @test length(x * p) == 1
    end

    @testset "Idempotency of normal_order and simplify" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @test isequal(simplify(simplify(σ12 * σ23)), simplify(σ12 * σ23))
        @test isequal(normal_order(normal_order(a * ad)), normal_order(a * ad))

        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        @test isequal(normal_order(normal_order(Sy * Sx)), normal_order(Sy * Sx))

        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @test isequal(normal_order(normal_order(p * x)), normal_order(p * x))
    end

    # ====================================================================
    #  normal_to_symmetric / symmetric_to_normal
    # ====================================================================

    @testset "normal_to_symmetric: basic Fock" begin
        # a†a + 1 → a†a + 1/2
        no = normal_order(a * ad)
        weyl = normal_to_symmetric(no)
        @test length(weyl) == 2
        @test weyl[Op[]] == _to_cnum(1 // 2)
        @test weyl[Op[ad, a]] == _to_cnum(1)
    end

    @testset "symmetric_to_normal: basic Fock" begin
        # a†a + 1/2 → a†a + 1
        weyl = _single_qadd(_to_cnum(1), Op[ad, a]) + _single_qadd(_to_cnum(1 // 2), Op[])
        no = symmetric_to_normal(weyl)
        @test length(no) == 2
        @test no[Op[]] == _to_cnum(1)
        @test no[Op[ad, a]] == _to_cnum(1)
    end

    @testset "normal_to_symmetric: higher powers" begin
        no = normal_order((a * ad)^2)
        weyl = normal_to_symmetric(no)
        back = symmetric_to_normal(weyl)
        @test isequal(back, no)

        no3 = normal_order((a * ad)^3)
        weyl3 = normal_to_symmetric(no3)
        back3 = symmetric_to_normal(weyl3)
        @test isequal(back3, no3)
    end

    @testset "Roundtrip: normal → symmetric → normal (Fock)" begin
        for n in 1:5
            no = normal_order((a * ad)^n)
            @test isequal(symmetric_to_normal(normal_to_symmetric(no)), no)
        end

        # a†·a (already normal ordered, single term)
        no = normal_order(ad * a)
        @test isequal(symmetric_to_normal(normal_to_symmetric(no)), no)
    end

    @testset "Roundtrip: symmetric → normal → symmetric (Fock)" begin
        no = normal_order(a * ad)
        weyl = normal_to_symmetric(no)
        @test isequal(normal_to_symmetric(symmetric_to_normal(weyl)), weyl)

        no2 = normal_order((a * ad)^2)
        weyl2 = normal_to_symmetric(no2)
        @test isequal(normal_to_symmetric(symmetric_to_normal(weyl2)), weyl2)
    end

    @testset "Roundtrip: Multi-mode Fock" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)

        no = normal_order(a1 * a2' * a1' * a2)
        @test isequal(symmetric_to_normal(normal_to_symmetric(no)), no)

        h3 = FockSpace(:a) ⊗ FockSpace(:b) ⊗ FockSpace(:c)
        b1 = Destroy(h3, :a, 1)
        b2 = Destroy(h3, :b, 2)
        b3 = Destroy(h3, :c, 3)
        no = normal_order(b1 * b1' * b2 * b2' * b3 * b3')
        @test isequal(symmetric_to_normal(normal_to_symmetric(no)), no)
    end

    @testset "Roundtrip: Jaynes-Cummings (Fock + NLevel)" begin
        hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a_jc = Destroy(hjc, :a, 1)
        σ12 = Transition(hjc, :σ, 1, 2, 2)
        σ21 = σ12'

        H = a_jc' * σ12 + a_jc * σ21
        H2 = H * H

        no = normal_order(H2)
        @test isequal(symmetric_to_normal(normal_to_symmetric(no)), no)

        comm = commutator(H, a_jc' * a_jc)
        no_comm = normal_order(comm)
        @test isequal(symmetric_to_normal(normal_to_symmetric(no_comm)), no_comm)
    end

    @testset "normal_to_symmetric: non-Fock operators unchanged" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        no = normal_order(σ12 * σ23)
        @test isequal(normal_to_symmetric(no), no)

        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        σy = Pauli(hp, :σ, 2)
        no = normal_order(σx * σy)
        @test isequal(normal_to_symmetric(no), no)

        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        no = normal_order(Sy * Sx)
        @test isequal(normal_to_symmetric(no), no)
    end

    @testset "normal_to_symmetric: QSym input" begin
        @test isequal(normal_to_symmetric(a), normal_to_symmetric(normal_order(a)))
        @test isequal(symmetric_to_normal(a), symmetric_to_normal(normal_order(a)))
    end

    @testset "Type stability: Weyl ordering" begin
        # Fock
        @inferred normal_to_symmetric(a)
        @inferred normal_to_symmetric(a' * a)
        @inferred normal_to_symmetric((a * a')^2)
        @inferred symmetric_to_normal(a)
        @inferred symmetric_to_normal(a' * a)
        @inferred symmetric_to_normal(normal_to_symmetric(a' * a))

        # PhaseSpace
        hps = PhaseSpace(:q)
        x = Position(hps, :x)
        p = Momentum(hps, :p)
        @inferred normal_to_symmetric(x)
        @inferred normal_to_symmetric(x * p)
        @inferred normal_to_symmetric(x * x * p * p)
        @inferred symmetric_to_normal(x * p)
        @inferred symmetric_to_normal(normal_to_symmetric(x * p))

        # Mixed Fock + PhaseSpace
        hmix = FockSpace(:c) ⊗ PhaseSpace(:osc)
        ac = Destroy(hmix, :a, 1)
        xx = Position(hmix, :x, 2)
        pp = Momentum(hmix, :p, 2)
        @inferred normal_to_symmetric(ac' * ac * xx * pp)
        @inferred symmetric_to_normal(ac' * ac * xx * pp)

        # Non-pair operators (no Heisenberg pair to expand)
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)
        σ23 = Transition(hn, :σ, 2, 3)
        @inferred normal_to_symmetric(σ12 * σ23)
        @inferred symmetric_to_normal(σ12 * σ23)

        hpa = PauliSpace(:p)
        σx = Pauli(hpa, :σ, 1)
        @inferred normal_to_symmetric(σx)
        @inferred symmetric_to_normal(σx)

        hs = SpinSpace(:s)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        @inferred normal_to_symmetric(Sx * Sy)
        @inferred symmetric_to_normal(Sx * Sy)
    end

    # ====================================================================
    #  Complex multi-space systems
    # ====================================================================

    @testset "Jaynes-Cummings" begin
        hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a_jc = Destroy(hjc, :a, 1)
        σ12 = Transition(hjc, :σ, 1, 2, 2)
        σ21 = σ12'

        H = a_jc' * σ12 + a_jc * σ21
        @test H isa QAdd

        H2 = H * H
        result = simplify(H2)
        @test result isa QAdd
        @test length(result) >= 2

        no_result = normal_order(H2)
        @test no_result isa QAdd
        @test isequal(normal_order(no_result), no_result)

        σ11 = Transition(hjc, :σ, 1, 1, 2)
        gs = expand_completeness(normal_order(σ11))
        @test gs isa QAdd
        @test length(gs) == 2
    end

    @testset "Mixed Fock + Spin + Pauli" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p) ⊗ SpinSpace(:s)
        a_op = Destroy(h, :a, 1)
        σx = Pauli(h, :σ, 1, 2)
        Sz = Spin(h, :S, 3, 3)

        prod = a_op * σx * Sz
        @test length(prod) == 1

        pauli_sq = σx * σx
        result = simplify(pauli_sq)
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))

        fock_prod = a_op * a_op'
        no_result = normal_order(fock_prod)
        @test length(no_result) == 2
    end
end
