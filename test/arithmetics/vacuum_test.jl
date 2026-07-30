using SecondQuantizedAlgebra
using QuantumOpticsBase
using Latexify: latexify
using Symbolics: Symbolics, @variables
using Test

# Reduce a `Complex{Num}` result to a plain number for comparison against the
# numeric backend.
function num_value(z)
    re = Symbolics.value(Symbolics.unwrap(real(z)))
    im_ = Symbolics.value(Symbolics.unwrap(imag(z)))
    return ComplexF64(re, im_)
end

@testset "vacuum expectation" begin
    @testset "Fock" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)

        @test expect(a * a', Vacuum()) == 1
        @test expect(a' * a, Vacuum()) == 0
        @test expect(a, Vacuum()) == 0
        @test expect(a', Vacuum()) == 0
        @test expect(one(a), Vacuum()) == 1

        for n in 1:5
            @test expect(a^n * (a')^n, Vacuum()) == factorial(n)
        end
        # Off-diagonal Fock matrix elements ⟨m|a^k|n⟩ vanish unless m + k == n.
        @test expect(a^2 * (a')^3, Vacuum()) == 0
        @test expect(a^3 * a * (a')^4, Vacuum()) == factorial(4)
    end

    @testset "leaves and scalars" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        @test expect(3, Vacuum()) == 3
        @test expect(2im, Vacuum()) == 2im
        @test expect(zero(a), Vacuum()) == 0
        @test expect(a * a', Vacuum()) isa Complex{Symbolics.Num}
    end

    @testset "symbolic coefficients" begin
        @variables g Δ
        h = FockSpace(:cav) ⊗ NLevelSpace(:atom, 2)
        a = Destroy(h, :a)
        σ(i, j) = Transition(h, :σ, i, j)
        H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
        @test expect(H, Vacuum()) == 0
        @test isequal(expect(H + Δ * a * a', Vacuum()), Complex(Δ, Symbolics.Num(0)))
    end

    @testset "NLevel" begin
        h = NLevelSpace(:atom, 3)
        σ(i, j) = Transition(h, :σ, i, j)

        @test expect(σ(1, 1), Vacuum()) == 1
        @test expect(σ(2, 2), Vacuum()) == 0
        @test expect(σ(3, 3), Vacuum()) == 0
        @test expect(σ(1, 3), Vacuum()) == 0
        @test expect(σ(1, 2) * σ(2, 1), Vacuum()) == 1
        @test expect(σ(2, 1) * σ(1, 2), Vacuum()) == 0
        @test expect(σ(1, 2) * σ(3, 1), Vacuum()) == 0

        # The reference state follows the declared ground state, not level 1.
        hg = NLevelSpace(:atom, 3, 2)
        τ(i, j) = Transition(hg, :τ, i, j)
        @test expect(τ(2, 2), Vacuum()) == 1
        @test expect(τ(1, 1), Vacuum()) == 0
    end

    @testset "Pauli" begin
        h = PauliSpace(:s)
        σx, σy, σz = Pauli(h, :σ, 1), Pauli(h, :σ, 2), Pauli(h, :σ, 3)

        @test expect(σx, Vacuum()) == 0
        @test expect(σy, Vacuum()) == 0
        @test expect(σz, Vacuum()) == 1
        @test expect(σx * σx, Vacuum()) == 1
        @test expect(σx * σy, Vacuum()) == im
        @test expect(σy * σx, Vacuum()) == -im
        @test expect(σx * σz * σx, Vacuum()) == -1
    end

    @testset "PhaseSpace" begin
        h = PhaseSpace(:osc)
        x, p = Position(h, :x), Momentum(h, :p)

        @test expect(x, Vacuum()) == 0
        @test expect(p, Vacuum()) == 0
        @test expect(x^2, Vacuum()) == 1 // 2
        @test expect(p^2, Vacuum()) == 1 // 2
        @test expect(x * p, Vacuum()) == im / 2
        @test expect(p * x, Vacuum()) == -im / 2
        @test expect(x^4, Vacuum()) == 3 // 4
        @test expect(x^6, Vacuum()) == 15 // 8
        @test expect((x^2 + p^2) / 2, Vacuum()) == 1 // 2   # zero-point energy
        @test expect(im * (x * p - p * x), Vacuum()) == -1
    end

    @testset "Spin" begin
        h = SpinSpace(:S)
        Sx, Sy, Sz = Spin(h, :S, 1), Spin(h, :S, 2), Spin(h, :S, 3)

        for S in (1 // 2, 1, 3 // 2, 2)
            @test expect(Sz, Vacuum(); spin = S) == S
            @test expect(Sx, Vacuum(); spin = S) == 0
            @test expect(Sy, Vacuum(); spin = S) == 0
            @test expect(Sx^2, Vacuum(); spin = S) == S / 2
            @test expect(Sx^2 + Sy^2 + Sz^2, Vacuum(); spin = S) == S * (S + 1)
            @test expect(Sx * Sy, Vacuum(); spin = S) == im * S / 2
        end

        # A `Dict` keyed by subspace position covers mixed spin sizes.
        hp = SpinSpace(:A) ⊗ SpinSpace(:B)
        SA, SB = Spin(hp, :SA, 3, 1), Spin(hp, :SB, 3, 2)
        @test expect(SA * SB, Vacuum(); spin = Dict(1 => 1 // 2, 2 => 3 // 2)) == 3 // 4
    end

    @testset "collective transitions" begin
        h = CollectiveNLevelSpace(:c, 3)
        S(i, j) = CollectiveTransition(h, :S, i, j)

        for N in (1, 4, 7)
            @test expect(S(1, 1), Vacuum(); particles = N) == N
            @test expect(S(2, 2), Vacuum(); particles = N) == 0
            @test expect(S(1, 2), Vacuum(); particles = N) == 0
            # Σᵢ Sⁱⁱ = N·I collectively, not the identity.
            @test expect(S(1, 1) + S(2, 2) + S(3, 3), Vacuum(); particles = N) == N
            # S²¹|g⟩ has norm² N; S¹²|g⟩ vanishes.
            @test expect(S(1, 2) * S(2, 1), Vacuum(); particles = N) == N
            @test expect(S(2, 1) * S(1, 2), Vacuum(); particles = N) == 0
        end
    end

    @testset "indexed sums" begin
        ha = NLevelSpace(:atom, 2)
        i = Index(ha, :i, 5, ha)
        σ(x, y, k) = IndexedOperator(Transition(ha, :σ, x, y), k)
        @test expect(Σ(σ(1, 1, i), i), Vacuum()) == 5
        @test expect(Σ(σ(2, 2, i), i), Vacuum()) == 0
        @test expect(Σ(σ(1, 1, i) + σ(2, 2, i), i), Vacuum()) == 5

        hf = FockSpace(:m)
        j, k = Index(hf, :j, 7, hf), Index(hf, :k, 7, hf)
        aj, ak = IndexedOperator(Destroy(hf, :a), j), IndexedOperator(Destroy(hf, :a), k)
        @test expect(Σ(aj' * aj, j), Vacuum()) == 0
        # The `+1` residual already carries the range; it must not be scaled twice.
        @test expect(Σ(aj * aj', j), Vacuum()) == 7
        @test expect(Σ(Σ(aj * ak', j), k), Vacuum()) == 7
    end

    @testset "sums that keep a symbolic scope" begin
        ha = NLevelSpace(:atom, 2)
        i, j = Index(ha, :i, 5, ha), Index(ha, :j, 5, ha)
        σ(x, y, k) = IndexedOperator(Transition(ha, :σ, x, y), k)

        # An index-dependent coefficient cannot collapse to a range factor.
        Δ = IndexedVariable(:Δ, i)
        r = expect(Σ(Δ * σ(1, 1, i), i), Vacuum())
        @test is_indexed_sum(real(r))
        @test isequal(SecondQuantizedAlgebra.get_sum_indices(real(r)), [i])
        # ...but a vanishing summand still short-circuits to zero.
        @test expect(Σ(Δ * σ(2, 2, i), i), Vacuum()) == 0

        # A `≠` constraint means the count is not the full range.
        c = expect(Σ(σ(1, 1, i), i, [j]), Vacuum())
        @test is_indexed_sum(real(c))
        @test isequal(SecondQuantizedAlgebra.get_sum_non_equal(real(c)), [(i, j)])
    end

    @testset "unresolved free indices" begin
        ha = NLevelSpace(:atom, 3)
        i, j = Index(ha, :i, 5, ha), Index(ha, :j, 5, ha)
        σ(x, y, k) = IndexedOperator(Transition(ha, :σ, x, y), k)

        # Both readings agree, so no delta survives.
        @test expect(σ(1, 1, i) * σ(1, 1, j), Vacuum()) == 1
        @test expect(σ(2, 2, i) * σ(2, 2, j), Vacuum()) == 0
        @test expect(σ(1, 2, i) * σ(2, 3, j), Vacuum()) == 0

        # They disagree, so the exact case split is returned.
        q = σ(1, 2, i) * σ(2, 1, j)
        r = expect(q, Vacuum())
        @test isequal(r, Complex(kronecker_delta(i, j), Symbolics.Num(0)))
        # The two branches it encodes are the two explicit readings.
        @test isequal(
            Symbolics.substitute(r, Dict(index_sym(j) => index_sym(i))),
            expect(change_index(q, j, i), Vacuum()),
        )
        @test expect(assume_distinct_index(q, [(i, j)]), Vacuum()) == 0

        @test isequal(
            expect(q + 2, Vacuum()),
            Complex(2 + kronecker_delta(i, j), Symbolics.Num(0)),
        )
    end

    @testset "kronecker_delta" begin
        h = NLevelSpace(:atom, 2)
        i, j = Index(h, :i, 5, h), Index(h, :j, 5, h)

        @test kronecker_delta(i, i) == 1
        @test is_kronecker_delta(kronecker_delta(i, j))
        @test !is_kronecker_delta(kronecker_delta(i, i))
        @test !is_kronecker_delta(1)
        # Symmetric: one canonical node per unordered pair, so like terms collect.
        @test isequal(kronecker_delta(i, j), kronecker_delta(j, i))
        @test isequal(kronecker_delta(i, j) - kronecker_delta(j, i), 0)
        # Substituting one index for the other folds the node away.
        @test isequal(
            Symbolics.substitute(kronecker_delta(i, j), Dict(index_sym(j) => index_sym(i))),
            1,
        )
        @test sprint(show, kronecker_delta(i, j)) == "δ(i, j)"
        @test occursin("\\delta_{i,j}", string(latexify(kronecker_delta(i, j))))
    end

    @testset "matrix elements and overlaps" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        # ⟨ψ|A|φ⟩ = ⟨0|Oψ† A Oφ|0⟩ with |n⟩ ∝ (a†)ⁿ|0⟩.
        norm2(n) = expect(a^n * (a')^n, Vacuum())
        @test norm2(3) == 6
        @test expect(a^2 * (a' * a) * (a')^2, Vacuum()) / norm2(2) == 2   # ⟨2|a†a|2⟩
        @test expect(a^1 * (a')^2, Vacuum()) == 0                        # ⟨1|2⟩
    end

    @testset "missing or invalid sizes" begin
        hs = SpinSpace(:s)
        @test_throws ArgumentError expect(Spin(hs, :S, 3), Vacuum())
        @test_throws ArgumentError expect(Spin(hs, :S, 3), Vacuum(); spin = 1 // 3)
        @test_throws ArgumentError expect(Spin(hs, :S, 3), Vacuum(); spin = -1)
        @test_throws ArgumentError expect(Spin(hs, :S, 3), Vacuum(); spin = Dict(2 => 1 // 2))

        hc = CollectiveNLevelSpace(:c, 2)
        S11 = CollectiveTransition(hc, :S, 1, 1)
        @test_throws ArgumentError expect(S11, Vacuum())
        @test_throws ArgumentError expect(S11, Vacuum(); particles = 3 // 2)
        @test_throws ArgumentError expect(S11, Vacuum(); particles = -1)

        # A bare operator against an indexed one on the same subspace has no site
        # label to identify, so there is nothing to case-split on.
        hf = FockSpace(:m)
        k = Index(hf, :k, 5, hf)
        a = Destroy(hf, :a)
        @test_throws ArgumentError expect(a' * IndexedOperator(a, k), Vacuum())
    end

    @testset "agrees with the numeric backend" begin
        h = FockSpace(:cav) ⊗ NLevelSpace(:atom, 2) ⊗ PauliSpace(:s) ⊗ PhaseSpace(:osc)
        a = Destroy(h, :a)
        σ(i, j) = Transition(h, :σ, i, j)
        σx, σz = Pauli(h, :σp, 1), Pauli(h, :σp, 3)
        x, p = Position(h, :x), Momentum(h, :p)

        b = FockBasis(14) ⊗ NLevelBasis(2) ⊗ SpinBasis(1 // 2) ⊗ FockBasis(14)
        ψ = basisstate(b, [1, 1, 1, 1])

        ops = [
            a * a', a' * a, σ(1, 2) * σ(2, 1), σz, σx, x^2, x * p,
            a * a' * σ(1, 2) * σ(2, 1) * σz * x^2, (a + a')^4, (x + p)^2,
            a^3 * (a')^3, x^4 * p^2, σ(1, 2) * a' + σ(2, 1) * a, σx * σz * σx, (x + p)^6,
        ]
        for op in ops
            @test num_value(expect(op, Vacuum())) ≈ expect(to_numeric(op, b), ψ) atol = 1.0e-8
        end
    end

    @testset "agrees with the numeric backend (collective)" begin
        nlev = 3
        h = CollectiveNLevelSpace(:c, nlev)
        S(i, j) = CollectiveTransition(h, :S, i, j)
        ops = [
            S(1, 1), S(2, 2), S(1, 2), S(2, 1), S(1, 2) * S(2, 1), S(2, 1) * S(1, 2),
            S(1, 2) * S(2, 1) * S(1, 2) * S(2, 1), S(1, 3) * S(3, 1), S(2, 3) * S(3, 2),
            S(2, 1) * S(1, 3) * S(3, 2), S(1, 2) * S(2, 3) * S(3, 1),
            (S(1, 2) + S(2, 1))^2, (S(1, 2) + S(2, 1))^4, S(1, 1)^3,
        ]
        ob = NLevelBasis(nlev)
        for N in (1, 2, 3, 5)
            b = ManyBodyBasis(ob, bosonstates(ob, N))
            ψ = basisstate(b, findfirst(o -> o == [N, 0, 0], b.occupations))
            for op in ops
                @test num_value(expect(op, Vacuum(); particles = N)) ≈
                    expect(to_numeric(op, b), ψ) atol = 1.0e-8
            end
        end
    end

    @testset "agrees with the numeric backend (spin)" begin
        h = SpinSpace(:S)
        Sx, Sy, Sz = Spin(h, :S, 1), Spin(h, :S, 2), Spin(h, :S, 3)
        ops = [
            Sx, Sy, Sz, Sx^2, Sy^2, Sz^2, Sx^2 + Sy^2 + Sz^2, Sx * Sy, Sy * Sx,
            Sx^4, Sx * Sz * Sx, (Sx + Sy)^3, Sx * Sy * Sz, (Sx + Sz)^5,
        ]
        for S in (1 // 2, 1, 3 // 2, 2, 5 // 2)
            b = SpinBasis(S)
            ψ = basisstate(b, 1)
            for op in ops
                @test num_value(expect(op, Vacuum(); spin = S)) ≈
                    expect(to_numeric(op, b), ψ) atol = 1.0e-8
            end
        end
    end
end
