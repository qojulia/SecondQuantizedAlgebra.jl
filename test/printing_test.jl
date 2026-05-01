using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QTermDict, _to_cnum

@testset "printing" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = a'

    @testset "HilbertSpaces" begin
        input = [
            FockSpace(:cavity),
            FockSpace(:a) ⊗ FockSpace(:b),
            NLevelSpace(:atom, 3, 1),
            PauliSpace(:p),
            SpinSpace(:s),
            PhaseSpace(:q),
        ]
        output = [
            "ℋ(cavity)",
            "ℋ(a) ⊗ ℋ(b)",
            "ℋ(atom)",
            "ℋ(p)",
            "ℋ(s)",
            "ℋ(q)",
        ]
        for (i, o) in zip(input, output)
            @test repr(i) == o
        end
    end

    @testset "Operators" begin
        hn = NLevelSpace(:atom, 3, 1)
        hp = PauliSpace(:p)
        hs = SpinSpace(:s)
        hps = PhaseSpace(:q)

        input = [
            Destroy(h, :a),
            Create(h, :a),
            Transition(hn, :σ, 1, 2),
            Transition(hn, :σ, 3, 1),
            Pauli(hp, :σ, 1),
            Pauli(hp, :σ, 2),
            Pauli(hp, :σ, 3),
            Spin(hs, :S, 1),
            Spin(hs, :S, 2),
            Spin(hs, :S, 3),
            Position(hps, :x),
            Momentum(hps, :p),
        ]
        output = [
            "a",
            "a'",
            "σ₁₂",
            "σ₃₁",
            "σx",
            "σy",
            "σz",
            "Sx",
            "Sy",
            "Sz",
            "x",
            "p",
        ]
        for (i, o) in zip(input, output)
            @test repr(i) == o
        end
    end

    @testset "Single-term QAdd" begin
        input = [
            1 * ad * a,
            3 * ad * a,
            -1 * a,
            -3 * a,
            QAdd(QTermDict(QSym[] => _to_cnum(5)), Index[], Tuple{Index, Index}[]),
            0.5 * a,
        ]
        output = [
            "a' * a",
            "3 * a' * a",
            "-a",
            "-3 * a",
            "5",
            "0.5 * a",
        ]
        for (i, o) in zip(input, output)
            @test repr(i) == o
        end
    end

    @testset "QAdd" begin
        input = [
            a + ad,
            2 * a + 3 * ad,
        ]
        output = [
            "a + a'",
            "2 * a + 3 * a'",
        ]
        for (i, o) in zip(input, output)
            @test repr(i) == o
        end
    end

    @testset "Simplify display" begin
        @test repr(simplify(a * ad)) == "1 + a' * a"
    end

    @testset "Indexed operators" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))

        input = [
            IndexedOperator(Transition(h2, :σ, 1, 2, 2), i),
            IndexedOperator(Destroy(h2, :a, 1), i),
            IndexedOperator(Destroy(h2, :a, 1), i)',
        ]
        output = [
            "σ_i₁₂",
            "a_i",
            "a_i'",
        ]
        for (inp, o) in zip(input, output)
            @test repr(inp) == o
        end
    end

    @testset "Indexed product" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        @test repr(IndexedVariable(:g, i) * b' * σ_i) == "g(i) * b' * σ_i₁₂"
    end

    @testset "Sum display" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        j = Index(h2, :j, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        H = Σ(IndexedVariable(:g, i) * b' * σ_i, i)
        @test repr(H) == "Σ(i=1:N) g(i) * b' * σ_i₁₂"

        # Sum with non-equal constraint
        σ_j = IndexedOperator(Transition(h2, :σ, 2, 1, 2), j)
        S = Σ(σ_i * σ_j, i, [j])
        @test repr(S) == "Σ(i=1:N)(i≠j) σ_i₁₂ * σ_j₂₁"
    end

    @testset "Fraction prefactor gets brackets" begin
        @variables g Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, (:g, :e))
        @qnumbers c::Destroy(h2, 1)
        σee = Transition(h2, :σ, 2, 2, 2)

        @test repr((g^2 / Δ) * c' * c) == "((g^2) / Δ) * c' * c"
        @test repr((g^2 / Δ + Δ) * σee) == "((g^2) / Δ + Δ) * σ₂₂"
        @test repr(g^2 * c' * c) == "g^2 * c' * c"
    end

    @testset "Sum separates index-independent terms" begin
        @variables N Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        gi = IndexedVariable(:g, i)

        # Independent term outside sum, multiple indexed terms in parens
        @test repr(Δ * b' * b + Σ(gi * (b * σ_i + b' * σ_i'), i)) ==
            "Δ * b' * b + Σ(i=1:N) (g(i) * b * σ_i₁₂ + g(i) * b' * σ_i₂₁)"

        # Single indexed term — no parentheses
        @test repr(Δ * b' * b + Σ(gi * b' * σ_i, i)) ==
            "Δ * b' * b + Σ(i=1:N) g(i) * b' * σ_i₁₂"

        # All terms indexed — no separation needed
        @test repr(Σ(gi * b' * σ_i + gi * b * σ_i', i)) ==
            "Σ(i=1:N) (g(i) * b * σ_i₂₁ + g(i) * b' * σ_i₁₂)"
    end

    @testset "Edge cases — no crash" begin
        @test repr(QAdd(QTermDict(), Index[], Tuple{Index, Index}[])) == "0"
    end

    @testset "Type inference" begin
        s = IOBuffer(sizehint = 0)
        @inferred show(s, a)
        @inferred show(s, ad)
    end
end
