using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, _single_qadd, _to_cnum

@testset "latexify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Operators" begin
        hn = NLevelSpace(:atom, 3, 1)
        hp = PauliSpace(:p)
        hs = SpinSpace(:s)
        hps = PhaseSpace(:q)

        input = [
            Destroy(h, :a),
            Create(h, :a),
            Transition(hn, :σ, 1, 2),
            Pauli(hp, :σ, 1),
            Pauli(hp, :σ, 2),
            Pauli(hp, :σ, 3),
            Spin(hs, :S, 1),
            Spin(hs, :S, 3),
            Position(hps, :x),
            Momentum(hps, :p),
        ]
        output = [
            L"a",
            L"a^{\dagger}",
            L"{\sigma}^{{12}}",
            L"{\sigma}_{{x}}",
            L"{\sigma}_{{y}}",
            L"{\sigma}_{{z}}",
            L"{S}_{{x}}",
            L"{S}_{{z}}",
            L"\hat{x}",
            L"\hat{p}",
        ]
        for (i, o) in zip(input, output)
            @test latexify(i) == o
        end
    end

    @testset "Transition superscript toggle" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)

        transition_superscript(true)
        @test latexify(σ12) == L"{\sigma}^{{12}}"

        transition_superscript(false)
        @test latexify(σ12) == L"{\sigma}_{{12}}"

        transition_superscript(true)  # reset
    end

    @testset "Products and scalars" begin
        input = [
            1 * ad * a,
            3 * ad * a,
            -1 * a,
            _single_qadd(_to_cnum(5), QSym[]),
        ]
        output = [
            L"a^{\dagger}a",
            L"3 a^{\dagger}a",
            L"-a",
            L"5",
        ]
        for (i, o) in zip(input, output)
            @test latexify(i) == o
        end
    end

    @testset "QAdd" begin
        @test latexify(a + ad) == L"a + a^{\dagger}"
    end

    @testset "Simplify display" begin
        @test latexify(simplify(a * ad)) == L"1 + a^{\dagger}a"
    end

    @testset "Symbolic prefactors" begin
        @variables g
        @test latexify(g * a) == L"g a"
    end

    @testset "Complex prefactors" begin
        result = simplify(Pauli(PauliSpace(:p), :σ, 1) * Pauli(PauliSpace(:p), :σ, 2))
        @test latexify(result) == L"\mathit{i} {\sigma}_{{z}}"
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
            L"{\sigma}_{i}^{{12}}",
            L"a_{i}",
            L"a_{i}^{\dagger}",
        ]
        for (inp, o) in zip(input, output)
            @test latexify(inp) == o
        end
    end

    @testset "Indexed product" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        @test latexify(IndexedVariable(:g, i) * b' * σ_i) ==
            L"g\left( i \right) b^{\dagger}{\sigma}_{i}^{{12}}"
    end

    @testset "Sum" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        H = Σ(IndexedVariable(:g, i) * b' * σ_i, i)
        @test latexify(H) ==
            L"\underset{i}{\overset{N}{\sum}}g\left( i \right) b^{\dagger}{\sigma}_{i}^{{12}}"
    end

    @testset "Sum separates index-independent terms" begin
        @variables N Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        gi = IndexedVariable(:g, i)

        H = Δ * b' * b + Σ(gi * (b * σ_i + b' * σ_i'), i)
        latex_str = latexify(H)
        # Δb†b should be outside the sum, indexed terms inside
        @test latex_str ==
            L"\Delta b^{\dagger}b + \underset{i}{\overset{N}{\sum}}\left( g\left( i \right) b{\sigma}_{i}^{{12}} + g\left( i \right) b^{\dagger}{\sigma}_{i}^{{21}} \right)"
    end

    @testset "Sum with non_equal constraint uses \\neq" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        j = Index(h2, :j, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        Γij = DoubleIndexedVariable(:Γ, i, j)

        s = Σ(Γij * σ_i, i, [j])
        @test latexify(s) ==
            L"\underset{i{\neq}i{\neq}j}{\overset{N}{\sum}}\Gamma\left( i, j \right) {\sigma}_{i}^{{12}}"
    end

    @testset "Fraction prefactor gets brackets" begin
        @variables g Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, (:g, :e))
        @qnumbers b::Destroy(h2, 1)
        σee = Transition(h2, :σ, 2, 2, 2)

        # Fraction prefactor: (g²/Δ) a†a
        @test latexify((g^2 / Δ) * b' * b) ==
            L"\left(\frac{g^{2}}{\Delta}\right) b^{\dagger}b"

        # Sum prefactor: (g²/Δ + Δ) σ₂₂
        @test latexify((g^2 / Δ + Δ) * σee) ==
            L"\left(\frac{g^{2}}{\Delta} + \Delta\right) {\sigma}^{{22}}"

        # Power prefactor: no brackets needed
        @test latexify(g^2 * b' * b) ==
            L"g^{2} b^{\dagger}b"

        # Simple symbolic: no brackets needed
        @test latexify(g * b) == L"g b"
    end

    @testset "Sum all-indexed terms" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        gi = IndexedVariable(:g, i)

        # All terms depend on index — no separation needed
        @test latexify(Σ(gi * b' * σ_i + gi * b * σ_i', i)) ==
            L"\underset{i}{\overset{N}{\sum}}\left( g\left( i \right) b{\sigma}_{i}^{{21}} + g\left( i \right) b^{\dagger}{\sigma}_{i}^{{12}} \right)"
    end

    @testset "Sum single indexed + independent without parens" begin
        @variables N Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        gi = IndexedVariable(:g, i)

        # Single indexed term — no parentheses around the sum body
        @test latexify(Δ * b' * b + Σ(gi * b' * σ_i, i)) ==
            L"\Delta b^{\dagger}b + \underset{i}{\overset{N}{\sum}}g\left( i \right) b^{\dagger}{\sigma}_{i}^{{12}}"
    end

    @testset "MIME text/latex" begin
        @test repr(MIME"text/latex"(), a) == latexify(a)
        @test repr(MIME"text/latex"(), ad) == latexify(ad)
        @test repr(MIME"text/latex"(), a * ad) == latexify(a * ad)
        @test repr(MIME"text/latex"(), a + ad) == latexify(a + ad)
    end
end
