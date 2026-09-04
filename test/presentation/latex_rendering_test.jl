using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, transition_superscript, make_time_dependent

hf = FockSpace(:c)
af = Destroy(hf, :a)
adf = af'

@testset "LaTeX (latexify)" begin
    @testset "Operators" begin
        hn = NLevelSpace(:atom, 3, 1)
        hp = PauliSpace(:p)
        hs = SpinSpace(:s)
        hps = PhaseSpace(:q)

        cases = [
            (Destroy(hf, :a), L"a"),
            (Create(hf, :a), L"a^{\dagger}"),
            (Transition(hn, :σ, 1, 2), L"{\sigma}^{{12}}"),
            (Pauli(hp, :σ, 1), L"{\sigma}_{{x}}"),
            (Pauli(hp, :σ, 2), L"{\sigma}_{{y}}"),
            (Pauli(hp, :σ, 3), L"{\sigma}_{{z}}"),
            (Spin(hs, :S, 1), L"{S}_{{x}}"),
            (Spin(hs, :S, 3), L"{S}_{{z}}"),
            (Position(hps, :x), L"\hat{x}"),
            (Momentum(hps, :p), L"\hat{p}"),
        ]
        for (input, out) in cases
            @test latexify(input) == out
        end
    end

    @testset "Compound names render as subscript label" begin
        hn = NLevelSpace(:atom, 3, 1)
        hp = PauliSpace(:p)
        hs = SpinSpace(:s)
        hps = PhaseSpace(:q)

        # Raw `a_pol` would be parsed by KaTeX as `a` with subscript `p`,
        # leaving `ol` as stray text. We emit `a_{\mathrm{pol}}` instead.
        cases = [
            (Destroy(hf, :a_pol), L"a_{\mathrm{pol}}"),
            (Create(hf, :c_bog), L"c_{\mathrm{bog}}^{\dagger}"),
            (Destroy(hf, :b_long_name), L"b_{\mathrm{long\_name}}"),
            (Transition(hn, :σ_lab, 1, 2), L"{\sigma_{\mathrm{lab}}}^{{12}}"),
            (Pauli(hp, :τ_atom, 2), L"{\tau_{\mathrm{atom}}}_{{y}}"),
            (Spin(hs, :S_lab, 3), L"{S_{\mathrm{lab}}}_{{z}}"),
            (Position(hps, :x_com), L"\hat{x_{\mathrm{com}}}"),
            (Momentum(hps, :p_rel), L"\hat{p_{\mathrm{rel}}}"),
        ]
        for (input, out) in cases
            @test latexify(input) == out
        end
    end

    @testset "Index slot suffixes render as comma subscript" begin
        # `i_2_1` (accumulated slot suffixes) must render as `i_{2,1}`, not the
        # invalid double subscript `_{i_2_1}`; bare names are unchanged.
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        iu = Index(h2, Symbol("i_2_1"), N, NLevelSpace(:atom, 2, 1))

        cases = [
            (IndexedOperator(Transition(h2, :σ, 1, 2, 2), i), L"{\sigma}_{i}^{{12}}"),
            (IndexedOperator(Transition(h2, :σ, 1, 2, 2), iu), L"{\sigma}_{i_{2,1}}^{{12}}"),
        ]
        for (input, out) in cases
            @test latexify(input) == out
        end
    end

    @testset "Transition superscript toggle" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)

        transition_superscript(true)
        @test latexify(σ12) == L"{\sigma}^{{12}}"

        transition_superscript(false)
        @test latexify(σ12) == L"{\sigma}_{{12}}"

        transition_superscript(true)
    end

    @testset "Products and scalars" begin
        cases = [
            (1 * adf * af, L"a^{\dagger}a"),
            (3 * adf * af, L"3 a^{\dagger}a"),
            (-1 * af, L"-a"),
            (5 * commutator(af, adf), L"5"),
        ]
        for (input, out) in cases
            @test latexify(input) == out
        end
    end

    @testset "QAdd" begin
        @test latexify(af + adf) == L"a + a^{\dagger}"
    end

    @testset "Simplify display" begin
        @test latexify(simplify(af * adf)) == L"1 + a^{\dagger}a"
    end

    @testset "Symbolic prefactors" begin
        @variables g
        @test latexify(g * af) == L"g a"
    end

    @testset "Complex prefactors" begin
        result = simplify(Pauli(PauliSpace(:p), :σ, 1) * Pauli(PauliSpace(:p), :σ, 2))
        @test latexify(result) == L"\mathit{i} {\sigma}_{{z}}"
    end

    @testset "Symbolic imaginary and complex prefactors" begin
        # Regression guard: a `Coeff` with a symbolic imaginary part must lower
        # to `Complex{Num}` in `latex_prefactor`, so rendering (including the
        # text/latex MIME path Documenter uses for `@example` output) never hits
        # `needs_pf_brackets(::Coeff)`.
        @variables g κ
        hL = FockSpace(:f)
        aL = Destroy(hL, :a)

        # native scaled pure-imaginary: the complex(false, i_val) branch
        @test latexify(2im * aL) == L"2\mathit{i} a"
        # symbolic pure-imaginary prefactor (the r_is_zero, non-Real-imag fallback)
        @test latexify((im * g) * aL) == L"g ~ \mathit{i} a"
        @test repr(MIME"text/latex"(), (im * g) * aL) == latexify((im * g) * aL)
        # symbolic complex prefactor: real and imaginary parts both symbolic
        @test latexify((g + im * κ) * aL) == L"g + \kappa ~ \mathit{i} a"

        # the exact shape that reached docs CI: a Spin commutator yields i * g * Sz
        Sx = Spin(SpinSpace(:S), :S, 1)
        Sy = Spin(SpinSpace(:S), :S, 2)
        @test repr(MIME"text/latex"(), commutator(g * Sx, Sy)) == L"g ~ \mathit{i} {S}_{{z}}"
    end

    @testset "Indexed operators" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))

        cases = [
            (IndexedOperator(Transition(h2, :σ, 1, 2, 2), i), L"{\sigma}_{i}^{{12}}"),
            (IndexedOperator(Destroy(h2, :a, 1), i), L"a_{i}"),
            (IndexedOperator(Destroy(h2, :a, 1), i)', L"a_{i}^{\dagger}"),
        ]
        for (input, out) in cases
            @test latexify(input) == out
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
        @test latexify(H) ==
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
            L"\underset{i{\neq}j}{\overset{N}{\sum}}\Gamma\left( i, j \right) {\sigma}_{i}^{{12}}"

        kk = Index(h2, :k, N, NLevelSpace(:atom, 2, 1))
        ll = Index(h2, :l, N, NLevelSpace(:atom, 2, 1))
        σ_k = IndexedOperator(Transition(h2, :σ, 2, 2, 2), kk)
        σ_l = IndexedOperator(Transition(h2, :σ, 2, 1, 2), ll)
        s2 = Σ(σ_k, kk) * σ_l
        @test latexify(s2) ==
            L"{\sigma}_{l}^{{21}} + \underset{k{\neq}l}{\overset{N}{\sum}}{\sigma}_{k}^{{22}}{\sigma}_{l}^{{21}}"
    end

    @testset "Fraction prefactor gets brackets" begin
        @variables g Δ
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, (:g, :e))
        @qnumbers b::Destroy(h2, 1)
        σee = Transition(h2, :σ, 2, 2, 2)

        @test latexify((g^2 / Δ) * b' * b) ==
            L"\left(\frac{g^{2}}{\Delta}\right) b^{\dagger}b"
        @test latexify((g^2 / Δ + Δ) * σee) ==
            L"\left(\frac{g^{2}}{\Delta} + \Delta\right) {\sigma}^{{22}}"
        @test latexify(g^2 * b' * b) == L"g^{2} b^{\dagger}b"
        @test latexify(g * b) == L"g b"
    end

    @testset "Sum all-indexed terms" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
        gi = IndexedVariable(:g, i)

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

        @test latexify(Δ * b' * b + Σ(gi * b' * σ_i, i)) ==
            L"\Delta b^{\dagger}b + \underset{i}{\overset{N}{\sum}}g\left( i \right) b^{\dagger}{\sigma}_{i}^{{12}}"
    end

    @testset "Per-term Σ scope: LaTeX distinct-index sums stay separate" begin
        @variables N
        hf = FockSpace(:f)
        af = Destroy(hf, :a)
        i = Index(hf, :i, N, hf)
        j = Index(hf, :j, N, hf)
        ai = IndexedOperator(af, i)
        aj_dag = IndexedOperator(af', j)

        expr = Σ(ai, i) + Σ(aj_dag, j)
        @test latexify(expr) ==
            L"\underset{i}{\overset{N}{\sum}}a_{i} + \underset{j}{\overset{N}{\sum}}a_{j}^{\dagger}"
    end

    @testset "MIME text/latex" begin
        @test repr(MIME"text/latex"(), af) == latexify(af)
        @test repr(MIME"text/latex"(), adf) == latexify(adf)
        @test repr(MIME"text/latex"(), af * adf) == latexify(af * adf)
        @test repr(MIME"text/latex"(), af + adf) == latexify(af + adf)
    end

    @testset "Averages" begin
        avg_a = average(af)
        avg_prod = average(adf * af)
        @test string(latexify(avg_a)) ==
            "\\begin{equation}\n\\langle a \\rangle\n\\end{equation}\n"
        @test string(latexify(avg_prod)) ==
            "\\begin{equation}\n\\langle a^{\\dagger}a \\rangle\n\\end{equation}\n"

        h = FockSpace(:site)
        @variables N
        i = Index(h, :i, N, h)
        ai = IndexedOperator(af, i)
        avg_sum = average(Σ(ai, i))
        @test string(latexify(avg_sum)) ==
            "\\begin{equation}\n\\underset{i}{\\overset{N}{\\sum}}\\langle a_{i} \\rangle\n\\end{equation}\n"
    end

    @testset "Lifted (time-dependent) averages" begin
        @variables t N
        h = FockSpace(:site)
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        ai = IndexedOperator(Destroy(h, :a), i)
        σ = Transition(NLevelSpace(:atom, 3, 1), :σ, 1, 2)

        @test string(latexify(make_time_dependent(average(af), t))) ==
            "\\begin{equation}\n\\langle a \\rangle\\left( t \\right)\n\\end{equation}\n"
        @test string(latexify(make_time_dependent(average(adf * af), t))) ==
            "\\begin{equation}\n\\langle a^{\\dagger}a \\rangle\\left( t \\right)\n\\end{equation}\n"
        @test string(latexify(make_time_dependent(average(σ), t))) ==
            "\\begin{equation}\n\\langle {\\sigma}^{{12}} \\rangle\\left( t \\right)\n\\end{equation}\n"
        @test string(latexify(make_time_dependent(average(Σ(ai, i)), t))) ==
            "\\begin{equation}\n\\underset{i}{\\overset{N}{\\sum}}\\langle a_{i} \\rangle\\left( t \\right)\n\\end{equation}\n"
        @test string(latexify(make_time_dependent(average(Σ(ai, i, [j])), t))) ==
            "\\begin{equation}\n\\underset{i{\\neq}j}{\\overset{N}{\\sum}}\\langle a_{i} \\rangle\\left( t \\right)\n\\end{equation}\n"
    end
end
