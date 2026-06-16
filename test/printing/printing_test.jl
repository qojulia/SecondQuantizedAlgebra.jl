using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, _single_qadd, _zero_qadd, _to_cnum,
    transition_superscript, make_time_dependent

@testset "Rendering" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    ad = a'

    hf = FockSpace(:c)
    af = Destroy(hf, :a)
    adf = af'

    @testset "Unicode (repr)" begin
        @testset "HilbertSpaces" begin
            cases = [
                (FockSpace(:cavity), "ℋ(cavity)"),
                (FockSpace(:a) ⊗ FockSpace(:b), "ℋ(a) ⊗ ℋ(b)"),
                (NLevelSpace(:atom, 3, 1), "ℋ(atom)"),
                (PauliSpace(:p), "ℋ(p)"),
                (SpinSpace(:s), "ℋ(s)"),
                (PhaseSpace(:q), "ℋ(q)"),
            ]
            for (input, out) in cases
                @test repr(input) == out
            end
        end

        @testset "Operators" begin
            hn = NLevelSpace(:atom, 3, 1)
            hp = PauliSpace(:p)
            hs = SpinSpace(:s)
            hps = PhaseSpace(:q)

            cases = [
                (Destroy(h, :a), "a"),
                (Create(h, :a), "a'"),
                (Transition(hn, :σ, 1, 2), "σ₁₂"),
                (Transition(hn, :σ, 3, 1), "σ₃₁"),
                (Pauli(hp, :σ, 1), "σx"),
                (Pauli(hp, :σ, 2), "σy"),
                (Pauli(hp, :σ, 3), "σz"),
                (Spin(hs, :S, 1), "Sx"),
                (Spin(hs, :S, 2), "Sy"),
                (Spin(hs, :S, 3), "Sz"),
                (Position(hps, :x), "x"),
                (Momentum(hps, :p), "p"),
            ]
            for (input, out) in cases
                @test repr(input) == out
            end
        end

        @testset "Single-term QAdd" begin
            cases = [
                (1 * ad * a, "a' * a"),
                (3 * ad * a, "3 * a' * a"),
                (-1 * a, "-a"),
                (-3 * a, "-3 * a"),
                (_single_qadd(_to_cnum(5), QSym[]), "5"),
                (0.5 * a, "0.5 * a"),
            ]
            for (input, out) in cases
                @test repr(input) == out
            end
        end

        @testset "QAdd" begin
            cases = [
                (a + ad, "a + a'"),
                (2 * a + 3 * ad, "2 * a + 3 * a'"),
            ]
            for (input, out) in cases
                @test repr(input) == out
            end
        end

        @testset "Simplify display" begin
            @test repr(simplify(a * ad)) == "1 + a' * a"
        end

        @testset "Indexed operators" begin
            @variables N
            h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
            i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))

            cases = [
                (IndexedOperator(Transition(h2, :σ, 1, 2, 2), i), "σ_i₁₂"),
                (IndexedOperator(Destroy(h2, :a, 1), i), "a_i"),
                (IndexedOperator(Destroy(h2, :a, 1), i)', "a_i'"),
            ]
            for (input, out) in cases
                @test repr(input) == out
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

            @test repr(Δ * b' * b + Σ(gi * (b * σ_i + b' * σ_i'), i)) ==
                "Δ * b' * b + Σ(i=1:N) (g(i) * b * σ_i₁₂ + g(i) * b' * σ_i₂₁)"
            @test repr(Δ * b' * b + Σ(gi * b' * σ_i, i)) ==
                "Δ * b' * b + Σ(i=1:N) g(i) * b' * σ_i₁₂"
            @test repr(Σ(gi * b' * σ_i + gi * b * σ_i', i)) ==
                "Σ(i=1:N) (g(i) * b * σ_i₂₁ + g(i) * b' * σ_i₁₂)"
        end

        @testset "Average and averaged sums" begin
            @variables N
            h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
            @qnumbers b::Destroy(h2, 1)
            i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
            j = Index(h2, :j, N, NLevelSpace(:atom, 2, 1))
            σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)
            σ_j = IndexedOperator(Transition(h2, :σ, 1, 2, 2), j)

            @test repr(average(b)) == "⟨b⟩"
            @test repr(average(b' * b)) == "⟨b' * b⟩"

            @test repr(average(Σ(b' * σ_i, i))) ==
                "Σ(i=1:N) ⟨b' * σ_i₁₂⟩"

            @test repr(average(Σ(σ_i, i, [j]))) ==
                "Σ(i=1:N)(i≠j) ⟨σ_i₁₂⟩"

            @test repr(average(Σ(σ_i * σ_j, i, j))) ==
                "Σ(i=1:N)Σ(j=1:N)(i≠j) ⟨σ_i₁₂ * σ_j₁₂⟩"

            @test repr(2 * average(Σ(b' * σ_i, i))) ==
                "2Σ(i=1:N) ⟨b' * σ_i₁₂⟩"

            @test repr(1 + average(Σ(b' * σ_i, i))) ==
                "1 + Σ(i=1:N) ⟨b' * σ_i₁₂⟩"
        end

        @testset "Lifted (time-dependent) averages" begin
            @variables t N
            hsite = FockSpace(:site)
            asite = Destroy(hsite, :a)
            i = Index(hsite, :i, N, hsite)
            ai = IndexedOperator(asite, i)

            @test repr(make_time_dependent(average(af), t)) == "⟨a⟩(t)"
            @test repr(make_time_dependent(average(adf * af), t)) == "⟨a' * a⟩(t)"
            # the lifted operator carries its own Σ scope; the prefix is not doubled
            @test repr(make_time_dependent(average(Σ(ai, i)), t)) == "⟨Σ(i=1:N) a_i⟩(t)"
        end

        @testset "Scoped constraints stay in separate sum groups" begin
            @variables N
            h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
            @qnumbers b::Destroy(h2, 1)
            i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
            j = Index(h2, :j, N, NLevelSpace(:atom, 2, 1))
            bi = IndexedOperator(b, i)

            expr = Σ(bi, i) + Σ(bi, i, [j])
            @test repr(expr) == "Σ(i=1:N) b_i + Σ(i=1:N)(i≠j) b_i"
        end

        @testset "Per-term Σ scope: distinct-index sums stay separate" begin
            # A QAdd whose .indices = [i, j] but whose terms each depend on only
            # one of them must NOT print as a global Σ_i Σ_j over the whole sum.
            @variables N
            hf = FockSpace(:f)
            af = Destroy(hf, :a)
            i = Index(hf, :i, N, hf)
            j = Index(hf, :j, N, hf)
            ai = IndexedOperator(af, i)
            aj_dag = IndexedOperator(af', j)

            expr = Σ(ai, i) + Σ(aj_dag, j)
            @test repr(expr) == "Σ(i=1:N) a_i + Σ(j=1:N) a_j'"

            # Same-index addition must group (single Σ over the sum), since both
            # terms share index i.
            ai_dag = IndexedOperator(af', i)
            same = Σ(ai, i) + Σ(ai_dag, i)
            @test repr(same) == "Σ(i=1:N) (a_i + a_i')"

            # Mixed: i-group and j-group split, i-group keeps its own pair.
            mixed = Σ(ai, i) + Σ(ai_dag, i) + Σ(aj_dag, j)
            @test repr(mixed) == "Σ(i=1:N) (a_i + a_i') + Σ(j=1:N) a_j'"

            # Diagonal-pin residual from Σa_i * Σa'_j: the bound diagonal lives
            # in scope of the surviving index only.
            prod = Σ(ai, i) * Σ(aj_dag, j)
            @test repr(prod) ==
                "N + Σ(i=1:N)Σ(j=1:N)(i≠j) a_i * a_j' + Σ(i=1:N) a_i' * a_i"
        end

        @testset "Edge cases, no crash" begin
            @test repr(_zero_qadd()) == "0"
        end

        @testset "Type inference" begin
            s = IOBuffer(sizehint = 0)
            @inferred show(s, a)
            @inferred show(s, ad)
        end
    end

    @testset "Multi-term display: real-negative coefficient uses ' - '" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        s = string(a' + a - 3 * a' * a)
        @test occursin(" - 3", s)
        @test !occursin("+ -", s)
    end

    @testset "_show_prefactor pure-imag and mixed branches" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        @variables x y

        # Pure-imag with isone(imag): expect literal "im *"
        @test occursin("im", string(im * a))

        # Pure-imag with imag == -1: expect "-im *"
        @test occursin("-im", string(-im * a))

        # Pure-imag with numeric |imag| != 1: prints "<value>im"
        @test occursin("2", string(2im * a))
        @test occursin("im", string(2im * a))

        # Pure-imag with symbolic imag: prints "<symbol>im"
        s_sym_im = string((im * x) * a)
        @test occursin("x", s_sym_im)
        @test occursin("im", s_sym_im)

        # Mixed real+imag: prints "(<re> + <im>im) *"
        s_mixed = string((x + im * y) * a)
        @test occursin("x", s_mixed)
        @test occursin("y", s_mixed)
        @test occursin("im", s_mixed)
    end

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
                (_single_qadd(_to_cnum(5), QSym[]), L"5"),
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

            σ_k = IndexedOperator(Transition(h2, :σ, 2, 2, 2), Index(h2, :k, N, NLevelSpace(:atom, 2, 1)))
            σ_l = IndexedOperator(Transition(h2, :σ, 2, 1, 2), Index(h2, :l, N, NLevelSpace(:atom, 2, 1)))
            kk = σ_k.index
            ll = σ_l.index
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
            ai = IndexedOperator(Destroy(h, :a), i)

            @test string(latexify(make_time_dependent(average(af), t))) ==
                "\\begin{equation}\n\\langle a \\rangle\\left( t \\right)\n\\end{equation}\n"
            @test string(latexify(make_time_dependent(average(adf * af), t))) ==
                "\\begin{equation}\n\\langle a^{\\dagger}a \\rangle\\left( t \\right)\n\\end{equation}\n"
            # the lifted operator carries its own Σ scope; the prefix is not doubled
            @test string(latexify(make_time_dependent(average(Σ(ai, i)), t))) ==
                "\\begin{equation}\n\\langle \\underset{i}{\\overset{N}{\\sum}}a_{i} \\rangle\\left( t \\right)\n\\end{equation}\n"
        end
    end
end
