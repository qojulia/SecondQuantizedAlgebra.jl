using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, transition_superscript, make_time_dependent,
    expim, exponential_form, trigonometric_form

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

            @test repr(Index(h, :i, 3, h)) == "i"
            digits = NLevelSpace(:digits, 12)
            @test repr(Transition(digits, :σ, 10, 11)) == "σ₁₀₁₁"

            @variables x y
            exact_display = string(((1 // 4) * sqrt(x * y)) * a)
            @test occursin("1//4", exact_display)
            @test !occursin("0.25", exact_display)
            float_display = string((0.25 * sqrt(x * y)) * a)
            @test occursin("0.25", float_display)
            @test !occursin("1//4", float_display)
        end

        @testset "Symbolic rational coefficients stay exact after simplification" begin
            @variables ω_exact t_exact
            simplified_quartic = string(simplify((a + ad)^4 / 4))
            @test occursin("3//4", simplified_quartic)
            @test !occursin("0.75", simplified_quartic)

            exponential_cosine = string(exponential_form(cos(ω_exact * t_exact)))
            @test occursin("1//2", exponential_cosine)
            @test !occursin("0.5", exponential_cosine)
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
                (5 * commutator(a, ad), "5"),
                (a / 4, "1//4 * a"),
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

        @testset "Unitary transforms" begin
            @variables θ t
            static = Rotation(a, θ)
            @test repr(static) ==
                "UnitaryTransform(a ↦ exp(-im*θ) * a, a' ↦ exp(im*θ) * a')"
            @test repr("text/plain", static) ==
                "UnitaryTransform with 2 rules:\n" *
                "  a ↦ exp(-im*θ) * a\n" *
                "  a' ↦ exp(im*θ) * a'"

            # `:limit` is the display context used by the REPL. Large transforms
            # are summarized rather than showing an arbitrary prefix of their rules.
            levels = NLevelSpace(:levels, 3)
            σ12 = Transition(levels, :σ, 1, 2)
            large = Rotation(σ12, [1 0 0; 0 1 0; 0 0 1])
            limited = sprint(show, large; context = :limit => true)
            @test limited == "UnitaryTransform(9 rules)"
            limited_plain = sprint(show, MIME("text/plain"), large; context = :limit => true)
            @test limited_plain == "UnitaryTransform with 9 rules"
            @test !occursin("more", repr(large))

            moving = Rotation(a, θ * t, t)
            @test repr(moving) ==
                "UnitaryTransform(a ↦ exp(-im*t*θ) * a, a' ↦ exp(im*t*θ) * a'; time = t, gauge = -θ * a' * a)"
            @test repr("text/plain", moving) ==
                "Time-dependent UnitaryTransform in t with 2 rules:\n" *
                "  a ↦ exp(-im*t*θ) * a\n" *
                "  a' ↦ exp(im*t*θ) * a'\n\n" *
                "Gauge:\n" *
                "  -θ * a' * a"

            # Time metadata is meaningful even when the selected parameter is constant
            # and therefore produces a zero gauge.
            constant_timed = Rotation(a, θ, t)
            @test repr(constant_timed) ==
                "UnitaryTransform(a ↦ exp(-im*θ) * a, a' ↦ exp(im*θ) * a'; time = t)"
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

        @testset "Heterogeneous operator averages sort without error" begin
            h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
            @qnumbers b::Destroy(h2, 1)
            σ = Transition(h2, :σ, 2, 2, 2)

            rendered = repr(average(b' * b) + average(σ) + average(b))
            @test all(
                occursin(term, rendered) for
                    term in ("⟨b⟩", "⟨σ₂₂⟩", "⟨b' * b⟩")
            )
        end

        @testset "Lifted (time-dependent) averages" begin
            @variables t N
            hsite = FockSpace(:site)
            asite = Destroy(hsite, :a)
            i = Index(hsite, :i, N, hsite)
            j = Index(hsite, :j, N, hsite)
            ai = IndexedOperator(asite, i)
            σ = Transition(NLevelSpace(:atom, 3, 1), :σ, 1, 2)

            @test repr(make_time_dependent(average(af), t)) == "⟨a⟩(t)"
            @test repr(make_time_dependent(average(adf * af), t)) == "⟨a' * a⟩(t)"
            @test repr(make_time_dependent(average(σ), t)) == "⟨σ₁₂⟩(t)"
            # Lifting descends into the sum body: each per-site moment becomes
            # time-dependent and the Σ stays outside, matching the non-lifted form
            # `Σ(i) ⟨a_i⟩` rather than collapsing to one collective ⟨Σ a_i⟩(t).
            @test repr(make_time_dependent(average(Σ(ai, i)), t)) == "Σ(i=1:N) ⟨a_i⟩(t)"
            @test repr(make_time_dependent(average(Σ(ai, i, [j])), t)) ==
                "Σ(i=1:N)(i≠j) ⟨a_i⟩(t)"
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

            prod = Σ(ai, i) * Σ(aj_dag, j)
            @test repr(prod) ==
                "N + Σ(i=1:N)Σ(j=1:N)(i≠j) a_i * a_j' + Σ(i=1:N) a_i' * a_i"
        end

        @testset "Edge cases, no crash" begin
            @test repr(a - a) == "0"
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

        raw_negative = string(a + (-big"1e100") * ad)
        @test occursin(" - ", raw_negative)
        @test !occursin("+ -", raw_negative)
    end

    @testset "explicit phase display" begin
        @variables θ g
        phase_cos = (expim(θ) + expim(-θ)) * a
        phase_sin = -im * (expim(θ) - expim(-θ)) * a
        grouped = g * (expim(θ) + expim(-θ)) * a
        @test count("exp(", string(phase_cos)) == 2
        @test count("exp(", string(phase_sin)) == 2
        @test count("exp(", string(grouped)) == 2
        @test !occursin("cos", string(phase_cos))
        @test !occursin("sin", string(phase_sin))
        @test !occursin("cos", string(latexify(phase_cos)))
        @test string(trigonometric_form(phase_cos)) == "2cos(θ) * a"
        @test occursin("cos", string(latexify(trigonometric_form(phase_cos))))
    end

    @testset "show_prefactor pure-imag and mixed branches" begin
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

    @testset "show_prefactor braces composite parts" begin
        # A composite part must be parenthesized, or the `im` suffix binds to its last
        # factor only: `x + yim * a` reads as `x + y*im*a`. Anchored patterns, so the
        # summand order Symbolics picks does not matter.
        h = FockSpace(:f)
        a = Destroy(h, :a)
        @variables x y g

        @test occursin(r"^\(.+\)im \* a$", string((im * (x + y)) * a))
        @test occursin(r"^\(.+\)im \* a$", string((im * (x / y)) * a))
        # The `im` suffix is juxtaposition, so every call has to be braced, not just the
        # loose heads: `x^2im` reads as `x^(2im)`, `x*yim` merges into an identifier, and
        # `sqrt(x)im` is a syntax error.
        @test string((im * x^2) * a) == "(x^2)im * a"
        @test occursin(r"^\(.+\)im \* a$", string((im * x * y) * a))
        @test occursin(r"^\(.+\)im \* a$", string((im * sqrt(x)) * a))
        @test occursin(r"^\(g \+ \(x\^2\)im\) \* a$", string((g + im * x^2) * a))
        @test occursin(r"^\(\(.+\) \+ gim\) \* a$", string(((x + y) + im * g) * a))
        @test occursin(r"^\(g \+ \(.+\)im\) \* a$", string((g + im * (x + y)) * a))
        @test occursin(r"^\(\(.+\) \+ \(.+\)im\) \* a$", string(((x + y) + im * (x - y)) * a))
        @test string((im * x) * a) == "xim * a"
        @test occursin(r"^\(.+\) \* a$", string((x + y) * a))
        @test occursin(r"^\(.+\) \* a$", string((x / y) * a))
        @test string((x + y) * commutator(a, ad)) == string(x + y)
    end
end
