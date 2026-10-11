using SecondQuantizedAlgebra
using Test
using JET
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis, basisstate

report_text(report) = sprint(show, MIME("text/plain"), report)

# Prints the offending reports on failure: a bare `@test isempty(...)` shows only the count.
function test_no_reports(report)
    reports = map(report_text, JET.get_reports(report))
    isempty(reports) || @info "unexpected JET reports" reports
    @test isempty(reports)
    return
end

@static if isempty(VERSION.prerelease)
    @testset "JET report_package (correctness)" begin
        # report_package finds actual type errors (method-not-found, etc.).
        # `ignore_missing_comparison = true
        # silences `Missing` union-split warnings that bubble up from Symbolics.
        # report_opt is checked separately (below) on selected hot paths.
        result = JET.report_package(
            SecondQuantizedAlgebra;
            target_modules = (SecondQuantizedAlgebra,),
            ignore_missing_comparison = true,
        )
        test_no_reports(result)
    end

    @testset "JET report_call on entry points (no errors)" begin
        # @report_call catches actual MethodError/UndefVar issues at runtime.
        # Entries are grouped by what they exercise: per-family products and
        # commutators, multi-body and mixed-space products, and the canonicalisation /
        # average / substitute / expand_completeness pipelines.
        #
        # `a^n` is skipped because JET descends through Symbolics' `Num^Integer`
        # and surfaces hundreds of macro-generated UndefVarError reports from
        # SymbolicUtils' BasicSymbolic internals — outside our control.
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = Create(hf, :a)
        hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)

        hjc = FockSpace(:f) ⊗ NLevelSpace(:atom, 2)
        ajc = Destroy(hjc, :a, 1); σjc = Transition(hjc, :σ, 1, 2, 2); σjc_p = Transition(hjc, :σ, 2, 1, 2)
        @variables θ::Real ω::Real t::Real
        U = Rotation(a, θ)
        Ut = Rotation(a, ω * t, t)

        @variables u::Number v::Number
        bogo_matrix = [u v; conj(v) conj(u)]
        B = Bogoliubov(a, bogo_matrix)

        for (name, expr) in [
                # Per-family binary products
                ("a' * a", () -> ad * a),
                ("σ12 * σ21", () -> σ12 * σ21),
                ("px * py", () -> px * py),
                ("sx * sy", () -> sx * sy),
                ("x * p", () -> xx * pp),
                # Per-family commutators
                ("commutator(a, a')", () -> commutator(a, ad)),
                ("commutator(σ12, σ21)", () -> commutator(σ12, σ21)),
                ("commutator(px, py)", () -> commutator(px, py)),
                ("commutator(sx, sy)", () -> commutator(sx, sy)),
                ("commutator(x, p)", () -> commutator(xx, pp)),
                # QAdd-QAdd commutator and multi-body products
                ("commutator(a*a', a'*a)", () -> commutator(a * ad, ad * a)),
                ("a * a' * a", () -> a * ad * a),
                ("(a + a') * (a - a')", () -> (a + ad) * (a - ad)),
                # Mixed ProductSpace (Jaynes-Cummings setup)
                ("a' * σ12 (Fock ⊗ NLevel)", () -> ajc * σjc),
                ("a*σ + a*σ' (Fock ⊗ NLevel)", () -> ajc * σjc + ajc * σjc_p),
                # Canonicalisation pipelines
                ("normal_order(a'*a*a')", () -> normal_order(ad * a * ad)),
                ("normal_order(σ21*σ12*σ21)", () -> normal_order(σ21 * σ12 * σ21)),
                ("normal_order(px*py)", () -> normal_order(px * py)),
                ("simplify(a*a + a'*a')", () -> simplify(a * a + ad * ad)),
                # Substitute, expand_completeness, average round-trip
                ("substitute(a, Dict(a=>a'))", () -> substitute(a, Dict(a => ad))),
                ("expand_completeness(σ12*σ21)", () -> expand_completeness(σ12 * σ21)),
                ("average(a*a')", () -> average(a * ad)),
                ("undo_average(average(a*a'))", () -> undo_average(average(a * ad))),
                # MutableArithmetics additive reductions
                ("sum([a'*a, a*a', a'*a])", () -> sum([ad * a, a * ad, ad * a])),
                ("reduce(+, [a'*a, a*a', a'*a])", () -> reduce(+, [ad * a, a * ad, ad * a])),
                # Exact unitary-transform public entry points.
                ("Rotation(a, θ)", () -> Rotation(a, θ)),
                ("conjugate(a, U)", () -> conjugate(a, U)),
                ("transform(a'*a, Ut)", () -> transform(ad * a, Ut)),
                ("inv(U)", () -> inv(U)),
                ("U * Ut", () -> U * Ut),
                (
                    "exponential_form(cos(ω*t)*a)",
                    () -> SecondQuantizedAlgebra.exponential_form(cos(ω * t) * a),
                ),
                (
                    "trigonometric_form(expim(ω*t)*a)",
                    () -> SecondQuantizedAlgebra.trigonometric_form(
                        SecondQuantizedAlgebra.expim(ω * t) * a,
                    ),
                ),
                (
                    "phase_terms(cos(ω*t))",
                    () -> SecondQuantizedAlgebra.phase_terms(
                        SecondQuantizedAlgebra.exponential_form(cos(ω * t)),
                    ),
                ),
                (
                    "expim(ω*t) * expim(θ)",
                    () -> SecondQuantizedAlgebra.expim(ω * t) *
                        SecondQuantizedAlgebra.expim(θ),
                ),
                (
                    "substitute(expim(ω*t), Dict(ω=>θ))",
                    () -> substitute(
                        SecondQuantizedAlgebra.expim(ω * t), Dict(ω => θ),
                    ),
                ),
                ("real(expim(θ))", () -> real(SecondQuantizedAlgebra.expim(θ))),
                ("abs2(expim(θ))", () -> abs2(SecondQuantizedAlgebra.expim(θ))),
                # Bogoliubov transforms and the affine substitution that resolves them.
                ("Bogoliubov(a, integer matrix)", () -> Bogoliubov(a, [1 0; 0 1])),
                (
                    "Bogoliubov(a, complex matrix)",
                    () -> Bogoliubov(a, ComplexF64[im 0; 0 -im]),
                ),
                ("Bogoliubov(a, symbolic matrix)", () -> Bogoliubov(a, bogo_matrix)),
                (
                    "substitute(Bogoliubov, Dict(u=>5//3, v=>4//3))",
                    () -> substitute(B, Dict(u => 5 // 3, v => 4 // 3)),
                ),
            ]
            rep = JET.@report_call target_modules = (SecondQuantizedAlgebra,) ignore_missing_comparison = true expr()
            @testset "$name" begin
                @test isempty(JET.get_reports(rep))
            end
        end
    end

    @testset "JET report_opt on selected hot paths" begin
        # Every entry is strict: zero runtime dispatch inside the package. The symbolic
        # walks stay inferred because coefficient atoms are typed `RawExpression` and a
        # constant's `::Any` payload is read only through the narrowing helpers in `cnum.jl`.
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = Create(hf, :a)
        hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)

        b = FockBasis(7)
        bn = NLevelBasis(3)
        bs = SpinBasis(1 // 2)
        ψ = basisstate(b, 1)

        for (name, thunk) in [
                # Per-family to_numeric on a leaf
                ("to_numeric(a, b)", () -> to_numeric(a, b)),
                ("to_numeric(σ12, bn)", () -> to_numeric(σ12, bn)),
                ("to_numeric(px, bs)", () -> to_numeric(px, bs)),
                ("to_numeric(sx, bs)", () -> to_numeric(sx, bs)),
                ("to_numeric(x, b)", () -> to_numeric(xx, b)),
                # Per-family average on a leaf
                ("average(a)", () -> average(a)),
                ("average(σ12)", () -> average(σ12)),
                ("average(px)", () -> average(px)),
                ("average(sx)", () -> average(sx)),
                ("average(x)", () -> average(xx)),
                # Leaf-level predicates and accessors
                ("acts_on(a)", () -> acts_on(a)),
                ("acts_on(σ12)", () -> acts_on(σ12)),
                ("is_average(a)", () -> is_average(a)),
                ("is_average(average(a))", () -> is_average(average(a))),
                ("has_sum_metadata(average(a))", () -> SecondQuantizedAlgebra.has_sum_metadata(average(a))),
                # Leaf-pair commutators with a `commute_pair` fast path.
                ("commutator(a, a')", () -> commutator(a, ad)),
                ("commutator(x, p)", () -> commutator(xx, pp)),
                # Commutators that fall through to `a*b - b*a`
                ("commutator(σ12, σ21)", () -> commutator(σ12, σ21)),
                ("commutator(px, py)", () -> commutator(px, py)),
                ("commutator(sx, sy)", () -> commutator(sx, sy)),
                # QAdd-QAdd commutator and multi-body products
                ("commutator(a*a', a'*a)", () -> commutator(a * ad, ad * a)),
                ("a*a'*a", () -> a * ad * a),
                ("(a + a')*(a - a')", () -> (a + ad) * (a - ad)),
                # Canonicalisation pipelines on QAdd
                ("normal_order(a)", () -> normal_order(a)),
                ("simplify(a*a + a'*a')", () -> simplify(a * a + ad * ad)),
                ("expand_completeness(σ12*σ21)", () -> expand_completeness(σ12 * σ21)),
                # Average/undo_average round trip
                ("average(a*a')", () -> average(a * ad)),
                ("undo_average(average(a))", () -> undo_average(average(a))),
                # Numeric conversion of a QAdd and of an average
                ("to_numeric(a' * a, b)", () -> to_numeric(ad * a, b)),
                ("numeric_average(average(a), ψ)", () -> numeric_average(average(a), ψ)),
            ]
            rep = JET.@report_opt target_modules = (SecondQuantizedAlgebra,) thunk()
            @testset "$name" begin
                test_no_reports(rep)
            end
        end
    end
end
