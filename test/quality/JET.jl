using SecondQuantizedAlgebra
using Test
using JET
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis, basisstate

report_text(report) = sprint(show, MIME("text/plain"), report)

# Collected rather than asserted one by one so the failure can name the offender: a bare
# `@test any(occursin, allowed)` prints only the predicate.
matches_allowed(text, needle::AbstractString) = occursin(needle, text)
matches_allowed(text, needle::Regex) = occursin(needle, text)

unmatched(report, allowed::AbstractVector) = filter(
    text -> !any(needle -> matches_allowed(text, needle), allowed),
    map(report_text, JET.get_reports(report)),
)

function test_allowed_only(report, allowed::AbstractVector)
    left = unmatched(report, allowed)
    isempty(left) || @info "unmatched JET reports" left
    @test isempty(left)
    return
end

@static if isempty(VERSION.prerelease)
    @testset "JET report_package (correctness)" begin
        # report_package finds actual type errors (method-not-found, etc.).
        # `ignore_missing_comparison = true
        # silences `Missing` union-split warnings that bubble up from Symbolics.
        # report_opt is checked separately (below) on selected hot paths; collapsing
        # the operators to a concrete `Vector{Op}` removed the per-leaf dynamic
        # dispatch the abstract hierarchy used to force, leaving only the `Coeff`
        # materialization and numeric/average conversion boundaries.
        result = JET.report_package(
            SecondQuantizedAlgebra;
            target_modules = (SecondQuantizedAlgebra,),
            ignore_missing_comparison = true,
        )
        # `Symbolics` defines `Num` only for the `SymReal` variant, so every
        # `Num(x::BasicSymbolic)` reads as a possible `MethodError` on the `TreeReal` arm.
        # Nothing here produces one: `unwrap(r)`, `unwrap(cosh(r))` and `unwrap(ω*t)` are all
        # `SymReal`. Narrowing the callees only relocates the report, since the widening
        # starts at `recognize(::BasicSymbolic)`, and `Symbolics.wrap` infers `Any`. Spelled
        # out in full so a SymbolicUtils rename re-fires the gate — both renderings, since
        # JET qualifies the type parameter only when the running session lacks the binding.
        test_allowed_only(
            result,
            [
                "Num(::SymbolicUtils.BasicSymbolicImpl.var\"typeof(BasicSymbolicImpl)\"" *
                    "{SymbolicUtils.TreeReal})",
                "Num(::SymbolicUtils.BasicSymbolicImpl.var\"typeof(BasicSymbolicImpl)\"" *
                    "{TreeReal})",
            ],
        )
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
        # Three buckets:
        #
        #  1. Strict (leaf-level operations). Fully inferred, zero dispatch reports.
        #     With the concrete `Vector{Op}` storage the per-operator hooks and
        #     leaf-level `to_numeric` infer statically, so these stay clean.
        #
        #  2. Hot-path allowed. QAdd-level pipelines whose only residual dispatch is
        #     the `Coeff` materialization boundary, plus `undo_average`'s rebuild from
        #     a Symbolics `average` node (see the allowlists below).
        #
        #  3. Numeric allowed. `to_numeric` / `numeric_average` on QAdd inputs, which
        #     additionally cross the `Any → ComplexF64` numeric-conversion boundary
        #     and call into QuantumOpticsBase.
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = Create(hf, :a)
        hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)

        b = FockBasis(7)
        bn = NLevelBasis(3)
        bs = SpinBasis(1 // 2)
        ψ = basisstate(b, 1)

        # Bucket 1: strict (zero dispatch reports)
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
            ]
            rep = JET.@report_opt target_modules = (SecondQuantizedAlgebra,) thunk()
            @testset "$name" begin
                @test isempty(JET.get_reports(rep))
            end
        end

        # Collapsing the operators to a concrete `Vector{Op}` removed the former
        # per-leaf operator dispatch (`site_compare`/`can_commute`/`commute_pair`),
        # so those reports are gone. Two residual `report_opt` boundaries remain, both
        # independent of the operator type:
        #
        # (a) Coefficient materialization. Constructing an operator builds its
        #     `Coeff`; the recognizer (`recognize`) and the numeric fold read a
        #     heterogeneous Symbolics expression whose `BasicSymbolic.val` is `::Any`,
        #     so the reduction (`+`/`*`/`==`/`ComplexF64`/`convert`) is dynamic. The
        #     Hermitian-symtype probe (`avg_symtype`) reaches the same boundary from
        #     the other side: `adjoint` conjugates each coefficient via `conj_atom`,
        #     which reads that same `::Any` symbolic tail.
        # The extra entries below are deliberately type-shaped: a new function name alone
        # must not silence every future report from that function.
        allowed_coeff_reports = Any[
            "SecondQuantizedAlgebra.recognize(",
            r"SecondQuantizedAlgebra\.raw_complex\(%\d+::Any, %\d+::Any\)",
            r"Core\.kwcall\(%\d+::@NamedTuple\{(?:normalize::Bool, )?real_slot::Bool\},",
            r"SecondQuantizedAlgebra\.to_cnum\(%\d+::Number\)::Any",
            r"SecondQuantizedAlgebra\.to_cnum\(%\d+::Any\)::Any",
            r"SecondQuantizedAlgebra\.to_cnum\(%\d+::Rational\)::(?:SecondQuantizedAlgebra\.)?Coeff",
            r"SecondQuantizedAlgebra\.to_cnum\(%\d+::Complex\)::(?:SecondQuantizedAlgebra\.)?Coeff",
            r"SecondQuantizedAlgebra\.normalize_phase\(%\d+::Any\)::Any",
            r"SecondQuantizedAlgebra\.Complex\(%\d+::Any, %\d+::Any\)::Complex",
            r"SecondQuantizedAlgebra\.imag\(%\d+::Number\)::Real",
            r"\(%\d+::Any SecondQuantizedAlgebra\.:/ %\d+::Any\)::Any",
            r"\(%\d+::SymbolicUtils\.BasicSymbolicImpl.*SecondQuantizedAlgebra\.:/ %\d+::Any\)::Any",
            r"\(%\d+::Any SecondQuantizedAlgebra\.:- %\d+::Any\)::Any",
            r"get_variables\(%\d+::Real\)::Any",
            r"SecondQuantizedAlgebra\.any\(%\d+::.*::Any\)::Any",
            r"SecondQuantizedAlgebra\.:// %\d+::Integer\)::Rational",
            r"SecondQuantizedAlgebra\.collect_trig!\(%\d+::Vector\{SymbolicUtils\.BasicSymbolicImpl.*?, %\d+::Any\)::Any",
            r"SecondQuantizedAlgebra\.append!\(%\d+::Vector\{SymbolicUtils\.BasicSymbolicImpl.*?, %\d+::Any\)::Any",
            r"SecondQuantizedAlgebra\.has_symbolic_trig\(%\d+::Any\)::Bool",
            "SecondQuantizedAlgebra.ComplexF64(",
            "convert(SecondQuantizedAlgebra.Coeff",
            "SecondQuantizedAlgebra.:*",
            "SecondQuantizedAlgebra.:+",
            "SecondQuantizedAlgebra.:(==)",
            "SecondQuantizedAlgebra.conj_atom(",
            ".val::Any",
        ]
        # (a2) The same `::Any` boundary further downstream. `BasicSymbolic` is a UnionAll, so
        #     `operation`/`arguments` on one return `Any` and every call after that inherits
        #     it. Two cold walks reach it: `qadjoint(::Num)` unwraps to the field before its
        #     phase and sign arms run, and the trig discovery in `reduce.jl` iterates
        #     `Monomial.syms`, whose eltype is that same UnionAll. Neither is fixable here:
        #     the type is already lost at the caller, and the eltype is SymbolicUtils' choice.
        allowed_symbolic_walk_reports = Any[
            "BasicSymbolicImpl)\"{T} where T)",
            "operation(",   # the two SymbolicUtils accessors that start the widening
            "arguments(",
            "SecondQuantizedAlgebra.strip_conj(",
            "SecondQuantizedAlgebra.only(",
            "SecondQuantizedAlgebra.expim_symbolic(",
            "SecondQuantizedAlgebra.sign(",
            "SecondQuantizedAlgebra.exp(",
            "SecondQuantizedAlgebra.:!(",
            "SecondQuantizedAlgebra.:<",
            "SecondQuantizedAlgebra.:-(",
            "SecondQuantizedAlgebra.unary_arg(",
            "SecondQuantizedAlgebra.find_partner(",
            "SecondQuantizedAlgebra.TrigHead(",
            "SecondQuantizedAlgebra.length(",
            "SecondQuantizedAlgebra.isequal(",
            "Dict{SecondQuantizedAlgebra.TrigKey",
            "Dict{Symbolics.Num, Symbolics.Num}",
            "Dict{Num, Num}",
            "unwrap(",
            ")[1]::Any",
            "SecondQuantizedAlgebra.Num(",
            "SecondQuantizedAlgebra.ParamRelation(",
            r"SecondQuantizedAlgebra\.simplify_raw_component\(%\d+::Any\)::Any",
        ]
        # (b) `undo_average` rebuilds a QAdd from a Symbolics `average` node, reading
        #     its `Any`-typed arguments/metadata (`average`, `to_qadd`) and folding
        #     the rebuilt terms through generic iteration.
        allowed_hotpath_reports = vcat(
            allowed_coeff_reports, allowed_symbolic_walk_reports, [
                "SecondQuantizedAlgebra.average(",
                "SecondQuantizedAlgebra.to_qadd(",
                "SecondQuantizedAlgebra.iszero(",
                "SecondQuantizedAlgebra.:^",
                "Base.indexed_iterate",
                "iterate(",
            ],
        )
        # (c) Numeric conversion seals `Any → ComplexF64` (`to_complex` /
        #     `reduce_const` / `fold_const` walking a Symbolics tree) and calls into
        #     QuantumOpticsBase (`to_numeric` / `expect` / `numeric_average`). This
        #     boundary is independent of the operator collapse and was always
        #     allowlisted.
        # (d) The `reduce_const` fallback. A variable-free symbolic constant that
        #     `fold_const` cannot reduce (e.g. `sqrt(2)`, `exp` of a constant) is
        #     evaluated through `Symbolics.symbolic_to_float`, whose result is typed
        #     `Any`, so the final `convert`/`ComplexF64` to a concrete `ComplexF64` is
        #     dynamic. This is the same `Any → ComplexF64` boundary as the entries
        #     above, reached through `compile_const`.
        allowed_numeric_reports = vcat(
            allowed_hotpath_reports, [
                "SecondQuantizedAlgebra.to_complex(",
                "SecondQuantizedAlgebra.reduce_const(",
                "SecondQuantizedAlgebra.fold_const(",
                "SecondQuantizedAlgebra.compile_const(",
                "SecondQuantizedAlgebra.numeric_average(",
                "SecondQuantizedAlgebra.numeric_average_impl(",
                "SecondQuantizedAlgebra.Complex(",
                "SecondQuantizedAlgebra.convert(",
                "convert(SecondQuantizedAlgebra.ComplexF64",
                "to_numeric(",
                "expect(",
                "string(",
            ],
        )

        # Bucket 2: hot-path allowed
        for (name, thunk) in [
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
            ]
            rep = JET.@report_opt target_modules = (SecondQuantizedAlgebra,) thunk()
            @testset "$name" begin
                test_allowed_only(rep, allowed_hotpath_reports)
            end
        end

        # Bucket 3: numeric allowed
        for (name, thunk) in [
                ("to_numeric(a' * a, b)", () -> to_numeric(ad * a, b)),
                ("numeric_average(average(a), ψ)", () -> numeric_average(average(a), ψ)),
            ]
            rep = JET.@report_opt target_modules = (SecondQuantizedAlgebra,) thunk()
            @testset "$name" begin
                test_allowed_only(rep, allowed_numeric_reports)
            end
        end
    end
end
