using SecondQuantizedAlgebra
using Test
using JET
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis, basisstate

_report_text(report) = sprint(show, MIME("text/plain"), report)

function _test_allowed_only(report, allowed::Vector{String})
    reports = JET.get_reports(report)
    for rep in reports
        text = _report_text(rep)
        @test any(needle -> occursin(needle, text), allowed)
    end
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
        @test isempty(JET.get_reports(result))
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
        # per-leaf operator dispatch (`_site_compare`/`_can_commute`/`_commute_pair`),
        # so those reports are gone. Two residual `report_opt` boundaries remain, both
        # independent of the operator type:
        #
        # (a) Coefficient materialization. Constructing an operator builds its
        #     `Coeff`; the recognizer (`_rec`) and the numeric fold read a
        #     heterogeneous Symbolics expression whose `BasicSymbolic.val` is `::Any`,
        #     so the reduction (`+`/`*`/`==`/`ComplexF64`/`convert`) is dynamic. The
        #     Hermitian-symtype probe (`_avg_symtype`) reaches the same boundary from
        #     the other side: `adjoint` conjugates each coefficient via `_conj_atom`,
        #     which reads that same `::Any` symbolic tail.
        allowed_coeff_reports = [
            "SecondQuantizedAlgebra._rec(",
            "SecondQuantizedAlgebra.ComplexF64(",
            "convert(SecondQuantizedAlgebra.Coeff",
            "SecondQuantizedAlgebra.:*",
            "SecondQuantizedAlgebra.:+",
            "SecondQuantizedAlgebra.:(==)",
            "SecondQuantizedAlgebra._conj_atom(",
            ".val::Any",
        ]
        # (b) `undo_average` rebuilds a QAdd from a Symbolics `average` node, reading
        #     its `Any`-typed arguments/metadata (`_average`, `_to_qadd`) and folding
        #     the rebuilt terms through generic iteration.
        allowed_hotpath_reports = vcat(
            allowed_coeff_reports, [
                "SecondQuantizedAlgebra._average(",
                "SecondQuantizedAlgebra._to_qadd(",
                "SecondQuantizedAlgebra.iszero(",
                "SecondQuantizedAlgebra.:^",
                "Base.indexed_iterate",
                "iterate(",
            ],
        )
        # (c) Numeric conversion seals `Any → ComplexF64` (`_to_complex` /
        #     `_reduce_const` / `_fold_const` walking a Symbolics tree) and calls into
        #     QuantumOpticsBase (`to_numeric` / `expect` / `_numeric_average`). This
        #     boundary is independent of the operator collapse and was always
        #     allowlisted.
        # (d) The `_reduce_const` fallback. A variable-free symbolic constant that
        #     `_fold_const` cannot reduce (e.g. `sqrt(2)`, `exp` of a constant) is
        #     evaluated through `Symbolics.symbolic_to_float`, whose result is typed
        #     `Any`, so the final `convert`/`ComplexF64` to a concrete `ComplexF64` is
        #     dynamic. This is the same `Any → ComplexF64` boundary as the entries
        #     above, reached through `_compile_const`.
        allowed_numeric_reports = vcat(
            allowed_hotpath_reports, [
                "SecondQuantizedAlgebra._to_complex(",
                "SecondQuantizedAlgebra._reduce_const(",
                "SecondQuantizedAlgebra._fold_const(",
                "SecondQuantizedAlgebra._compile_const(",
                "SecondQuantizedAlgebra._numeric_average(",
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
                # Leaf-pair commutators with a `_commute_pair` fast path.
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
                _test_allowed_only(rep, allowed_hotpath_reports)
            end
        end

        # Bucket 3: numeric allowed
        for (name, thunk) in [
                ("to_numeric(a' * a, b)", () -> to_numeric(ad * a, b)),
                ("numeric_average(average(a), ψ)", () -> numeric_average(average(a), ψ)),
            ]
            rep = JET.@report_opt target_modules = (SecondQuantizedAlgebra,) thunk()
            @testset "$name" begin
                _test_allowed_only(rep, allowed_numeric_reports)
            end
        end
    end
end
