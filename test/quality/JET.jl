using SecondQuantizedAlgebra
using Test
using JET
using QuantumOpticsBase: FockBasis, NLevelBasis, SpinBasis, LazyKet, basisstate, coherentstate, dense, destroy, fockstate, nlevelstate, tensor
using Latexify: latexify
using QuantumToolbox: QuantumToolbox
import SciMLOperators as SO
using SecondQuantizedAlgebra: constraint_pairs, exponential_form, expim, get_sum_body, get_sum_indices, get_sum_non_equal, indexed_sum, order_key, phase_terms, qadd_order_key, rename, set_acts_on, term_order_key, to_cnum, to_num, transition_superscript, trigonometric_form

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
        result = JET.report_package(
            SecondQuantizedAlgebra;
            target_modules = (SecondQuantizedAlgebra,),
        )
        test_no_reports(result)
    end

    @testset "JET report_package ($ext)" for ext in (
            :SecondQuantizedAlgebraQuantumOpticsBaseExt,
            :SecondQuantizedAlgebraQuantumToolboxExt,
        )
        m = Base.get_extension(SecondQuantizedAlgebra, ext)
        @test m isa Module
        test_no_reports(JET.report_package(m; target_modules = (m,)))
    end

    @testset "JET report_call on entry points (no errors)" begin
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
        A = a + ad
        F = cos(A)

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
                ("a^3", () -> a^3),
                ("a'^2 * a", () -> ad^2 * a),
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
                # Formal operator-expression cold path.
                ("cos(a + a')", () -> cos(A)),
                ("sin(a + a')", () -> sin(A)),
                ("expim(a + a')", () -> expim(A)),
                ("substitute(cos(A), Dict(a=>a'))", () -> substitute(F, Dict(a => ad))),
                ("taylor(cos(A), 0:4)", () -> taylor(F, 0:4)),
                # Exact unitary-transform public entry points.
                ("Rotation(a, θ)", () -> Rotation(a, θ)),
                ("conjugate(a, U)", () -> conjugate(a, U)),
                ("conjugate(cos(A), U)", () -> conjugate(F, U)),
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
            rep = JET.@report_call target_modules = (SecondQuantizedAlgebra,) expr()
            @testset "$name" begin
                test_no_reports(rep)
            end
        end
    end

    @testset "JET report_opt (no runtime dispatch)" begin
        # Every entry is strict: zero runtime dispatch inside the package and both backend
        # extensions.
        targets = (
            SecondQuantizedAlgebra,
            Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumOpticsBaseExt),
            Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumToolboxExt),
        )
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = Create(hf, :a)
        hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)

        b = FockBasis(7)
        bn = NLevelBasis(3)
        bs = SpinBasis(1 // 2)
        ψ = basisstate(b, 1)

        hjc = FockSpace(:c) ⊗ NLevelSpace(:atom, 2)
        ajc = Destroy(hjc, :a, 1); σjc = Transition(hjc, :σ, 1, 2, 2)
        Hjc = ajc' * σjc + ajc * σjc'
        ψjc = basisstate(FockBasis(3) ⊗ NLevelBasis(2), 1)
        ψq = QuantumToolbox.fock(8, 1)
        ψqjc = QuantumToolbox.tensor(QuantumToolbox.fock(4, 0), QuantumToolbox.basis(2, 1))

        @variables g::Real Δ::Real θ::Real α::Number
        Hp = Δ * ad * a + g * (a + ad)
        U = Rotation(a, θ)
        D = Displace(a, α)

        hi = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        ai = Destroy(hi, :a, 1)
        σi(k, l, m) = IndexedOperator(Transition(hi, :σ, k, l, 2), m)
        i = Index(hi, :i, 3, NLevelSpace(:atom, 2)); j = Index(hi, :j, 3, NLevelSpace(:atom, 2))
        Hi = -Δ * ai' * ai + Σ(g * (ai' * σi(1, 2, i) + ai * σi(2, 1, i)), i)
        bi = FockBasis(3) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2)
        sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3, 4])
        Hin = substitute(Hi, Dict(g => 1.0, Δ => 0.5))
        custom = 2.0 * destroy(b)
        @variables t::Real ω::Real r::Real ϕ::Real u::Number v::Number
        H0 = ω * ajc' * ajc
        hc = CollectiveNLevelSpace(:atom, 3); S12 = CollectiveTransition(hc, :S, 1, 2)
        B = Bogoliubov(a, [u v; conj(v) conj(u)])
        @variables κ::Real
        Lq = to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => t -> 1 + im * t)).data
        Lp = to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => (p, t) -> p[1] * cos(t))).data
        Ls = to_numeric(ad * a + a, 8; op_type = identity).data
        c1 = to_numeric(sqrt(κ) * a, 8; time_parameter = Dict(κ => (p, t) -> p[1]))
        Ll = QuantumToolbox.liouvillian(to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => t -> cos(t))), [c1]).data
        vq = rand(ComplexF64, 9); wq = similar(vq)
        vl = rand(ComplexF64, 81); wl = similar(vl)
        ψl = LazyKet(FockBasis(3) ⊗ NLevelBasis(2), (fockstate(FockBasis(3), 1), nlevelstate(NLevelBasis(2), 1)))

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
                ("a^3", () -> a^3),
                ("(a + a')*(a - a')", () -> (a + ad) * (a - ad)),
                # Canonicalisation pipelines on QAdd
                ("normal_order(a)", () -> normal_order(a)),
                ("simplify(a*a + a'*a')", () -> simplify(a * a + ad * ad)),
                ("expand_completeness(σ12*σ21)", () -> expand_completeness(σ12 * σ21)),
                # Average/undo_average round trip
                ("average(a*a')", () -> average(a * ad)),
                ("undo_average(average(a))", () -> undo_average(average(a))),
                # Substitution, unitary transforms, Weyl ordering
                ("substitute(a, Dict(a => a'))", () -> substitute(a, Dict(a => ad))),
                ("substitute(Hp, Dict(g => 1.0, Δ => 2.0))", () -> substitute(Hp, Dict(g => 1.0, Δ => 2.0))),
                ("substitute(Hp, Dict(a => a + 1))", () -> substitute(Hp, Dict(a => a + 1))),
                ("substitute(U, Dict(θ => 0.5))", () -> substitute(U, Dict(θ => 0.5))),
                ("transform(a'*a, Displace(a, α))", () -> transform(ad * a, D)),
                ("conjugate(a, Rotation(a, θ))", () -> conjugate(a, U)),
                ("normal_to_symmetric(a'*a)", () -> normal_to_symmetric(ad * a)),
                ("symmetric_to_normal(a'*a)", () -> symmetric_to_normal(ad * a)),
                ("transform(a'*a, Squeeze(a, r, ϕ))", () -> transform(ad * a, Squeeze(a, r, ϕ))),
                ("transform(a'*a, Bogoliubov)", () -> transform(ad * a, B)),
                ("RotatingFrame(H0, t)", () -> RotatingFrame(H0, t)),
                ("DisplacementFrame(a, ω a'a)", () -> DisplacementFrame(a, ω * ad * a)),
                ("gauge_term(RotatingFrame(H0, t))", () -> gauge_term(RotatingFrame(H0, t))),
                (
                    "trigonometric_form(expim(ω*t)*a)",
                    () -> SecondQuantizedAlgebra.trigonometric_form(SecondQuantizedAlgebra.expim(ω * t) * a),
                ),
                ("anticommutator(a, a')", () -> anticommutator(a, ad)),
                ("find_operators(hjc, 2)", () -> find_operators(hjc, 2)),
                ("S12 * S12'", () -> S12 * S12'),
                # Printing
                ("repr(Hp)", () -> repr(Hp)),
                ("latexify(Hp)", () -> latexify(Hp)),
                ("latexify(average(indexed H))", () -> latexify(average(Hi))),
                # Indexed sums
                ("Σ(σ21_i) * Σ(σ12_j)", () -> Σ(σi(2, 1, i), i) * Σ(σi(1, 2, j), j)),
                ("simplify(indexed H)", () -> simplify(Hi)),
                ("average(indexed H)", () -> average(Hi)),
                ("change_index(indexed H, i, j)", () -> change_index(Hi, i, j)),
                ("assume_distinct_index(σ21_i σ12_j)", () -> assume_distinct_index(σi(2, 1, i) * σi(1, 2, j), [(i, j)])),
                ("make_time_dependent(average(indexed H), t)", () -> make_time_dependent(average(Hi), t)),
                # QuantumOptics numeric conversion and expectation values
                ("to_numeric(a' * a, b)", () -> to_numeric(ad * a, b)),
                ("to_numeric(a' * a, b; op_type = identity)", () -> to_numeric(ad * a, b; op_type = identity)),
                ("to_numeric(Hp, b; parameter)", () -> to_numeric(Hp, b; parameter = Dict(g => 1.0, Δ => 0.5))),
                ("to_numeric(2.0, b)", () -> to_numeric(2.0, b)),
                ("to_numeric([a, a'], b)", () -> to_numeric([a, ad], b)),
                (
                    "to_numeric(a' * a, hf, 7; QuantumOpticsBackend)",
                    () -> to_numeric(ad * a, hf, 7; backend = QuantumOpticsBackend()),
                ),
                ("to_numeric(Hjc, bjc)", () -> to_numeric(Hjc, FockBasis(3) ⊗ NLevelBasis(2))),
                ("to_numeric(indexed H, b, sites)", () -> to_numeric(Hin, bi, sites)),
                ("to_numeric(a' * a, b; operators)", () -> to_numeric(ad * a, b; operators = Dict(a => custom))),
                (
                    "to_numeric(α * a + h.c., b; time_parameter)",
                    () -> to_numeric(α * a + conj(α) * ad, b; time_parameter = Dict(α => t -> 1 + im * t)),
                ),
                ("numeric_average(average(a), ψ)", () -> numeric_average(average(a), ψ)),
                ("numeric_average(a' * a, [ψ, ψ])", () -> numeric_average(ad * a, [ψ, ψ])),
                ("numeric_average(Hjc, ψjc)", () -> numeric_average(Hjc, ψjc)),
                ("numeric_average(Hjc, LazyKet)", () -> numeric_average(Hjc, ψl)),
                # QuantumToolbox numeric conversion and expectation values
                ("to_numeric(a' * a, 8)", () -> to_numeric(ad * a, 8)),
                ("to_numeric(Hjc, (4, 2))", () -> to_numeric(Hjc, (4, 2))),
                ("solver apply: time-dependent QobjEvo", () -> Lq(wq, vq, nothing, nothing, 0.3)),
                ("solver apply: p-aware QobjEvo", () -> Lp(wq, vq, nothing, [0.5], 0.3)),
                ("solver apply: static lazy QobjEvo", () -> Ls(wq, vq, nothing, nothing, 0.3)),
                ("solver apply: Liouvillian with time-dependent collapse", () -> Ll(wl, vl, nothing, [0.5], 0.3)),
                (
                    "to_numeric(Hjc, hjc, (3, 2); QuantumToolboxBackend)",
                    () -> to_numeric(Hjc, hjc, (3, 2); backend = QuantumToolboxBackend()),
                ),
                ("numeric_average(a' * a, ψq)", () -> numeric_average(ad * a, ψq)),
                ("numeric_average(average(a), ψq)", () -> numeric_average(average(a), ψq)),
                ("numeric_average(Hjc, ψqjc)", () -> numeric_average(Hjc, ψqjc)),
                ("numeric_average(a' * a, [ψq, ψq])", () -> numeric_average(ad * a, [ψq, ψq])),
            ]
            rep = JET.@report_opt target_modules = targets thunk()
            @testset "$name" begin
                test_no_reports(rep)
            end
        end
    end

    @testset "JET report_opt (feature coverage)" begin
        targets = (
            SecondQuantizedAlgebra,
            Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumOpticsBaseExt),
            Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumToolboxExt),
        )
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = a'
        hn = NLevelSpace(:atom, 3); σ12 = Transition(hn, :σ, 1, 2); σ11 = Transition(hn, :σ, 1, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2); pz = Pauli(hp, :σ, 3)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sz = Spin(hs, :S, 3)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)
        hjc = FockSpace(:c) ⊗ NLevelSpace(:atom, 2)
        ajc = Destroy(hjc, :a, 1); σjc = Transition(hjc, :σ, 1, 2, 2)
        hc = CollectiveNLevelSpace(:atom, 3); S12 = CollectiveTransition(hc, :S, 1, 2)
        @variables θ::Real ω::Real t::Real r::Real ϕ::Real g::Real Δ::Real
        hi = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        ai = Destroy(hi, :a, 1)
        σi(k, l, m) = IndexedOperator(Transition(hi, :σ, k, l, 2), m)
        i = Index(hi, :i, 3, NLevelSpace(:atom, 2)); j = Index(hi, :j, 3, NLevelSpace(:atom, 2))
        gi = IndexedVariable(:g, i); Jij = DoubleIndexedVariable(:J, i, j)
        Hi = -Δ * ai' * ai + Σ(gi * (ai' * σi(1, 2, i) + ai * σi(2, 1, i)), i)
        Hne = Σ(σi(2, 1, i) * σi(1, 2, j), i, [j])
        avgHi = average(Hi)
        U = Rotation(a, θ)
        Hin = substitute(Σ(g * (ai' * σi(1, 2, i) + ai * σi(2, 1, i)), i), Dict(g => 1.0))
        sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3, 4])
        ψsites = tensor(fockstate(FockBasis(3), 0), nlevelstate(NLevelBasis(2), 1), nlevelstate(NLevelBasis(2), 1), nlevelstate(NLevelBasis(2), 1))
        b = FockBasis(5); ψ = coherentstate(b, 0.3)
        ψq = QuantumToolbox.coherent(6, 0.3); ρq = QuantumToolbox.ket2dm(ψq)
        Hq = to_numeric(ad * a + a, 5; op_type = identity)

        for (name, thunk) in [
                # Construction and accessors
                ("Index", () -> Index(hi, :k, 3, NLevelSpace(:atom, 2))),
                ("IndexedVariable", () -> IndexedVariable(:g, i)),
                ("DoubleIndexedVariable", () -> DoubleIndexedVariable(:J, i, j)),
                ("Σ with non-equal indices", () -> Σ(σi(2, 1, i) * σi(1, 2, j), i, [j])),
                ("∑", () -> ∑(σi(2, 1, i), i)),
                ("Σ Jij", () -> Σ(Jij * σi(2, 1, i) * σi(1, 2, j), i, j)),
                ("FockSpace ⊗ NLevelSpace", () -> FockSpace(:a) ⊗ NLevelSpace(:b, 2)),
                ("UnitaryTransform", () -> UnitaryTransform(ad * a, θ)),
                ("generators", () -> generators(U)),
                ("get_prefactor", () -> get_prefactor(2g * a)),
                ("get_sum_body", () -> get_sum_body(avgHi)),
                ("get_sum_indices", () -> get_sum_indices(average(Σ(σi(2, 1, i), i)))),
                ("get_sum_non_equal", () -> get_sum_non_equal(average(Hne))),
                ("indexed_sum", () -> indexed_sum(get_sum_body(average(Σ(σi(2, 1, i), i))), [i])),
                ("constraint_pairs", () -> constraint_pairs(Hne)),
                ("operator_name / operator_index", () -> (operator_name(a), operator_index(σi(1, 2, i)))),
                ("index accessors", () -> (index_name(i), index_range(i), index_slot(i), index_sym(i), has_index(i))),
                ("order keys", () -> (order_key(a), qadd_order_key(ad * a + a), term_order_key(first(keys((ad * a).arguments))))),
                ("transition_superscript", () -> transition_superscript(true)),
                ("to_num", () -> to_num(to_cnum(2g))),
                ("rename / set_acts_on", () -> (rename(a, :b), set_acts_on(a, 2))),
                ("qadjoint / qconj / dagger", () -> (qadjoint(ad * a + a), qconj(a * ad), dagger(a * ad))),
                ("acts_on(product)", () -> acts_on(ajc' * σjc)),
                ("fundamental_operators / unique_up_to_adjoint", () -> unique_up_to_adjoint(fundamental_operators(hn))),
                # Algebra across operator families
                ("anticommutator(px, py)", () -> anticommutator(px, py)),
                ("px * py * pz", () -> px * py * pz),
                ("sx * sz * sx", () -> sx * sz * sx),
                ("x * p * x", () -> xx * pp * xx),
                ("normal_order(px * sx * px)", () -> normal_order(px * sx * px)),
                ("normal_to_symmetric(x * p * x)", () -> normal_to_symmetric(xx * pp * xx)),
                ("collective S12 * S12'", () -> commutator(S12, S12')),
                ("expand_completeness(σ11 * a)", () -> expand_completeness(σ11 * a)),
                ("exact and symbolic coefficients", () -> (g + Δ)^2 * a * ad + sqrt(2) * g * a + (1 // 2) * ad),
                ("exponential_form / phase_terms", () -> phase_terms(exponential_form(cos(ω * t)))),
                ("expim(ω*t) * expim(θ)", () -> expim(ω * t) * expim(θ)),
                # Indexed and averaged expressions
                ("simplify(non-equal sum)", () -> simplify(Hne)),
                ("commutator(indexed H, a)", () -> commutator(Hi, ai)),
                ("average(non-equal sum)", () -> average(Hne)),
                ("undo_average(indexed average)", () -> undo_average(avgHi)),
                ("substitute(indexed average)", () -> substitute(avgHi, Dict(Δ => 1.0))),
                ("inner_adjoint(compound average)", () -> inner_adjoint(g * average(a) + im * average(ad * a) * average(a))),
                ("inner_adjoint(time-dependent average)", () -> inner_adjoint(make_time_dependent(average(a), t))),
                ("inner_adjoint(indexed average)", () -> inner_adjoint(avgHi)),
                ("get_variables / get_operators / get_indices", () -> (get_variables(Hi), get_operators(Hi), get_indices(Hi))),
                # Printing
                ("repr(indexed H)", () -> repr(Hi)),
                ("repr(indexed average)", () -> repr(avgHi)),
                ("latexify(complex coefficients)", () -> latexify((1 // 2 + 2im) * ad * a - a + (g + im * Δ) * ad)),
                ("latexify(expim(ω*t) * a)", () -> latexify(expim(ω * t) * a)),
                # Numeric backends
                ("numeric_basis / numeric_backend", () -> (numeric_basis(ψ), numeric_backend(ψ), numeric_basis(ψq))),
                ("to_numeric(pauli product)", () -> to_numeric(px * py, SpinBasis(1 // 2))),
                ("to_numeric(spin product)", () -> to_numeric(sx * sz, SpinBasis(1))),
                ("to_numeric(phase-space product)", () -> to_numeric(xx * pp, FockBasis(5))),
                ("to_numeric(a' * a, b; op_type = dense)", () -> to_numeric(ad * a, b; op_type = dense)),
                ("to_numeric(composite; op_type = identity)", () -> to_numeric(ajc' * σjc, FockBasis(3) ⊗ NLevelBasis(2); op_type = identity)),
                ("expect(a' * a, ψ)", () -> expect(ad * a, ψ)),
                ("numeric_average(indexed H, product state)", () -> numeric_average(Hin, ψsites)),
                ("QuantumToolbox to_numeric(nlevel)", () -> to_numeric(σ12 + σ12', 3)),
                ("QuantumToolbox to_numeric(indexed H, dims, sites)", () -> to_numeric(Hin, (4, 2, 2, 2), sites)),
                ("QuantumToolbox op_type = concretize", () -> to_numeric(ad * a, 5; op_type = SO.concretize)),
                ("QuantumToolbox op_type = to_dense", () -> to_numeric(ad * a, 5; op_type = QuantumToolbox.to_dense)),
                ("QuantumToolbox expect(ket)", () -> expect(ad * a, ψq)),
                ("QuantumToolbox expect(VecSum, density)", () -> QuantumToolbox.expect(Hq, ρq)),
                ("QuantumToolbox numeric_average(density)", () -> numeric_average(ad * a, ρq)),
                ("QuantumToolbox concretize(VecSum)", () -> SO.concretize(Hq.data)),
            ]
            rep = JET.@report_opt target_modules = targets thunk()
            @testset "$name" begin
                test_no_reports(rep)
            end
        end
    end
end
