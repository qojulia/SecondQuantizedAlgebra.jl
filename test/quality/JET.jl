using SecondQuantizedAlgebra
using Test
using JET
using QuantumOpticsBase: FockBasis, ManyBodyBasis, NLevelBasis, SpinBasis, LazyKet, basisstate, bosonstates, coherentstate, dense, destroy, fockstate, nlevelstate, tensor
using LinearAlgebra: det, dot, norm
using Symbolics: Equation, Num
using Latexify: latexify
using QuantumToolbox: QuantumToolbox
import SciMLOperators as SO
using SecondQuantizedAlgebra: Coeff, coefficient_terms, SuperOperator, SuperExpr, sprepost, spre, spost, liouvillian, to_numeric_master, constraint_pairs, expect, exponential_form, expim, get_sum_body, get_sum_indices, get_sum_non_equal, indexed_sum, order_key, phase_terms, qadd_order_key, rename, set_acts_on, term_order_key, to_cnum, to_num, transition_superscript, trigonometric_form

const TARGETS = (
    SecondQuantizedAlgebra,
    Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumOpticsBaseExt),
    Base.get_extension(SecondQuantizedAlgebra, :SecondQuantizedAlgebraQuantumToolboxExt),
)

report_text(report) = sprint(show, MIME("text/plain"), report)

# Prints the offending reports on failure: a bare `@test isempty(...)` shows only the count.
function test_no_reports(report)
    reports = map(report_text, JET.get_reports(report))
    isempty(reports) || @info "unexpected JET reports" reports
    @test isempty(reports)
    return
end

# Every entry must be free of errors (`report_call`) and of runtime dispatch (`report_opt`)
# inside the package and both backend extensions.
function test_entries(entries)
    for (name, thunk) in entries
        @testset "$name" begin
            test_no_reports(JET.@report_call target_modules = TARGETS thunk())
            test_no_reports(JET.@report_opt target_modules = TARGETS thunk())
        end
    end
    return
end

# Calls whose runtime dispatch comes from a backend's own type design. Errors must still be
# absent; `report_opt` is marked broken so an upstream fix surfaces as an unexpected pass.
function test_known_dispatch(entries)
    for (name, reason, thunk) in entries
        @testset "$name ($reason)" begin
            test_no_reports(JET.@report_call target_modules = TARGETS thunk())
            @test_broken isempty(JET.get_reports(JET.@report_opt target_modules = TARGETS thunk()))
        end
    end
    return
end

@static if isempty(VERSION.prerelease)
    @testset "JET report_package (correctness)" begin
        test_no_reports(
            JET.report_package(SecondQuantizedAlgebra; target_modules = (SecondQuantizedAlgebra,)),
        )
    end

    @testset "JET report_package ($ext)" for ext in (
            :SecondQuantizedAlgebraQuantumOpticsBaseExt,
            :SecondQuantizedAlgebraQuantumToolboxExt,
        )
        m = Base.get_extension(SecondQuantizedAlgebra, ext)
        @test m isa Module
        test_no_reports(JET.report_package(m; target_modules = (m,)))
    end

    @testset "JET per feature" begin
        # Single-family leaves
        hf = FockSpace(:f); a = Destroy(hf, :a); ad = a'
        hn = NLevelSpace(:atom, 3)
        σ12 = Transition(hn, :σ, 1, 2); σ21 = Transition(hn, :σ, 2, 1); σ11 = Transition(hn, :σ, 1, 1)
        hp = PauliSpace(:s); px = Pauli(hp, :σ, 1); py = Pauli(hp, :σ, 2); pz = Pauli(hp, :σ, 3)
        hs = SpinSpace(:s); sx = Spin(hs, :S, 1); sy = Spin(hs, :S, 2); sz = Spin(hs, :S, 3)
        hph = PhaseSpace(:osc); xx = Position(hph, :x); pp = Momentum(hph, :p)
        hc = CollectiveNLevelSpace(:atom, 3); S12 = CollectiveTransition(hc, :S, 1, 2)
        hl = NLevelSpace(:atom, (:g, :e)); σge = Transition(hl, :σ, :g, :e)

        # Two coupled cavities
        hab = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(hab, :a, 1); b1 = Destroy(hab, :b, 2)

        # Jaynes-Cummings
        hjc = FockSpace(:c) ⊗ NLevelSpace(:atom, 2)
        ajc = Destroy(hjc, :a, 1); σjc = Transition(hjc, :σ, 1, 2, 2)
        Hjc = ajc' * σjc + ajc * σjc'

        @variables g::Real Δ::Real θ::Real ω::Real t::Real r::Real ϕ::Real κ::Real
        @variables α::Number u::Number v::Number η::Real ωd::Real dx::Real dp::Real
        Hp = Δ * ad * a + g * (a + ad)
        H0 = ω * ajc' * ajc

        # Unitary transforms
        U = Rotation(a, θ)
        Ut = Rotation(a, ω * t, t)
        D = Displace(a, α)
        bogo_matrix = [u v; conj(v) conj(u)]
        B = Bogoliubov(a, bogo_matrix)
        fock_reference = ω * ad * a - im * κ * cos(ωd * t) * (a - ad)
        quadrature_reference = (ω / 2) * (xx^2 + pp^2) + η * cos(ωd * t) * xx
        c = to_cnum(g)
        C = [to_cnum(g) to_cnum(1); to_cnum(Δ) to_cnum(2)]

        # Tavis-Cummings with indexed atoms
        hi = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        ai = Destroy(hi, :a, 1)
        σi(k, l, m) = IndexedOperator(Transition(hi, :σ, k, l, 2), m)
        i = Index(hi, :i, 3, NLevelSpace(:atom, 2)); j = Index(hi, :j, 3, NLevelSpace(:atom, 2))
        gi = IndexedVariable(:g, i); Jij = DoubleIndexedVariable(:J, i, j)
        Hi = -Δ * ai' * ai + Σ(g * (ai' * σi(1, 2, i) + ai * σi(2, 1, i)), i)
        Hgi = -Δ * ai' * ai + Σ(gi * (ai' * σi(1, 2, i) + ai * σi(2, 1, i)), i)
        Hne = Σ(σi(2, 1, i) * σi(1, 2, j), i, [j])
        avgHi = average(Hgi)
        Hin = substitute(Hi, Dict(g => 1.0, Δ => 0.5))
        sites = Dict{Int, Vector{Int}}(1 => [1], 2 => [2, 3, 4])

        # Transmon with a Josephson cosine
        htr = PhaseSpace(:transmon) ⊗ FockSpace(:resonator)
        φ = Position(htr, :φ, 1); n = Momentum(htr, :n, 1); ares = Destroy(htr, :a, 2)
        hb = FockSpace(:transmon_boson) ⊗ FockSpace(:resonator)
        bt = Destroy(hb, :b, 1); ct = Destroy(hb, :c, 2)
        @variables E_C::Real E_J::Real φ_zpf::Real n_zpf::Real
        Htr = 4 * E_C * n^2 - E_J * cos(φ) + ω * ares' * ares + g * n * (ares' + ares)
        ladder = Dict(φ => φ_zpf * (bt + bt'), n => im * n_zpf * (bt' - bt))
        Hb = substitute(Htr, ladder)
        F = cos(φ_zpf * (bt + bt')); G = sin(φ_zpf * (bt + bt')); P = expim(θ * bt' * bt)
        hm = FockSpace(:mode); im_ = Index(hm, :i, 4, hm); km = Index(hm, :k, 4, hm)
        am = IndexedOperator(Destroy(hm, :a), im_)
        Fi = cos(am + am')
        Ub = Rotation(bt, θ)
        Ubt = Rotation(bt, ω * t, t)

        # QuantumOptics states and bases
        b = FockBasis(7); bn = NLevelBasis(3); bs = SpinBasis(1 // 2)
        bjc = FockBasis(3) ⊗ NLevelBasis(2)
        bi = FockBasis(3) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2)
        ψ = basisstate(b, 1)
        ψjc = basisstate(bjc, 1)
        ψl = LazyKet(bjc, (fockstate(FockBasis(3), 1), nlevelstate(NLevelBasis(2), 1)))
        ψsites = tensor(
            fockstate(FockBasis(3), 0), nlevelstate(NLevelBasis(2), 1),
            nlevelstate(NLevelBasis(2), 1), nlevelstate(NLevelBasis(2), 1),
        )
        custom = 2.0 * destroy(b)
        bmb = ManyBodyBasis(NLevelBasis(3), bosonstates(NLevelBasis(3), 2))

        # QuantumToolbox states and solver operators
        ψq = QuantumToolbox.fock(8, 1)
        ψspin = QuantumToolbox.basis(2, 0)
        ψqjc = QuantumToolbox.tensor(QuantumToolbox.fock(4, 0), QuantumToolbox.basis(2, 1))
        ρq = QuantumToolbox.ket2dm(QuantumToolbox.coherent(6, 0.3))
        Hq = to_numeric(ad * a + a, 5; op_type = identity)
        Lq = to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => t -> 1 + im * t)).data
        Lp = to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => (p, t) -> p[1] * cos(t))).data
        Ls = to_numeric(ad * a + a, 8; op_type = identity).data
        c1 = to_numeric(sqrt(κ) * a, 8; time_parameter = Dict(κ => (p, t) -> p[1]))
        Ll = QuantumToolbox.liouvillian(
            to_numeric(α * a + conj(α) * ad, 8; time_parameter = Dict(α => t -> cos(t))), [c1],
        ).data
        vq = rand(ComplexF64, 8); wq = similar(vq)
        vl = rand(ComplexF64, 64); wl = similar(vl)

        @testset "Spaces, operators and indices" begin
            test_entries(
                [
                    ("FockSpace ⊗ NLevelSpace", () -> FockSpace(:a) ⊗ NLevelSpace(:b, 2)),
                    ("find_operators(hjc, 2)", () -> find_operators(hjc, 2)),
                    ("fundamental_operators / unique_up_to_adjoint", () -> unique_up_to_adjoint(fundamental_operators(hn))),
                    ("Index", () -> Index(hi, :k, 3, NLevelSpace(:atom, 2))),
                    ("IndexedVariable", () -> IndexedVariable(:g, i)),
                    ("DoubleIndexedVariable", () -> DoubleIndexedVariable(:J, i, j)),
                    ("operator_name / operator_index", () -> (operator_name(a), operator_index(σi(1, 2, i)))),
                    ("index accessors", () -> (index_name(i), index_range(i), index_slot(i), index_sym(i), has_index(i))),
                    ("order keys", () -> (order_key(a), qadd_order_key(ad * a + a), term_order_key(first(keys((ad * a).arguments))))),
                    ("transition_superscript", () -> transition_superscript(true)),
                    ("rename / set_acts_on", () -> (rename(a, :b), set_acts_on(a, 2))),
                    ("acts_on(leaf)", () -> (acts_on(a), acts_on(σ12))),
                    ("acts_on(product)", () -> acts_on(ajc' * σjc)),
                    ("@qnumbers", () -> @qnumbers q::Destroy(hf)),
                    ("Create", () -> (Create(hf, :a), Create(hjc, :a, 1))),
                    ("named levels", () -> (NLevelSpace(:atom, (:g, :e)), Transition(hl, :σ, :g, :e))),
                    ("explicit ground state", () -> NLevelSpace(:atom, 3, 2)),
                    ("collective named levels", () -> CollectiveTransition(CollectiveNLevelSpace(:atom, (:g, :e)), :S, :g, :e)),
                    ("operator predicates", () -> (is_destroy(a), is_create(ad), is_transition(σ12), is_collective_transition(S12), is_pauli(px), is_spin(sx), is_position(xx), is_momentum(pp))),
                    ("optype", () -> optype(a)),
                    ("concrete index", () -> i(2)),
                    ("rename(index)", () -> rename(i, :k)),
                    ("qadjoint / qconj / dagger", () -> (qadjoint(ad * a + a), qconj(a * ad), dagger(a * ad))),
                ],
            )
        end

        @testset "Operator algebra" begin
            test_entries(
                [
                    ("a' * a", () -> ad * a),
                    ("σ12 * σ21", () -> σ12 * σ21),
                    ("px * py", () -> px * py),
                    ("sx * sy", () -> sx * sy),
                    ("x * p", () -> xx * pp),
                    ("px * py * pz", () -> px * py * pz),
                    ("sx * sz * sx", () -> sx * sz * sx),
                    ("x * p * x", () -> xx * pp * xx),
                    ("a * a' * a", () -> a * ad * a),
                    ("a^3", () -> a^3),
                    ("a'^2 * a", () -> ad^2 * a),
                    ("(a + a') * (a - a')", () -> (a + ad) * (a - ad)),
                    ("a' * σ (Fock ⊗ NLevel)", () -> ajc' * σjc),
                    ("Jaynes-Cummings coupling", () -> ajc' * σjc + ajc * σjc'),
                    ("collective S12 * S12'", () -> S12 * S12'),
                    ("commutator(a, a')", () -> commutator(a, ad)),
                    ("commutator(σ12, σ21)", () -> commutator(σ12, σ21)),
                    ("commutator(px, py)", () -> commutator(px, py)),
                    ("commutator(sx, sy)", () -> commutator(sx, sy)),
                    ("commutator(x, p)", () -> commutator(xx, pp)),
                    ("commutator(S12, S12')", () -> commutator(S12, S12')),
                    ("commutator(a * a', a' * a)", () -> commutator(a * ad, ad * a)),
                    ("anticommutator(a, a')", () -> anticommutator(a, ad)),
                    ("anticommutator(px, py)", () -> anticommutator(px, py)),
                    ("sum([a'a, a a', a'a])", () -> sum([ad * a, a * ad, ad * a])),
                    ("reduce(+, [a'a, a a', a'a])", () -> reduce(+, [ad * a, a * ad, ad * a])),
                    ("sum(f, [a + a', a - a'])", () -> sum(x -> x^2, [a + ad, a - ad])),
                    ("scalar shifts", () -> (Hp + 1, 1 + Hp, Hp - g, g - Hp, a + g)),
                    ("scalar division", () -> (Hp / 2, Hp / g, Hp // 2, Hp // g, a / 2, a // 2)),
                    ("(a + a')^3", () -> (a + ad)^3),
                    ("zero / one / isone", () -> (zero(Hp), one(Hp), isone(Hp))),
                    ("term lookup", () -> (length(Hp), Hp[[ad, a]], haskey(Hp, [ad, a]))),
                    ("commutator(a, H)", () -> commutator(a, Hp)),
                ],
            )
        end

        @testset "Canonical forms" begin
            test_entries(
                [
                    ("normal_order(a)", () -> normal_order(a)),
                    ("normal_order(a' a a')", () -> normal_order(ad * a * ad)),
                    ("normal_order(σ21 σ12 σ21)", () -> normal_order(σ21 * σ12 * σ21)),
                    ("normal_order(px * py)", () -> normal_order(px * py)),
                    ("normal_order(px * sx * px)", () -> normal_order(px * sx * px)),
                    ("simplify(a a + a' a')", () -> simplify(a * a + ad * ad)),
                    ("expand_completeness(σ12 σ21)", () -> expand_completeness(σ12 * σ21)),
                    ("expand_completeness(σ11 a)", () -> expand_completeness(σ11 * a)),
                    ("normal_to_symmetric(a' a)", () -> normal_to_symmetric(ad * a)),
                    ("symmetric_to_normal(a' a)", () -> symmetric_to_normal(ad * a)),
                    ("normal_to_symmetric(x p x)", () -> normal_to_symmetric(xx * pp * xx)),
                    ("leaf canonical forms", () -> (simplify(a), expand(a), expand_completeness(σ12), normal_to_symmetric(a), symmetric_to_normal(a))),
                ],
            )
        end

        @testset "Coefficients and phases" begin
            test_entries(
                [
                    ("exact and symbolic coefficients", () -> (g + Δ)^2 * a * ad + sqrt(2) * g * a + (1 // 2) * ad),
                    ("to_num", () -> to_num(to_cnum(2g))),
                    ("get_prefactor", () -> get_prefactor(2g * a)),
                    ("exponential_form(cos(ωt) a)", () -> exponential_form(cos(ω * t) * a)),
                    ("trigonometric_form(expim(ωt) a)", () -> trigonometric_form(expim(ω * t) * a)),
                    ("phase_terms(exponential_form(cos(ωt)))", () -> phase_terms(exponential_form(cos(ω * t)))),
                    ("expim(ωt) * expim(θ)", () -> expim(ω * t) * expim(θ)),
                    ("substitute(expim(ωt), ω => θ)", () -> substitute(expim(ω * t), Dict(ω => θ))),
                    ("real / abs2 of a phase", () -> (real(expim(θ)), abs2(expim(θ)))),
                    ("Coeff arithmetic", () -> (c + 1, c - c, 2c, c / 2, c / c, c^2, -c, conj(c), isone(c), zero(c))),
                    ("Coeff linear algebra", () -> (det(C), norm(Coeff[to_cnum(3), to_cnum(4)]), norm(expim(θ)), dot(c, c))),
                    ("simplify(Coeff)", () -> simplify(c + c)),
                    ("phase forms of operators", () -> (exponential_form(a), trigonometric_form(a))),
                    ("recognized sum and product", () -> (to_cnum(g * Δ + Δ^2 + 2), to_cnum(3 * g * Δ * sqrt(2) * (g + Δ)))),
                    ("to_num with merged radicals", () -> to_num(to_cnum(sqrt(Num(3)) * sqrt(Num(5)) * g))),
                    ("coefficient_terms", () -> coefficient_terms(g^2 + 2g + 3)),
                    ("simplify trigonometric identity", () -> simplify((cos(θ)^2 + sin(θ)^2) * a)),
                ],
            )
        end

        @testset "Substitution" begin
            test_entries(
                [
                    ("substitute(a, a => a')", () -> substitute(a, Dict(a => ad))),
                    ("substitute(Hp, numeric parameters)", () -> substitute(Hp, Dict(g => 1.0, Δ => 2.0))),
                    ("substitute(Hp, a => a + 1)", () -> substitute(Hp, Dict(a => a + 1))),
                ],
            )
        end

        @testset "Unitary transforms and frames" begin
            test_entries(
                [
                    ("Rotation(a, θ)", () -> Rotation(a, θ)),
                    ("UnitaryTransform(a' a, θ)", () -> UnitaryTransform(ad * a, θ)),
                    ("generators(U)", () -> generators(U)),
                    ("inv(U)", () -> inv(U)),
                    ("U * Ut", () -> U * Ut),
                    ("substitute(U, θ => 0.5)", () -> substitute(U, Dict(θ => 0.5))),
                    ("conjugate(a, Rotation)", () -> conjugate(a, U)),
                    ("transform(a' a, timed Rotation)", () -> transform(ad * a, Ut)),
                    ("transform(a' a, Displace)", () -> transform(ad * a, D)),
                    ("transform(a' a, Squeeze)", () -> transform(ad * a, Squeeze(a, r, ϕ))),
                    ("Bogoliubov(a, integer matrix)", () -> Bogoliubov(a, [1 0; 0 1])),
                    ("Bogoliubov(a, complex matrix)", () -> Bogoliubov(a, ComplexF64[im 0; 0 -im])),
                    ("Bogoliubov(a, symbolic matrix)", () -> Bogoliubov(a, bogo_matrix)),
                    ("substitute(Bogoliubov, u, v => exact)", () -> substitute(B, Dict(u => 5 // 3, v => 4 // 3))),
                    ("transform(a' a, Bogoliubov)", () -> transform(ad * a, B)),
                    ("RotatingFrame(H0, t)", () -> RotatingFrame(H0, t)),
                    ("gauge_term(RotatingFrame)", () -> gauge_term(RotatingFrame(H0, t))),
                    ("DisplacementFrame(a, ω a' a)", () -> DisplacementFrame(a, ω * ad * a)),
                    ("DisplacementFrame(a, driven H, t)", () -> DisplacementFrame(a, fock_reference, t)),
                    ("DisplacementFrame(x, p, H)", () -> DisplacementFrame(xx, pp, (ω / 2) * (xx^2 + pp^2) + η * xx)),
                    ("DisplacementFrame(x, p, driven H, t)", () -> DisplacementFrame(xx, pp, quadrature_reference, t)),
                    ("timed Displace / Squeeze", () -> (Displace(a, η * t, t), Squeeze(a, r * t, 0, t))),
                    ("two-mode Rotation / Squeeze", () -> (Rotation(a1, b1, θ), Squeeze(a1, b1, r))),
                    ("timed two-mode Rotation / Squeeze", () -> (Rotation(a1, b1, ω * t, t), Squeeze(a1, b1, r * t, t))),
                    ("quadrature Rotation / Squeeze / Displace", () -> (Rotation(xx, pp, θ), Squeeze(xx, pp, r), Displace(xx, pp, dx, dp))),
                    ("timed quadrature Displace", () -> Displace(xx, pp, η * t, ω * t^2, t)),
                    ("spin Rotation", () -> (Rotation(sx, 3, θ), Rotation(sx, 3, t, t))),
                    ("N-level Rotation", () -> (Rotation(σge, [0 1; 1 0]), Rotation(σge, [0 1; 1 0], t))),
                    ("multimode Bogoliubov", () -> (Bogoliubov((a1, b1), [1 0; 0 1], [0 0; 0 0]), Bogoliubov(a, [1;;], [0;;]))),
                    ("leaf generators", () -> (UnitaryTransform(pz, θ), RotatingFrame(pz, t))),
                    ("adjoint(U)", () -> U'),
                    ("transform(a, U)", () -> transform(a, U)),
                ],
            )
        end

        @testset "Indexed sums" begin
            test_entries(
                [
                    ("Σ(σ21_i) * Σ(σ12_j)", () -> Σ(σi(2, 1, i), i) * Σ(σi(1, 2, j), j)),
                    ("Σ with non-equal indices", () -> Σ(σi(2, 1, i) * σi(1, 2, j), i, [j])),
                    ("∑", () -> ∑(σi(2, 1, i), i)),
                    ("Σ J_ij", () -> Σ(Jij * σi(2, 1, i) * σi(1, 2, j), i, j)),
                    ("simplify(indexed H)", () -> simplify(Hi)),
                    ("simplify(non-equal sum)", () -> simplify(Hne)),
                    ("commutator(indexed H, a)", () -> commutator(Hgi, ai)),
                    ("change_index(indexed H, i, j)", () -> change_index(Hi, i, j)),
                    ("assume_distinct_index(σ21_i σ12_j)", () -> assume_distinct_index(σi(2, 1, i) * σi(1, 2, j), [(i, j)])),
                    ("constraint_pairs", () -> constraint_pairs(Hne)),
                    ("get_variables / get_operators / get_indices", () -> (get_variables(Hgi), get_operators(Hgi), get_indices(Hgi))),
                    ("change_index(indexed H, Dict)", () -> change_index(Hi, Dict(i => j))),
                    ("get_indices(leaf)", () -> get_indices(σi(2, 1, i))),
                ],
            )
        end

        @testset "Averages" begin
            test_entries(
                [
                    ("average(leaf) per family", () -> (average(a), average(σ12), average(px), average(sx), average(xx))),
                    ("average(a a')", () -> average(a * ad)),
                    ("undo_average(average(a))", () -> undo_average(average(a))),
                    ("undo_average(average(a a'))", () -> undo_average(average(a * ad))),
                    ("is_average", () -> (is_average(a), is_average(average(a)))),
                    ("has_sum_metadata", () -> SecondQuantizedAlgebra.has_sum_metadata(average(a))),
                    ("average(indexed H)", () -> average(Hi)),
                    ("average(non-equal sum)", () -> average(Hne)),
                    ("undo_average(indexed average)", () -> undo_average(avgHi)),
                    ("substitute(indexed average)", () -> substitute(avgHi, Dict(Δ => 1.0))),
                    ("make_time_dependent(indexed average)", () -> make_time_dependent(average(Hi), t)),
                    ("get_sum_body", () -> get_sum_body(avgHi)),
                    ("get_sum_indices", () -> get_sum_indices(average(Σ(σi(2, 1, i), i)))),
                    ("get_sum_non_equal", () -> get_sum_non_equal(average(Hne))),
                    ("indexed_sum", () -> indexed_sum(get_sum_body(average(Σ(σi(2, 1, i), i))), [i])),
                    ("inner_adjoint(compound average)", () -> inner_adjoint(g * average(a) + im * average(ad * a) * average(a))),
                    ("inner_adjoint(time-dependent average)", () -> inner_adjoint(make_time_dependent(average(a), t))),
                    ("inner_adjoint(indexed average)", () -> inner_adjoint(avgHi)),
                    ("average(scalar)", () -> (average(2.0), average(g))),
                    ("undo_average(Equation)", () -> undo_average(Equation(average(a), average(ad)))),
                    ("make_time_dependent(sum of averages)", () -> make_time_dependent(average(a) + average(ad), t)),
                    ("acts_on / get_indices(average)", () -> (acts_on(average(a)), get_indices(average(σi(2, 1, i))))),
                    ("is_indexed_sum", () -> is_indexed_sum(average(Σ(σi(2, 1, i), i)))),
                ],
            )
        end

        @testset "Formal operator functions" begin
            test_entries(
                [
                    ("cos(φ)", () -> cos(φ)),
                    ("sin(φ_zpf (b + b'))", () -> sin(φ_zpf * (bt + bt'))),
                    ("expim(θ b' b)", () -> expim(θ * bt' * bt)),
                    ("transmon H", () -> 4 * E_C * n^2 - E_J * cos(φ) + ω * ares' * ares + g * n * (ares' + ares)),
                    ("substitute(H, φ, n => ladder)", () -> substitute(Htr, ladder)),
                    ("taylor(Hb, 0:4)", () -> taylor(Hb, 0:4)),
                    ("taylor(expim(θ b' b), 0:3)", () -> taylor(P, 0:3)),
                    ("cos + sin", () -> F + G),
                    ("cos - sin", () -> F - G),
                    ("cos * sin", () -> F * G),
                    ("cos^2", () -> F^2),
                    ("θ * cos", () -> θ * F),
                    ("cos * c (disjoint mode)", () -> F * ct),
                    ("c * cos (disjoint mode)", () -> ct * F),
                    ("(θ b + θ² b') * cos", () -> (θ * bt + θ^2 * bt') * F),
                    ("cos + b' b", () -> F + bt' * bt),
                    ("adjoint(cos * b)", () -> (F * bt)'),
                    ("-cos", () -> -F),
                    ("cos(cos)", () -> cos(F)),
                    ("isequal / hash / iszero", () -> (isequal(F, F), hash(F), iszero(F - F))),
                    ("substitute(θ cos, θ => 1)", () -> substitute(θ * F, Dict(θ => 1))),
                    ("conjugate(cos, Rotation)", () -> conjugate(F, Ub)),
                    ("transform(expim, Rotation)", () -> transform(P, Ub)),
                    ("change_index(cos(a_i + a_i'))", () -> change_index(Fi, im_, km)),
                    ("get_variables(θ cos)", () -> get_variables(θ * F)),
                    ("normal_order / simplify / expand", () -> (normal_order(F), simplify(F), expand(F))),
                    ("commutator(cos, b)", () -> commutator(F, bt)),
                    ("commutator(b, cos) / (cos, sin) / (cos, 2)", () -> (commutator(bt, F), commutator(F, G), commutator(F, 2))),
                    ("cos ± operator / scalar", () -> (F + bt, bt + F, F + 1, 1 + F, F - g)),
                    ("cos * polynomial / scalar", () -> (F * (bt + bt'), F * g)),
                    ("cos / scalar", () -> (F / 2, F / g, F // 2, F // g)),
                    ("zero / one / isone(cos)", () -> (zero(F), one(F), isone(F))),
                    ("expand_completeness / assume_distinct_index", () -> (expand_completeness(F), assume_distinct_index(Fi, [(im_, km)]))),
                    ("get_indices / get_operators(cos)", () -> (get_indices(Fi), get_operators(F))),
                    ("change_index(cos, Dict)", () -> change_index(Fi, Dict(im_ => km))),
                    ("transform(expim, timed Rotation)", () -> transform(P, Ubt)),
                    ("scalar * (cos * polynomial)", () -> 2 * (F * bt)),
                ],
            )
        end

        @testset "Printing" begin
            test_entries(
                [
                    ("repr(Hp)", () -> repr(Hp)),
                    ("repr(indexed H)", () -> repr(Hgi)),
                    ("repr(indexed average)", () -> repr(avgHi)),
                    ("repr(θ cos * c)", () -> repr(θ * F * ct)),
                    ("latexify(Hp)", () -> latexify(Hp)),
                    ("latexify(complex coefficients)", () -> latexify((1 // 2 + 2im) * ad * a - a + (g + im * Δ) * ad)),
                    ("latexify(expim(ωt) a)", () -> latexify(expim(ω * t) * a)),
                    ("latexify(average(indexed H))", () -> latexify(average(Hi))),
                    ("latexify(2 cos)", () -> latexify(2 * F)),
                    ("latexify(θ cos * c)", () -> latexify(θ * F * ct)),
                    ("repr(Hilbert spaces)", () -> (repr(hjc), repr(hn), repr(hc), repr(hp), repr(hs), repr(hph))),
                    ("repr(Index)", () -> repr(i)),
                    ("repr(UnitaryTransform)", () -> (repr(U), repr(MIME("text/plain"), Ut))),
                    ("repr(cos + sin)", () -> repr(F + G)),
                    ("latexify(cos + sin)", () -> latexify(F + G)),
                    ("text/latex MIME", () -> (repr(MIME("text/latex"), Hp), repr(MIME("text/latex"), [a, ad]))),
                ],
            )
        end

        @testset "Numeric: QuantumOptics" begin
            test_entries(
                [
                    ("to_numeric(leaf) per family", () -> (to_numeric(a, b), to_numeric(σ12, bn), to_numeric(px, bs), to_numeric(sx, bs), to_numeric(xx, b))),
                    ("to_numeric(a' a, b)", () -> to_numeric(ad * a, b)),
                    ("to_numeric(pauli product)", () -> to_numeric(px * py, SpinBasis(1 // 2))),
                    ("to_numeric(spin product)", () -> to_numeric(sx * sz, SpinBasis(1))),
                    ("to_numeric(phase-space product)", () -> to_numeric(xx * pp, FockBasis(5))),
                    ("to_numeric(2.0, b)", () -> to_numeric(2.0, b)),
                    ("to_numeric([a, a'], b)", () -> to_numeric([a, ad], b)),
                    ("op_type = identity", () -> to_numeric(ad * a, b; op_type = identity)),
                    ("op_type = dense", () -> to_numeric(ad * a, b; op_type = dense)),
                    ("parameter", () -> to_numeric(Hp, b; parameter = Dict(g => 1.0, Δ => 0.5))),
                    ("operators", () -> to_numeric(ad * a, b; operators = Dict(a => custom))),
                    ("time_parameter", () -> to_numeric(α * a + conj(α) * ad, b; time_parameter = Dict(α => t -> 1 + im * t))),
                    ("space + cutoff with QuantumOpticsBackend", () -> to_numeric(ad * a, hf, 7; backend = QuantumOpticsBackend())),
                    ("composite basis", () -> to_numeric(Hjc, bjc)),
                    ("composite basis, op_type = identity", () -> to_numeric(ajc' * σjc, bjc; op_type = identity)),
                    ("indexed H with sites", () -> to_numeric(Hin, bi, sites)),
                    ("numeric_average(average(a), ψ)", () -> numeric_average(average(a), ψ)),
                    ("numeric_average(a' a, [ψ, ψ])", () -> numeric_average(ad * a, [ψ, ψ])),
                    ("numeric_average(Hjc, ψjc)", () -> numeric_average(Hjc, ψjc)),
                    ("numeric_average(Hjc, LazyKet)", () -> numeric_average(Hjc, ψl)),
                    ("numeric_average(indexed H, product state)", () -> numeric_average(Hin, ψsites)),
                    ("expect(a' a, ψ)", () -> expect(ad * a, ψ)),
                    ("numeric_basis / numeric_backend", () -> (numeric_basis(ψ), numeric_backend(ψ))),
                    ("to_numeric(op, state)", () -> to_numeric(ad * a, ψ)),
                    ("to_numeric(scalar, state)", () -> to_numeric(2.0, ψ)),
                    ("to_numeric(symbolic scalar)", () -> to_numeric(g, b; parameter = Dict(g => 1.0))),
                    ("leaf per space with QuantumOpticsBackend", () -> (to_numeric(σ12, hn, 0; backend = QuantumOpticsBackend()), to_numeric(px, hp, 0; backend = QuantumOpticsBackend()), to_numeric(xx, hph, 6; backend = QuantumOpticsBackend()))),
                    ("collective transition on ManyBodyBasis", () -> to_numeric(S12, bmb)),
                    ("numeric_average(sum and product of averages)", () -> (numeric_average(average(a) + average(ad * a), ψ), numeric_average(average(a) * average(ad), ψ))),
                    ("backend hooks", () -> (numeric_operator(QuantumOpticsBackend(), a, b), numeric_subbasis(QuantumOpticsBackend(), bjc, 1), numeric_embed(QuantumOpticsBackend(), bjc, 1, destroy(FockBasis(3))), numeric_identity(QuantumOpticsBackend(), bjc), numeric_num_subsystems(QuantumOpticsBackend(), bjc), numeric_expect(QuantumOpticsBackend(), to_numeric(ad * a, b), ψ))),
                ],
            )
        end

        @testset "Numeric: QuantumToolbox" begin
            test_entries(
                [
                    ("to_numeric(a' a, 8)", () -> to_numeric(ad * a, 8)),
                    ("to_numeric(nlevel)", () -> to_numeric(σ12 + σ12', 3)),
                    ("to_numeric(Hjc, (4, 2))", () -> to_numeric(Hjc, (4, 2))),
                    ("space + dims with QuantumToolboxBackend", () -> to_numeric(Hjc, hjc, (3, 2); backend = QuantumToolboxBackend())),
                    ("indexed H with dims and sites", () -> to_numeric(Hin, (4, 2, 2, 2), sites)),
                    ("op_type = concretize", () -> to_numeric(ad * a, 5; op_type = SO.concretize)),
                    ("op_type = to_dense", () -> to_numeric(ad * a, 5; op_type = QuantumToolbox.to_dense)),
                    ("solver apply: time-dependent QobjEvo", () -> Lq(wq, vq, nothing, nothing, 0.3)),
                    ("solver apply: p-aware QobjEvo", () -> Lp(wq, vq, nothing, [0.5], 0.3)),
                    ("solver apply: static lazy QobjEvo", () -> Ls(wq, vq, nothing, nothing, 0.3)),
                    ("solver apply: Liouvillian with time-dependent collapse", () -> Ll(wl, vl, nothing, [0.5], 0.3)),
                    ("numeric_average(a' a, ψq)", () -> numeric_average(ad * a, ψq)),
                    ("numeric_average(average(a), ψq)", () -> numeric_average(average(a), ψq)),
                    ("numeric_average(Hjc, ψqjc)", () -> numeric_average(Hjc, ψqjc)),
                    ("numeric_average(a' a, [ψq, ψq])", () -> numeric_average(ad * a, [ψq, ψq])),
                    ("numeric_average(density)", () -> numeric_average(ad * a, ρq)),
                    ("expect(a' a, ket)", () -> expect(ad * a, ψq)),
                    ("expect(VecSum, density)", () -> QuantumToolbox.expect(Hq, ρq)),
                    ("concretize(VecSum)", () -> SO.concretize(Hq.data)),
                    ("numeric_basis", () -> numeric_basis(ψq)),
                    ("Pauli / Spin / phase-space averages", () -> (numeric_average(px, ψspin), numeric_average(sx, ψspin), numeric_average(xx, ψq))),
                    ("backend hooks", () -> (numeric_subbasis(QuantumToolboxBackend(), (4, 2), 1), numeric_embed(QuantumToolboxBackend(), (4, 2), 1, QuantumToolbox.destroy(4)), numeric_identity(QuantumToolboxBackend(), (4, 2)), numeric_num_subsystems(QuantumToolboxBackend(), [4, 2]), numeric_identity(QuantumToolboxBackend(), [4, 2]))),
                ],
            )
        end

        @testset "Superoperators and Lindbladians" begin
            Asup = sprepost(a, ad)
            Bsup = spre(ad * a) + spost(a)
            Lsup = hamiltonian_action(Δ * ad * a) + κ * dissipator(a)
            Lnum = hamiltonian_action(2 * ad * a) + 0.5 * dissipator(a)
            ha2 = NLevelSpace(:atom, 2)
            is_ = Index(ha2, :i, 3, ha2); js_ = Index(ha2, :j, 3, ha2)
            σs(x, y, q) = IndexedOperator(Transition(ha2, :σ, x, y), q)
            Lind = Σ(dissipator(σs(1, 2, is_)), is_)
            Lfor = hamiltonian_action(cos(xx))
            Xsup = SuperExpr(Lfor)
            Dsup = Displace(a, α)
            Lmix = κ * dissipator(cos(xx)) + dissipator(xx)
            test_entries(
                [
                    ("sprepost / spre / spost", () -> (sprepost(a, ad), spre(a), spost(a))),
                    ("SuperOperator(pairs, indices)", () -> SuperOperator(collect(Asup), Index[])),
                    ("hamiltonian_action / dissipator", () -> (hamiltonian_action(Hp), dissipator(a), dissipator(a, ad))),
                    ("collective dissipator", () -> dissipator(σs(1, 2, is_), DoubleIndexedVariable(:Γ, is_, js_), is_ => js_)),
                    ("liouvillian", () -> liouvillian(Hp, [a, ad]; rates = [κ, 1])),
                    ("map + / scalar *", () -> (Asup + Bsup, g * Asup, Lsup + Lsup, 2 * Lsup)),
                    ("composition", () -> (Asup ∘ Bsup, Lsup ∘ Lsup)),
                    ("apply", () -> (Asup(ad * a), Lsup(ad * a))),
                    ("transpose / adjoint", () -> (transpose(Asup), Asup', transpose(Lsup))),
                    ("SuperOperator(L)", () -> SuperOperator(Lsup)),
                    ("Σ / change_index", () -> (Σ(dissipator(σs(1, 2, is_)), is_), change_index(Lind, is_, js_))),
                    ("simplify maps", () -> (simplify(Lsup), simplify(Asup))),
                    ("transform maps", () -> (transform(Lsup, Dsup), transform(Lsup, Ut), transform(Asup, Dsup))),
                    ("SuperExpr", () -> (SuperExpr(Lfor), Xsup + Xsup, taylor(Xsup, 0:2))),
                    ("formal and polynomial jumps", () -> (Lmix, dissipator(cos(xx), xx + 1), SuperExpr(Lmix))),
                    ("hamiltonian_action(cos(x) + 1)", () -> hamiltonian_action(cos(xx) + 1)),
                    ("Lindbladian(SuperOperator)", () -> (Lindbladian(SuperOperator(Lsup)), Lindbladian(SuperOperator(Lnum)))),
                    ("to_numeric(L, basis)", () -> to_numeric(Lnum, FockBasis(4))),
                    ("to_numeric(L, dims)", () -> to_numeric(Lnum, 5)),
                    ("to_numeric_master QuantumOptics", () -> to_numeric_master(Lnum, hf, 4; backend = QuantumOpticsBackend())),
                    ("to_numeric_master QuantumToolbox", () -> to_numeric_master(Lnum, hf, 4; backend = QuantumToolboxBackend())),
                ],
            )
        end

        @testset "Known runtime dispatch" begin
            test_known_dispatch(
                [
                    (
                        "SpinSpace with runtime spin", "SpinBasis{S} carries the spin as a type parameter",
                        () -> to_numeric(sx, hs, 1 // 2; backend = QuantumOpticsBackend()),
                    ),
                    (
                        "QuantumToolbox Vector dims", "QuantumToolbox dimensions carry their length in the type",
                        () -> to_numeric(Hjc, [4, 2]),
                    ),
                ],
            )
        end
    end
end
