"""
    benchmark/quantumalgebra_comparison.jl

Fair head-to-head benchmark of SecondQuantizedAlgebra.jl (SQA) against
QuantumAlgebra.jl (QA), built as a `BenchmarkTools.BenchmarkGroup` suite like
the rest of `benchmark/`.

Fairness contract
-----------------
Both packages are timed *producing the same canonical result from the same
physical input*, written idiomatically in each:

  * SQA canonicalizes **eagerly** — every `*`/`+` already returns the
    normal-ordered, fully-simplified expression. So for SQA we time the
    construction/arithmetic itself.
  * QA canonicalizes **lazily** — products stay symbolic until `normal_form`
    is called. So for QA we time the same arithmetic wrapped in `normal_form`
    (and `normal_form ∘ comm` for commutators), which is the idiomatic way to
    reach the canonical form SQA hands you for free. (QA also offers
    `auto_normal_form(true)` for eager behaviour; see `auto_normal_form_demo`.)

Only operations expressible in *both* packages are compared head-to-head.
QA's two-level support is Pauli / spin-½ only (no general N-level transitions,
no arbitrary spin-S), so those SQA capabilities are shown separately in the
docs page rather than as rigged timing wins. The σⁱ basis also differs
(QA expands `σ⁺σ⁻` into σx/σy/σz), so results are physically equivalent but
not term-by-term identical.

The suite is laid out as `SUITE[group][name]["SQA" | "QA"]`, so SQA and QA are
paired leaves under each benchmark and can be tabulated side by side. The
"Scaling / large systems" group sweeps a size parameter and is plotted.

Usage:
    julia --project=benchmark benchmark/quantumalgebra_comparison.jl
"""

using BenchmarkTools
using SecondQuantizedAlgebra
import QuantumAlgebra as QA
import CairoMakie as Makie
using Symbolics: @variables

# A second bosonic species for the small two-mode benchmarks (QA's `a()` is
# the first); the many-mode sweep instead uses integer-indexed modes a(k).
QA.@boson_ops bb

# Eager product: normal_form after every multiply — the same incremental
# canonicalization SQA does on every `*` (equivalent to `auto_normal_form(true)`,
# but local). This is the FAIR QA workflow for powers/products: comparing SQA's
# eager `*` against QA's lazy build-then-normalize-once would penalise QA for a
# workflow choice rather than an algorithmic difference (a lazy `Hⁿ` must expand
# the whole product before collapsing, which blows up). `qa_lazy_pow` keeps the
# default lazy workflow around for the honest "blowup footgun" callout.
qa_eager_prod(factors) = foldl((acc, h) -> QA.normal_form(acc * h), factors)
qa_lazy_pow(H, n) = QA.normal_form(foldl(*, fill(H, n)))

# ── Registries ───────────────────────────────────────────────────────────────
# `ORDER` records insertion order so the report table is deterministic;
# `BenchmarkGroup` itself is dict-backed and unordered.
const ORDER = Tuple{String, String}[]
# Benchmarks skipped by the time cap: (group, name) => trial seconds.
const CAPPED = Dict{Tuple{String, String}, Float64}()
# Size sweeps for plotting: label => [(x, group, name), …] (insertion order).
const SWEEPS = Tuple{String, Vector{Tuple{Float64, String, String}}}[]

function _sweep_points(label)
    for (l, pts) in SWEEPS
        l == label && return pts
    end
    pts = Tuple{Float64, String, String}[]
    push!(SWEEPS, (label, pts))
    return pts
end

"""
    pair!(SUITE, group, name, sqa_thunk, qa_thunk)

Register the equivalent SQA and QA thunks as paired leaves
`SUITE[group][name]["SQA"]` and `SUITE[group][name]["QA"]`. The thunk value is
interpolated into `@benchmarkable` so neither side pays for a global lookup.
"""
function pair!(SUITE, group, name, sqa_thunk, qa_thunk)
    haskey(SUITE, group) || (SUITE[group] = BenchmarkGroup())
    SUITE[group][name] = BenchmarkGroup()
    SUITE[group][name]["SQA"] = @benchmarkable $sqa_thunk()
    SUITE[group][name]["QA"] = @benchmarkable $qa_thunk()
    push!(ORDER, (group, name))
    return true
end

# One warm-up call (compilation), then one timed call.
_trial(thunk) = (thunk(); @elapsed thunk())

"""
    safe_pair!(SUITE, group, name, sqa_thunk, qa_thunk; cap = 3.0)

Like [`pair!`](@ref) but first times a single evaluation of each side. If
either exceeds `cap` seconds the benchmark is *not* registered (so an eager
term-blowup degrades gracefully instead of stalling the suite); it is recorded
in `CAPPED` and shown as such in the table. Returns `true` if registered.
"""
function safe_pair!(SUITE, group, name, sqa_thunk, qa_thunk; cap = 3.0)
    ts = _trial(sqa_thunk)
    tq = _trial(qa_thunk)
    t = max(ts, tq)
    println(stderr, "  trial: $name — SQA $(round(ts * 1e3; digits = 2)) ms, QA $(round(tq * 1e3; digits = 2)) ms")
    flush(stderr)
    if t > cap
        push!(ORDER, (group, name))
        CAPPED[(group, name)] = t
        return false
    end
    return pair!(SUITE, group, name, sqa_thunk, qa_thunk)
end

function sweep_pair!(SUITE, label, x, group, name, sqa_thunk, qa_thunk; cap = 3.0)
    ok = safe_pair!(SUITE, group, name, sqa_thunk, qa_thunk; cap = cap)
    ok && push!(_sweep_points(label), (Float64(x), group, name))
    return ok
end

# ── Equivalence validation ────────────────────────────────────────────────────
# Representation-robust, operator-vs-operator known-answer identities, checked
# independently in each package. This validates that we are comparing two
# *correct* implementations — the premise of a fair comparison.
function validate_equivalence()
    hf = FockSpace(:f)
    c = Destroy(hf, :c)
    ht = NLevelSpace(:a, 2, 1)
    σ⁻ = Transition(ht, :σ, 1, 2)
    σee = Transition(ht, :σ, 2, 2)
    σgg = Transition(ht, :σ, 1, 1)

    checks = [
        (
            "[a, a†] · canonical (aa† = 1 + a†a)",
            isequal(c * c', 1 + c' * c),
            QA.normal_form(QA.a() * QA.a'()) == QA.normal_form(QA.QuExpr(1) + QA.a'() * QA.a()),
        ),
        (
            "cubic reorder (aaa† = a†aa + 2a)",
            isequal(c * c * c', c' * c * c + 2 * c),
            QA.normal_form(QA.a() * QA.a() * QA.a'()) ==
                QA.normal_form(QA.a'() * QA.a() * QA.a() + 2 * QA.a()),
        ),
        (
            "two-level completeness (σgg + σee = 1)",
            isequal(σgg + σee, σ⁻ * σ⁻' + σ⁻' * σ⁻),
            QA.normal_form(QA.σm() * QA.σp() + QA.σp() * QA.σm()) == QA.normal_form(QA.QuExpr(1)),
        ),
    ]

    println("Equivalence validation (known-answer identities):")
    allok = true
    for (desc, sqa_ok, qa_ok) in checks
        allok &= sqa_ok & qa_ok
        println("  ", (sqa_ok & qa_ok) ? "✓" : "✗", "  SQA:$(sqa_ok) QA:$(qa_ok)  $desc")
    end
    allok || @warn "Some equivalence identities failed — comparison fairness is suspect."
    return allok
end

# ── Suite ──────────────────────────────────────────────────────────────────────

function benchmark_qa_comparison!(SUITE)
    # ── Group 1 — Core algebra (point comparisons) ───────────────────────────
    let
        hf = FockSpace(:c)
        ha = NLevelSpace(:a, 2, 1)
        hjc = hf ⊗ ha
        a = Destroy(hjc, :a, 1)
        σge = Transition(hjc, :σ, 1, 2, 2)   # |g⟩⟨e| = σ⁻
        σee = Transition(hjc, :σ, 2, 2, 2)   # |e⟩⟨e|
        @variables ω_c ω_a g

        Hjc_s() = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')
        H_s = Hjc_s()

        ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
        Hjc_q_lazy = ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
            gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
        Hjc_q() = QA.normal_form(
            ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
                gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
        )

        pair!(SUITE, "Core algebra", "Jaynes–Cummings: build H", Hjc_s, Hjc_q)
        pair!(
            SUITE, "Core algebra", "Jaynes–Cummings: H²",
            () -> H_s * H_s, () -> QA.normal_form(Hjc_q_lazy * Hjc_q_lazy),
        )

        # Schrieffer–Wolff. A single symbolic prefactor χ on both sides isolates
        # operator algebra from scalar-rational handling.
        @variables ε χ
        H0_s = ω_c * a' * a - ε * σee
        V_s = g * (a' * σge' + a * σge)
        S_s = χ * (a' * σge' - a * σge)

        εq = QA.Pr"ε"; χq = QA.Pr"χ"
        H0_q = ωcq * QA.a'() * QA.a() - εq * QA.σp() * QA.σm()
        V_q = gq * (QA.a'() * QA.σp() + QA.a() * QA.σm())
        S_q = χq * (QA.a'() * QA.σp() - QA.a() * QA.σm())

        pair!(
            SUITE, "Core algebra", "Schrieffer–Wolff: [S, V]",
            () -> commutator(S_s, V_s), () -> QA.normal_form(QA.comm(S_q, V_q)),
        )
        pair!(
            SUITE, "Core algebra", "Schrieffer–Wolff: [S, [S, H₀]]",
            () -> commutator(S_s, commutator(S_s, H0_s)),
            () -> QA.normal_form(QA.comm(S_q, QA.normal_form(QA.comm(S_q, H0_q)))),
        )
    end

    let
        h2c = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2c, :a, 1)
        a2 = Destroy(h2c, :b, 2)
        @variables J ω1 ω2

        H2_s() = ω1 * a1' * a1 + ω2 * a2' * a2 + J * (a1' * a2 + a2' * a1)
        H_s = H2_s()

        ω1q = QA.Pr"ω1"; ω2q = QA.Pr"ω2"; Jq = QA.Pr"J"
        H2_q_lazy = ω1q * QA.a'() * QA.a() + ω2q * bb'() * bb() + Jq * (QA.a'() * bb() + bb'() * QA.a())
        H2_q() = QA.normal_form(
            ω1q * QA.a'() * QA.a() + ω2q * bb'() * bb() + Jq * (QA.a'() * bb() + bb'() * QA.a())
        )

        pair!(SUITE, "Core algebra", "Two cavities: build H", H2_s, H2_q)
        pair!(
            SUITE, "Core algebra", "Two cavities: H²",
            () -> H_s * H_s, () -> QA.normal_form(H2_q_lazy * H2_q_lazy),
        )
        pair!(
            SUITE, "Core algebra", "Multi-mode 6-op chain",
            () -> a1 * a2 * a1' * a2' * a1 * a2',
            () -> QA.normal_form(QA.a() * bb() * QA.a'() * bb'() * QA.a() * bb'()),
        )
    end

    # ── Group 2 — Indexed sums ───────────────────────────────────────────────
    @variables N

    let
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        a = Destroy(h, :a, 1)
        i = Index(h, :i, N, NLevelSpace(:atom, 2, 1))
        gi = IndexedVariable(:g, i)
        σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

        tc_s() = Σ(gi * (a' * σ(1, 2) + a * σ(2, 1)), i)
        tc_q() = QA.normal_form(QA.∑(:i, QA.Pr"g_i" * (QA.a'() * QA.σm(:i) + QA.a() * QA.σp(:i))))
        pair!(SUITE, "Indexed sums", "Tavis–Cummings Σ construction", tc_s, tc_q)
    end

    let
        hd = FockSpace(:c) ⊗ PauliSpace(:s)
        b = Destroy(hd, :b, 1)
        @variables ω λ
        id = Index(hd, :i, N, PauliSpace(:s))
        jd = Index(hd, :j, N, PauliSpace(:s))
        Sx(idx) = IndexedOperator(Pauli(hd, :σ, 1), idx)
        Sz(idx) = IndexedOperator(Pauli(hd, :σ, 3), idx)

        Hd_s() = ω * b' * b + Σ(ω * Sz(id) + λ * (b + b') * Sx(id), id)
        H_s = Hd_s()
        ωq = QA.Pr"ω"; λq = QA.Pr"λ"
        Hd_q() = QA.normal_form(
            ωq * QA.a'() * QA.a() + QA.∑(:i, ωq * QA.σz(:i) + λq * (QA.a() + QA.a'()) * QA.σx(:i))
        )
        H_q = Hd_q()

        pair!(SUITE, "Indexed sums", "Dicke Σ construction", Hd_s, Hd_q)
        pair!(
            SUITE, "Indexed sums", "Dicke diagonal-collapse [H, σ_j]",
            () -> commutator(H_s, Sz(jd)), () -> QA.normal_form(QA.comm(H_q, QA.σz(:j))),
        )
    end

    # Large indexed double-sum: Σᵢ Σⱼ Jᵢⱼ σˣᵢ σˣⱼ.
    let
        hds = FockSpace(:c) ⊗ SpinSpace(:s)
        id = Index(hds, :i, N, SpinSpace(:s))
        jd = Index(hds, :j, N, SpinSpace(:s))
        Sx(idx) = IndexedOperator(Spin(hds, :S, 1, 2), idx)
        Jij = DoubleIndexedVariable(:J, id, jd; identical = false)

        ds_s() = Σ(Σ(Jij * Sx(id) * Sx(jd), id), jd)
        ds_q() = QA.normal_form(QA.∑(:i, QA.∑(:j, QA.Pr"J_i,j" * QA.σx(:i) * QA.σx(:j))))
        pair!(SUITE, "Indexed sums", "Double Σ_ij spin–spin", ds_s, ds_q)
    end

    # ── Group 3 — Expectation values ─────────────────────────────────────────
    let
        h = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a = Destroy(h, :a, 1)
        σge = Transition(h, :σ, 1, 2, 2)
        σee = Transition(h, :σ, 2, 2, 2)
        @variables ω_c ω_a g
        H_s = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')

        ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
        H_q = QA.normal_form(
            ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
                gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
        )

        pair!(
            SUITE, "Expectation values", "⟨a† σ⁻⟩",
            () -> average(a' * σge), () -> QA.expval(QA.a'() * QA.σm()),
        )
        pair!(SUITE, "Expectation values", "⟨H_JC⟩", () -> average(H_s), () -> QA.expval(H_q))
    end

    # ── Group 4 — Scaling / large systems (swept + plotted) ──────────────────

    # Fock normal-ordering: (a · a†)ⁿ for growing n.
    let
        hf = FockSpace(:f)
        c = Destroy(hf, :c)
        for n in 2:10
            # QA uses its eager workflow here (fair match to SQA's eager `*`).
            sweep_pair!(
                SUITE, "Fock (a·a†)ⁿ", n, "Scaling / large systems", "Fock (a·a†)ⁿ n=$n",
                () -> (c * c')^n, () -> qa_eager_prod(fill(QA.a() * QA.a'(), n)),
            )
        end
    end

    # Nested commutator [H, [H, … σ⁻]] for growing depth.
    let
        hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a = Destroy(hjc, :a, 1)
        σge = Transition(hjc, :σ, 1, 2, 2)
        σee = Transition(hjc, :σ, 2, 2, 2)
        @variables ω_c ω_a g
        H_s = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')
        ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
        H_q = QA.normal_form(
            ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
                gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
        )
        nested_s(H, op, n) = n == 0 ? op : commutator(H, nested_s(H, op, n - 1))
        nested_q(H, op, n) = n == 0 ? op : QA.normal_form(QA.comm(H, nested_q(H, op, n - 1)))
        for d in 1:8
            sweep_pair!(
                SUITE, "Nested [H, σ] depth", d, "Scaling / large systems", "Nested [H, σ] depth=$d",
                () -> nested_s(H_s, σge, d), () -> nested_q(H_q, QA.σm(), d),
            )
        end
    end

    # Many-mode tight-binding chain: M coupled bosonic modes. QA uses
    # integer-indexed modes a(k) (distinct indices commute); SQA uses a product
    # of M Fock spaces.
    let
        @variables ωm Jm
        ωq = QA.Pr"ω"; Jq = QA.Pr"J"
        for M in (2, 4, 8, 16)
            spaces = reduce(⊗, (FockSpace(Symbol(:f, k)) for k in 1:M))
            as = [Destroy(spaces, Symbol(:a, k), k) for k in 1:M]
            buildH_s() =
                sum(ωm * as[k]' * as[k] for k in 1:M) +
                Jm * sum(as[k]' * as[k + 1] + as[k + 1]' * as[k] for k in 1:(M - 1))
            Hlazy_q =
                sum(ωq * QA.a'(k) * QA.a(k) for k in 1:M) +
                sum(Jq * (QA.a'(k) * QA.a(k + 1) + QA.a'(k + 1) * QA.a(k)) for k in 1:(M - 1))
            buildH_q() = QA.normal_form(
                sum(ωq * QA.a'(k) * QA.a(k) for k in 1:M) +
                    sum(Jq * (QA.a'(k) * QA.a(k + 1) + QA.a'(k + 1) * QA.a(k)) for k in 1:(M - 1)),
            )
            sweep_pair!(
                SUITE, "Many-mode chain: build H", M, "Scaling / large systems",
                "Many-mode build H M=$M", buildH_s, buildH_q,
            )
            H_s = buildH_s()
            sweep_pair!(
                SUITE, "Many-mode chain: H²", M, "Scaling / large systems",
                "Many-mode H² M=$M", () -> H_s * H_s, () -> QA.normal_form(Hlazy_q * Hlazy_q),
            )
        end
    end

    # Hⁿ term-blowup: Jaynes–Cummings and a finite (3-spin) Dicke.
    let
        hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
        a = Destroy(hjc, :a, 1)
        σge = Transition(hjc, :σ, 1, 2, 2)
        σee = Transition(hjc, :σ, 2, 2, 2)
        @variables ω_c ω_a g
        H_s = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')
        ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
        Hlazy_q = ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
            gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
        for n in 2:4
            # QA uses its eager workflow here (fair match to SQA's eager `*`).
            sweep_pair!(
                SUITE, "Jaynes–Cummings Hⁿ", n, "Scaling / large systems", "JC Hⁿ n=$n",
                () -> foldl(*, fill(H_s, n)), () -> qa_eager_prod(fill(Hlazy_q, n)),
            )
        end
    end

    let
        # Finite Dicke with 3 explicit spin-½ sites.
        hfd = FockSpace(:c) ⊗ PauliSpace(:s1) ⊗ PauliSpace(:s2) ⊗ PauliSpace(:s3)
        b = Destroy(hfd, :b, 1)
        @variables ωd λd
        sx(site) = Pauli(hfd, :σ, 1, site)
        sz(site) = Pauli(hfd, :σ, 3, site)
        Hd_s = ωd * b' * b + sum(ωd * sz(s) + λd * (b + b') * sx(s) for s in 2:4)

        ωq = QA.Pr"ω"; λq = QA.Pr"λ"
        Hd_lazy_q = ωq * QA.a'() * QA.a() +
            sum(ωq * QA.σz(k) + λq * (QA.a() + QA.a'()) * QA.σx(k) for k in 1:3)
        for n in 2:4
            # QA uses its eager workflow here (fair match to SQA's eager `*`).
            sweep_pair!(
                SUITE, "Dicke (3 spins) Hⁿ", n, "Scaling / large systems", "Dicke(3) Hⁿ n=$n",
                () -> foldl(*, fill(Hd_s, n)), () -> qa_eager_prod(fill(Hd_lazy_q, n)),
            )
        end
    end

    return SUITE
end

# ── auto_normal_form demo (eager QA workflow) ─────────────────────────────────
# QA's `auto_normal_form(true)` makes `*`/`+` eager, like SQA. This is a fairer
# match to SQA's model on a per-operation basis. Measured separately because it
# is a global mode toggle.
function auto_normal_form_demo(; seconds = 2.0)
    hf = FockSpace(:c)
    ha = NLevelSpace(:a, 2, 1)
    hjc = hf ⊗ ha
    a = Destroy(hjc, :a, 1)
    σge = Transition(hjc, :σ, 1, 2, 2)
    σee = Transition(hjc, :σ, 2, 2, 2)
    @variables ω_c ω_a g
    Hjc_s() = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')

    ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
    Hjc_auto_q() = ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
        gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())

    s = median(run((@benchmarkable $Hjc_s()); seconds = seconds))
    QA.auto_normal_form(true)
    q = try
        median(run((@benchmarkable $Hjc_auto_q()); seconds = seconds))
    finally
        QA.auto_normal_form(false)
    end
    return (sqa = s, qa_auto = q)
end

# QA's *default* lazy workflow blows up on high powers because the whole product
# is expanded before collapsing. This contrasts the lazy default against the
# eager workflow used in the head-to-head — an ergonomics point (SQA is eager by
# default, so users never hit this), not a raw-speed claim.
function workflow_blowup_demo(; n = 4, seconds = 2.0)
    ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
    HL = ωcq * QA.a'() * QA.a() + ωaq * QA.σp() * QA.σm() +
        gq * (QA.a'() * QA.σm() + QA.a() * QA.σp())
    lazy = median(run((@benchmarkable qa_lazy_pow($HL, $n)); seconds = seconds))
    eager = median(run((@benchmarkable qa_eager_prod(fill($HL, $n))); seconds = seconds))
    return (; n, lazy, eager)
end

# ── Report ────────────────────────────────────────────────────────────────────

fmt_ratio(t_sqa, t_qa) = t_qa == 0 ? "—" : string(round(t_qa / t_sqa; digits = 1), "×")

function markdown_table(med)
    io = IOBuffer()
    println(io, "| Benchmark | SQA | QuantumAlgebra | QA / SQA | SQA allocs | QA allocs |")
    println(io, "|---|---:|---:|---:|---:|---:|")
    current = ""
    for (group, name) in ORDER
        if group != current
            current = group
            println(io, "| **$(current)** | | | | | |")
        end
        if haskey(CAPPED, (group, name))
            println(io, "| ", name, " | _capped (>$(round(CAPPED[(group, name)]; digits = 1)) s/op)_ | | | | |")
            continue
        end
        s = med[group][name]["SQA"]
        q = med[group][name]["QA"]
        println(
            io,
            "| ", name,
            " | ", BenchmarkTools.prettytime(time(s)),
            " | ", BenchmarkTools.prettytime(time(q)),
            " | ", fmt_ratio(time(s), time(q)),
            " | ", allocs(s),
            " | ", allocs(q),
            " |",
        )
    end
    return String(take!(io))
end

function scaling_plot(med, path)
    sweeps = filter(s -> !isempty(s[2]), SWEEPS)
    isempty(sweeps) && return nothing
    ncol = 2
    nrow = cld(length(sweeps), ncol)
    fig = Makie.Figure(; size = (520 * ncol, 360 * nrow))
    for (idx, (label, pts)) in enumerate(sweeps)
        r, c = fldmod1(idx, ncol)
        ax = Makie.Axis(
            fig[r, c]; title = label, xlabel = "size", ylabel = "time (µs)",
            yscale = log10, xscale = log10,
        )
        xs = Float64[]; ys_s = Float64[]; ys_q = Float64[]
        for (x, group, name) in pts
            push!(xs, x)
            push!(ys_s, time(med[group][name]["SQA"]) / 1e3)
            push!(ys_q, time(med[group][name]["QA"]) / 1e3)
        end
        Makie.scatterlines!(ax, xs, ys_s; label = "SQA", color = :dodgerblue)
        Makie.scatterlines!(ax, xs, ys_q; label = "QuantumAlgebra", color = :darkorange)
        Makie.axislegend(ax; position = :lt, framevisible = false)
    end
    Makie.save(path, fig)
    return path
end

const SUITE = BenchmarkGroup()
println(stderr, "Building suite (runs one trial eval per benchmark)…"); flush(stderr)
const _tbuild = @elapsed benchmark_qa_comparison!(SUITE)
println(stderr, "Suite built in $(round(_tbuild; digits = 1)) s\n"); flush(stderr)

if abspath(PROGRAM_FILE) == @__FILE__
    println(stderr, "QA v$(pkgversion(QA)), SQA v$(pkgversion(SecondQuantizedAlgebra))\n"); flush(stderr)
    validate_equivalence(); flush(stdout)

    # NOTE: `tune!` is statistically the right thing for a proper benchmark, but
    # at this suite size it costs ~10 min. Commented out for fast iteration —
    # re-enable for publication-quality numbers.
    # t_tune = @elapsed tune!(SUITE)
    # println(stderr, "tune! done in $(round(t_tune; digits = 1)) s"); flush(stderr)
    t_run = @elapsed (results = run(SUITE; verbose = true, seconds = 2))
    println(stderr, "run done in $(round(t_run; digits = 1)) s"); flush(stderr)
    med = median(results)

    table = markdown_table(med)
    println("\n", table)

    auto = auto_normal_form_demo()
    println(
        "\nauto_normal_form(true) demo — JC build H:",
        "\n  SQA          ", BenchmarkTools.prettytime(time(auto.sqa)),
        "\n  QA (eager)   ", BenchmarkTools.prettytime(time(auto.qa_auto)),
    )

    bl = workflow_blowup_demo()
    println(
        "\nQA workflow on JC H^$(bl.n) (lazy default vs eager):",
        "\n  QA lazy      ", BenchmarkTools.prettytime(time(bl.lazy)), "  (", allocs(bl.lazy), " allocs)",
        "\n  QA eager     ", BenchmarkTools.prettytime(time(bl.eager)), "  (", allocs(bl.eager), " allocs)",
    )

    outdir = joinpath(@__DIR__, "results")
    mkpath(outdir)
    open(joinpath(outdir, "qa_comparison.md"), "w") do f
        println(f, "<!-- Generated by benchmark/quantumalgebra_comparison.jl -->")
        println(f, "<!-- SQA v$(pkgversion(SecondQuantizedAlgebra)), QA v$(pkgversion(QA)) -->\n")
        println(f, table)
    end

    plotpath = joinpath(outdir, "qa_comparison_scaling.png")
    scaling_plot(med, plotpath)
    println(stderr, "\nWrote $(joinpath(outdir, "qa_comparison.md")) and $plotpath")
end
