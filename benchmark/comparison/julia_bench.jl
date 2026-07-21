"""
    benchmark/comparison/julia_bench.jl

Cross-package comparison benchmarks for the two Julia packages:
SecondQuantizedAlgebra.jl (SQA, eager) and QuantumAlgebra.jl (QA, lazy).
Scenario definitions, canonical keys, fairness contract and output format
live in BENCHMARKS.md.

Usage:
    julia --project=benchmark benchmark/comparison/julia_bench.jl
"""

using BenchmarkTools
using SecondQuantizedAlgebra
import QuantumAlgebra as QA
using Symbolics: @variables
import JSON
using Dates: today

const CAP_SECONDS = 30.0
const RESULTS_DIR = joinpath(@__DIR__, "results")

# Eager fold for QA products of n > 2 factors (fairness contract, rule 2).
qa_eager_prod(factors) = foldl((acc, h) -> QA.normal_form(acc * h), factors)

# ── Known-answer equivalence checks (fairness contract, rule 3) ──────────────
function validate_equivalence()
    hf = FockSpace(:f)
    c = Destroy(hf, :c)
    ht = NLevelSpace(:a, 2, 1)
    σ⁻ = Transition(ht, :σ, 1, 2)
    σee = Transition(ht, :σ, 2, 2)
    σgg = Transition(ht, :σ, 1, 1)

    checks = [
        (
            "a a† = 1 + a†a",
            isequal(c * c', 1 + c' * c),
            QA.normal_form(QA.a() * QA.a'()) == QA.normal_form(QA.QuExpr(1) + QA.a'() * QA.a()),
        ),
        (
            "a a a† = a†a a + 2a",
            isequal(c * c * c', c' * c * c + 2 * c),
            QA.normal_form(QA.a() * QA.a() * QA.a'()) ==
                QA.normal_form(QA.a'() * QA.a() * QA.a() + 2 * QA.a()),
        ),
        (
            "σ⁻σ⁺ + σ⁺σ⁻ = 1",
            isequal(σgg + σee, σ⁻ * σ⁻' + σ⁻' * σ⁻),
            QA.normal_form(QA.σm() * QA.σp() + QA.σp() * QA.σm()) == QA.normal_form(QA.QuExpr(1)),
        ),
    ]

    allok = true
    for (desc, sqa_ok, qa_ok) in checks
        allok &= sqa_ok & qa_ok
        println(stderr, "  ", (sqa_ok & qa_ok) ? "✓" : "✗", "  SQA:$(sqa_ok) QA:$(qa_ok)  $desc")
    end
    return allok
end

# ── Scenario thunks ───────────────────────────────────────────────────────────
# Returns an ordered Vector of (key, sqa_thunk, qa_thunk).
function build_scenarios()
    scenarios = Tuple{String, Function, Function}[]

    # 1. Jaynes–Cummings family.
    hjc = FockSpace(:c) ⊗ NLevelSpace(:a, 2, 1)
    a = Destroy(hjc, :a, 1)
    σge = Transition(hjc, :σ, 1, 2, 2)
    σee = Transition(hjc, :σ, 2, 2, 2)
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
    H_q = Hjc_q()

    push!(scenarios, ("jc_build", Hjc_s, Hjc_q))
    push!(scenarios, ("jc_H2", () -> H_s * H_s, () -> QA.normal_form(Hjc_q_lazy * Hjc_q_lazy)))
    push!(
        scenarios,
        ("jc_heisenberg", () -> commutator(H_s, a), () -> QA.normal_form(QA.comm(H_q, QA.a()))),
    )

    nested_s(H, op, n) = n == 0 ? op : commutator(H, nested_s(H, op, n - 1))
    nested_q(H, op, n) = n == 0 ? op : QA.normal_form(QA.comm(H, nested_q(H, op, n - 1)))
    for d in 1:7
        push!(
            scenarios,
            (
                "jc_nested_d$d",
                () -> nested_s(H_s, σge, d),
                () -> nested_q(H_q, QA.σm(), d),
            ),
        )
    end

    for p in 2:7
        push!(
            scenarios,
            (
                "jc_pow_n$p",
                () -> foldl(*, fill(H_s, p)),
                () -> qa_eager_prod(fill(Hjc_q_lazy, p)),
            ),
        )
    end

    # 2. Bosonic normal ordering (a·a†)ⁿ.
    hf = FockSpace(:f)
    c = Destroy(hf, :c)
    for n in (2, 4, 6, 8, 10, 12, 14)
        push!(
            scenarios,
            (
                "fock_reorder_n$n",
                () -> foldl(*, fill(c * c', n)),
                () -> qa_eager_prod(fill(QA.a() * QA.a'(), n)),
            ),
        )
    end

    # 3. Bose–Hubbard chain. A build-cost sweep over the chain length M feeds the
    #    scaling figure (M = 8 reuses `bose_hubbard_build`); the fixed M = 8 pair
    #    (build + H²) feeds the docs table.
    @variables ω J U
    ωq = QA.Pr"ω"; Jq = QA.Pr"J"; Uq = QA.Pr"U"

    bh_ops_sqa(M) = [
        Destroy(reduce(⊗, (FockSpace(Symbol(:f, k)) for k in 1:M)), Symbol(:a, k), k)
            for k in 1:M
    ]
    bh_build_sqa(as, M) =
        sum(ω * as[k]' * as[k] for k in 1:M) +
        J * sum(as[k]' * as[k + 1] + as[k + 1]' * as[k] for k in 1:(M - 1)) +
        (U / 2) * sum(as[k]' * as[k]' * as[k] * as[k] for k in 1:M)
    bh_lazy_q(M) =
        sum(ωq * QA.a'(k) * QA.a(k) for k in 1:M) +
        sum(Jq * (QA.a'(k) * QA.a(k + 1) + QA.a'(k + 1) * QA.a(k)) for k in 1:(M - 1)) +
        sum(Uq / 2 * QA.a'(k) * QA.a'(k) * QA.a(k) * QA.a(k) for k in 1:M)

    for M in (2, 4, 16, 32)
        as = bh_ops_sqa(M)
        push!(
            scenarios,
            (
                "bh_chain_M$M",
                () -> bh_build_sqa(as, M),
                () -> QA.normal_form(bh_lazy_q(M)),
            ),
        )
    end

    M = 8
    as8 = bh_ops_sqa(M)
    Hbh_s = bh_build_sqa(as8, M)
    Hbh_q_lazy = bh_lazy_q(M)
    push!(
        scenarios,
        ("bose_hubbard_build", () -> bh_build_sqa(as8, M), () -> QA.normal_form(bh_lazy_q(M))),
    )
    push!(
        scenarios,
        (
            "bose_hubbard_H2",
            () -> Hbh_s * Hbh_s,
            () -> QA.normal_form(Hbh_q_lazy * Hbh_q_lazy),
        ),
    )

    # 4. Tavis–Cummings symbolic sum.
    @variables N
    htc = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
    atc = Destroy(htc, :a, 1)
    i = Index(htc, :i, N, NLevelSpace(:atom, 2, 1))
    j = Index(htc, :j, N, NLevelSpace(:atom, 2, 1))
    gi = IndexedVariable(:g, i)
    σi(α, β, idx) = IndexedOperator(Transition(htc, :σ, α, β, 2), idx)

    tc_build_s() = Σ(gi * (atc' * σi(1, 2, i) + atc * σi(2, 1, i)), i)
    Htc_s = tc_build_s()
    tc_build_q() = QA.normal_form(
        QA.∑(:i, QA.Pr"g_i" * (QA.a'() * QA.σm(:i) + QA.a() * QA.σp(:i)))
    )
    Htc_q = tc_build_q()

    push!(scenarios, ("tavis_cummings_build", tc_build_s, tc_build_q))
    push!(
        scenarios,
        (
            "tavis_cummings_comm",
            () -> commutator(Htc_s, σi(2, 2, j)),
            () -> QA.normal_form(QA.comm(Htc_q, QA.σp(:j) * QA.σm(:j))),
        ),
    )

    # 5. Mean-field expectation value.
    push!(scenarios, ("meanfield_average", () -> average(H_s), () -> QA.expval(H_q)))

    return scenarios
end

# ── Timing ────────────────────────────────────────────────────────────────────
# One warm-up call (compilation), then one timed call for the cap check.
_trial(thunk) = (thunk(); @elapsed thunk())

function run_side!(results, capped, key, thunk)
    t = _trial(thunk)
    if t > CAP_SECONDS
        capped[key] = t
        println(stderr, "  $key: capped ($(round(t; digits = 1)) s/op)")
        return nothing
    end
    b = @benchmarkable $thunk()
    # Up to a 10 s wall-clock budget, but stop early once enough samples are
    # collected (BenchmarkTools' default sample cap): fast scenarios finish in a
    # fraction of a second, expensive ones use the full budget for a tighter min.
    trial = run(b; seconds = 10, samples = 10_000)
    m = minimum(trial)
    results[key] = Dict("time_ns" => time(m), "allocs" => Int(allocs(m)))
    println(stderr, "  $key: $(BenchmarkTools.prettytime(time(m))) ($(allocs(m)) allocs)")
    return nothing
end

function write_json(path, package, version, results, capped)
    out = Dict(
        "meta" => Dict(
            "package" => package,
            "version" => string(version),
            "language" => "Julia",
            "runtime_version" => string(VERSION),
            "date" => string(today()),
        ),
        "results" => results,
        "capped" => capped,
    )
    open(path, "w") do io
        JSON.print(io, out, 2)
    end
    return path
end

function main()
    println(stderr, "Known-answer equivalence checks:")
    validate_equivalence() || error("equivalence checks failed; aborting")

    scenarios = build_scenarios()
    mkpath(RESULTS_DIR)

    for (package, pkgmod, side) in (
            ("SecondQuantizedAlgebra.jl", SecondQuantizedAlgebra, 2),
            ("QuantumAlgebra.jl", QA, 3),
        )
        println(stderr, "\n$package:")
        results = Dict{String, Any}()
        capped = Dict{String, Float64}()
        for sc in scenarios
            run_side!(results, capped, sc[1], sc[side])
        end
        name = package == "SecondQuantizedAlgebra.jl" ? "sqa" : "quantumalgebra"
        path = write_json(
            joinpath(RESULTS_DIR, "$name.json"), package, pkgversion(pkgmod), results, capped,
        )
        println(stderr, "wrote $path")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
