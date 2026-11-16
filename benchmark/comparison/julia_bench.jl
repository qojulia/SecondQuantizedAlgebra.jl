"""
    benchmark/comparison/julia_bench.jl

Cross-package comparison benchmarks for the two Julia packages:
SecondQuantizedAlgebra.jl (SQA, eager) and QuantumAlgebra.jl (QA, lazy).
Scenario definitions, canonical keys, fairness contract and output format
live in BENCHMARKS.md.

Usage:
    julia --project=benchmark benchmark/comparison/julia_bench.jl
"""

using BenchmarkTools: @benchmarkable, run, minimum, time, allocs, memory, prettytime
using SecondQuantizedAlgebra: SecondQuantizedAlgebra, FockSpace, Destroy, NLevelSpace, Transition, ⊗, @variables, commutator, Index, IndexedVariable, IndexedOperator, Σ, average, expand_completeness
import QuantumAlgebra as QA
import JSON
using Dates: today
using BenchmarkTools: median

QA.auto_normal_form(false)
QA.use_σpm(false)

const CAP_SECONDS = 30.0
const RESULTS_DIR = get(ENV, "SQA_COMPARE_OUTPUT", joinpath(@__DIR__, "results"))

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
            iszero(expand_completeness(σ⁻ * σ⁻' + σ⁻' * σ⁻) - 1),
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
    adag = a'
    σeg = σge'
    Hjc_s() = ω_c * adag * a + ω_a * σee + g * (adag * σge + a * σeg)
    H_s = Hjc_s()

    ωcq = QA.Pr"ω_c"; ωaq = QA.Pr"ω_a"; gq = QA.Pr"g"
    aq_jc = QA.a()
    adq_jc = QA.a'()
    smq_jc = QA.σm()
    spq_jc = QA.σp()
    eeq_jc = QA.normal_form(spq_jc * smq_jc)
    Hjc_q() = QA.normal_form(ωcq * adq_jc * aq_jc + ωaq * eeq_jc + gq * (adq_jc * smq_jc + aq_jc * spq_jc))
    H_q = Hjc_q()

    push!(scenarios, ("jc_build", Hjc_s, Hjc_q))
    push!(scenarios, ("jc_H2", () -> H_s * H_s, () -> QA.normal_form(H_q * H_q)))
    push!(
        scenarios,
        ("jc_heisenberg", () -> commutator(H_s, a), () -> QA.normal_form(QA.comm(H_q, aq_jc))),
    )

    nested_s(H, op, n) = n == 0 ? op : commutator(H, nested_s(H, op, n - 1))
    nested_q(H, op, n) = n == 0 ? op : QA.normal_form(QA.comm(H, nested_q(H, op, n - 1)))
    for d in 1:7
        push!(
            scenarios,
            (
                "jc_nested_d$d",
                () -> nested_s(H_s, σge, d),
                () -> nested_q(H_q, smq_jc, d),
            ),
        )
    end

    for p in 2:7
        push!(
            scenarios,
            (
                "jc_pow_n$p",
                () -> foldl(*, fill(H_s, p)),
                () -> qa_eager_prod(fill(H_q, p)),
            ),
        )
    end

    # 2. Bosonic normal ordering (a·a†)ⁿ.
    hf = FockSpace(:f)
    c = Destroy(hf, :c)
    factor_s = c * c'
    factor_q = QA.normal_form(QA.a() * QA.a'())
    for n in (2, 4, 6, 8, 10, 12, 14)
        push!(
            scenarios,
            (
                "fock_reorder_n$n",
                () -> foldl(*, fill(factor_s, n)),
                () -> qa_eager_prod(fill(factor_q, n)),
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
    bh_build_sqa(as, ads, M) =
        sum(ω * ads[k] * as[k] for k in 1:M) +
        J * sum(ads[k] * as[k + 1] + ads[k + 1] * as[k] for k in 1:(M - 1)) +
        (U / 2) * sum(ads[k] * ads[k] * as[k] * as[k] for k in 1:M)
    bh_ops_q(M) = [QA.a(k) for k in 1:M]
    bh_lazy_q(aq, adq, M) =
        sum(ωq * adq[k] * aq[k] for k in 1:M) +
        sum(Jq * (adq[k] * aq[k + 1] + adq[k + 1] * aq[k]) for k in 1:(M - 1)) +
        sum(Uq / 2 * adq[k] * adq[k] * aq[k] * aq[k] for k in 1:M)

    for M in (2, 4, 16, 32)
        as = bh_ops_sqa(M)
        ads = adjoint.(as)
        aqM = bh_ops_q(M)
        adqM = adjoint.(aqM)
        push!(
            scenarios,
            (
                "bh_chain_M$M",
                () -> bh_build_sqa(as, ads, M),
                () -> QA.normal_form(bh_lazy_q(aqM, adqM, M)),
            ),
        )
    end

    M = 8
    as8 = bh_ops_sqa(M)
    ads8 = adjoint.(as8)
    aq = bh_ops_q(M)
    adq = adjoint.(aq)
    Hbh_s = bh_build_sqa(as8, ads8, M)
    Hbh_q = QA.normal_form(bh_lazy_q(aq, adq, M))
    push!(
        scenarios,
        ("bose_hubbard_build", () -> bh_build_sqa(as8, ads8, M), () -> QA.normal_form(bh_lazy_q(aq, adq, M))),
    )
    push!(
        scenarios,
        (
            "bose_hubbard_H2",
            () -> Hbh_s * Hbh_s,
            () -> QA.normal_form(Hbh_q * Hbh_q),
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

    atc_dag = atc'
    sm_i = σi(1, 2, i)
    sp_i = σi(2, 1, i)
    population_s = σi(2, 2, j)
    sm_iq = QA.σm(:i)
    sp_iq = QA.σp(:i)
    giq = QA.Pr"g_i"
    population_q = QA.normal_form(QA.σp(:j) * QA.σm(:j))
    tc_build_s() = Σ(gi * (atc_dag * sm_i + atc * sp_i), i)
    Htc_s = tc_build_s()
    tc_build_q() = QA.normal_form(QA.∑(:i, giq * (adq_jc * sm_iq + aq_jc * sp_iq)))
    Htc_q = tc_build_q()

    push!(scenarios, ("tavis_cummings_build", tc_build_s, tc_build_q))
    push!(
        scenarios,
        (
            "tavis_cummings_comm",
            () -> commutator(Htc_s, population_s),
            () -> QA.normal_form(QA.comm(Htc_q, population_q)),
        ),
    )

    # 5. Mean-field expectation value.
    push!(scenarios, ("meanfield_average", () -> average(H_s), () -> QA.expval(H_q)))

    if get(ENV, "SQA_COMPARE_VALIDATE", "true") == "true"
        stirling(n, k) = k == 0 ? Int(n == 0) : n == 0 ? 0 : k * stirling(n - 1, k) + stirling(n - 1, k - 1)
        for n in (2, 4, 8, 14)
            expected_s = sum(stirling(n + 1, k + 1) * c'^k * c^k for k in 0:n)
            expected_q = QA.normal_form(sum(stirling(n + 1, k + 1) * QA.a'()^k * QA.a()^k for k in 0:n))
            @assert isequal(foldl(*, fill(c * c', n)), expected_s)
            @assert qa_eager_prod(fill(QA.normal_form(QA.a() * QA.a'()), n)) == expected_q
        end
        @assert isequal(commutator(H_s, a), -ω_c * a - g * σge)
        @assert QA.normal_form(QA.comm(H_q, QA.a())) == QA.normal_form(-ωcq * QA.a() - gq * QA.σm())
        for d in (2, 4)
            expected_s = sum((-1)^k * binomial(d, k) * H_s^(d - k) * σge * H_s^k for k in 0:d)
            power_q(n) = n == 0 ? QA.QuExpr(1) : qa_eager_prod(fill(H_q, n))
            expected_q = QA.normal_form(sum((-1)^k * binomial(d, k) * power_q(d - k) * QA.σm() * power_q(k) for k in 0:d))
            @assert nested_q(H_q, QA.σm(), d) == expected_q
            @assert isequal(nested_s(H_s, σge, d), expected_s)
        end
        @assert isequal(Hbh_s, Hbh_s')
        @assert Hbh_q == QA.normal_form(Hbh_q')
        expected_s = -ω * as8[1] - J * as8[2] - U * as8[1]' * as8[1]^2
        expected_q = QA.normal_form(-ωq * aq[1] - Jq * aq[2] - Uq * aq[1]' * aq[1]^2)
        @assert isequal(commutator(Hbh_s, as8[1]), expected_s)
        @assert QA.normal_form(QA.comm(Hbh_q, aq[1])) == expected_q
        gj = IndexedVariable(:g, j)
        @assert isequal(commutator(Htc_s, σi(2, 2, j)), gj * (atc' * σi(1, 2, j) - atc * σi(2, 1, j)))
        @assert QA.normal_form(QA.comm(Htc_q, QA.σp(:j) * QA.σm(:j))) == QA.normal_form(QA.Pr"g_j" * (QA.a'() * QA.σm(:j) - QA.a() * QA.σp(:j)))
        println(stderr, "Scenario checks passed: exact normal-order coefficients, JC drift/binomial commutators, Bose–Hubbard drift/Hermiticity, indexed TC contraction")
    end
    return scenarios
end

# ── Timing ────────────────────────────────────────────────────────────────────
# One warm-up call (compilation), then one timed call for the cap check.
trial(thunk) = (thunk(); @elapsed thunk())

function run_side!(results, capped, key, thunk)
    t = trial(thunk)
    if t > CAP_SECONDS
        capped[key] = t
        println(stderr, "  $key: capped ($(round(t; digits = 1)) s/op)")
        return nothing
    end
    b = @benchmarkable $thunk()
    # Up to a 10 s wall-clock budget, but stop early once enough samples are
    # collected (BenchmarkTools' default sample cap): fast scenarios finish in a
    # fraction of a second, expensive ones use the full budget for a tighter min.
    measurement = run(b; seconds = parse(Float64, get(ENV, "SQA_COMPARE_SECONDS", "10")), samples = parse(Int, get(ENV, "SQA_COMPARE_SAMPLES", "10000")))
    m = minimum(measurement)
    results[key] = Dict("time_ns" => time(m), "median_ns" => median(measurement.times), "samples_ns" => measurement.times, "allocs" => Int(allocs(m)), "bytes" => memory(m))
    println(stderr, "  $key: $(prettytime(time(m))) ($(allocs(m)) allocs)")
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
            "contract" => "canonical-symbolic-v2",
            "run_id" => get(ENV, "SQA_COMPARE_RUN_ID", "standalone"),
            "coefficient_domain" => "symbolic",
            "julia_threads" => Threads.nthreads(),
            "auto_normal_form" => false,
            "spin_basis" => "xyz",
            "validation" => "passed",
        ),
        "results" => results,
        "capped" => capped,
    )
    open(path * ".tmp", "w") do io
        JSON.print(io, out, 2)
    end
    mv(path * ".tmp", path; force = true)
    return path
end

function main()
    println(stderr, "Known-answer equivalence checks:")
    validate_equivalence() || error("equivalence checks failed; aborting")

    scenarios = build_scenarios()
    get(ENV, "SQA_COMPARE_VALIDATE_ONLY", "false") == "true" && return nothing
    selected = split(get(ENV, "SQA_COMPARE_KEYS", ""), ",")
    isempty(first(selected)) || filter!(sc -> sc[1] in selected, scenarios)
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
