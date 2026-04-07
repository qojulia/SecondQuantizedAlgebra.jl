"""
    benchmark/eager_vs_lazy.jl

Compare eager normal ordering (NormalOrder during multiplication) vs lazy ordering
(LazyOrder during multiplication, then normal_order() applied afterwards).

Benchmarks increasingly complex systems and produces a grouped bar plot.

Usage:
    julia --project=benchmark benchmark/eager_vs_lazy.jl
"""

using SecondQuantizedAlgebra
using Symbolics: @variables
using BenchmarkTools

# ── Helpers ─────────────────────────────────────────────────────────────────

"""
    bench_eager(f) -> BenchmarkTools.Trial

Run `f()` with `NormalOrder()` active (eager ordering during multiplication).
"""
function bench_eager(f)
    prev = get_ordering()
    set_ordering!(NormalOrder())
    try
        return @benchmark $f()
    finally
        set_ordering!(prev)
    end
end

"""
    bench_lazy(f) -> BenchmarkTools.Trial

Run `f()` with `LazyOrder()` active, then apply `normal_order` to the result.
"""
function bench_lazy(f)
    prev = get_ordering()
    set_ordering!(LazyOrder())
    try
        return @benchmark normal_order($f())
    finally
        set_ordering!(prev)
    end
end

# ── System constructors ─────────────────────────────────────────────────────
# Each returns a closure that builds the expression from scratch,
# so the benchmark measures construction + ordering.

function make_fock_power(n::Int)
    hf = FockSpace(:f)
    @qnumbers c::Destroy(hf)
    return () -> (c * c')^n
end

function make_jaynes_cummings_Hn(n::Int)
    hf = FockSpace(:c)
    ha = NLevelSpace(:a, 2, 1)
    h = hf ⊗ ha
    a = Destroy(h, :a, 1)
    σge = Transition(h, :σ, 1, 2, 2)
    σee = Transition(h, :σ, 2, 2, 2)
    @variables ω_c ω_a g
    H = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')
    return () -> H^n
end

function make_multimode_fock(n_modes::Int)
    spaces = [FockSpace(Symbol(:m, i)) for i in 1:n_modes]
    h = reduce(⊗, spaces)
    ops = [Destroy(h, Symbol(:a, i), i) for i in 1:n_modes]
    @variables J
    # Hopping chain: sum_i J * (a_i† * a_{i+1} + h.c.)
    return function ()
        H = J * ops[1]' * ops[1]
        for i in 1:(n_modes - 1)
            H += J * (ops[i]' * ops[i + 1] + ops[i + 1]' * ops[i])
        end
        return H * H
    end
end

function make_lambda_system_Hn(n::Int)
    h = FockSpace(:c) ⊗ NLevelSpace(:Λ, 3, 1)
    a = Destroy(h, :a, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Ω_p Ω_c Δ_p Δ_c
    H = Δ_p * σ(2, 2) + (Δ_p - Δ_c) * σ(3, 3) +
        Ω_p * (a' * σ(1, 2) + a * σ(2, 1)) +
        Ω_c * (σ(2, 3) + σ(3, 2))
    return () -> H^n
end

function make_nlevel_chain(N::Int)
    h = NLevelSpace(:atom, N, 1)
    t(i, j) = Transition(h, :σ, i, j)
    # Build a nearest-neighbour coupling Hamiltonian
    @variables Ω
    return function ()
        H = Ω * t(1, 1) # diagonal seed
        for i in 1:(N - 1)
            H += Ω * (t(i, i + 1) + t(i + 1, i))
        end
        return H * H
    end
end

# ── Benchmark suite ─────────────────────────────────────────────────────────

struct BenchCase
    label::String
    f::Any  # closure
end

cases = BenchCase[
    # Fock power scaling
    BenchCase("Fock (cc†)²", make_fock_power(2)),
    BenchCase("Fock (cc†)³", make_fock_power(3)),
    BenchCase("Fock (cc†)⁴", make_fock_power(4)),
    BenchCase("Fock (cc†)⁵", make_fock_power(5)),
    # Jaynes-Cummings power scaling
    BenchCase("JC H²", make_jaynes_cummings_Hn(2)),
    BenchCase("JC H³", make_jaynes_cummings_Hn(3)),
    # Multi-mode hopping H²
    BenchCase("2-mode H²", make_multimode_fock(2)),
    BenchCase("3-mode H²", make_multimode_fock(3)),
    BenchCase("4-mode H²", make_multimode_fock(4)),
    # Λ-system power scaling
    BenchCase("Λ H²", make_lambda_system_Hn(2)),
    BenchCase("Λ H³", make_lambda_system_Hn(3)),
    # N-level chain H²
    BenchCase("3-level H²", make_nlevel_chain(3)),
    BenchCase("4-level H²", make_nlevel_chain(4)),
    BenchCase("5-level H²", make_nlevel_chain(5)),
]

# ── Run benchmarks ──────────────────────────────────────────────────────────

println("Running eager vs lazy benchmarks ($(length(cases)) cases)...")
println("="^60)

labels = String[]
eager_times = Float64[]
lazy_times = Float64[]

for (i, c) in enumerate(cases)
    print("[$i/$(length(cases))] $(c.label) ... ")

    t_eager = bench_eager(c.f)
    t_lazy = bench_lazy(c.f)

    te = median(t_eager).time / 1.0e6  # ms
    tl = median(t_lazy).time / 1.0e6   # ms

    push!(labels, c.label)
    push!(eager_times, te)
    push!(lazy_times, tl)

    ratio = tl / te
    println("eager=$(round(te; digits = 3))ms  lazy=$(round(tl; digits = 3))ms  ratio=$(round(ratio; digits = 2))×")
end

# ── Print table ─────────────────────────────────────────────────────────────

println("\n", "="^70)
println(rpad("System", 20), rpad("Eager (ms)", 15), rpad("Lazy (ms)", 15), "Ratio (lazy/eager)")
println("-"^70)
for i in eachindex(labels)
    r = lazy_times[i] / eager_times[i]
    println(
        rpad(labels[i], 20),
        rpad(round(eager_times[i]; digits = 3), 15),
        rpad(round(lazy_times[i]; digits = 3), 15),
        round(r; digits = 2), "×",
    )
end
println("="^70)

# ── Plot ────────────────────────────────────────────────────────────────────

try
    using CairoMakie

    n = length(labels)
    ratios = lazy_times ./ eager_times
    xs = collect(1:n)

    fig = Figure(; size = (1000, 500), fontsize = 13)

    ax = Axis(
        fig[1, 1];
        title = "Eager vs Lazy Normal Ordering",
        ylabel = "Time (ms)",
        xticks = (1:n, labels),
        xticklabelrotation = π / 4,
        yscale = log10,
    )

    scatterlines!(ax, xs, eager_times;
        label = "Eager (NormalOrder)",
        color = :steelblue, marker = :circle, markersize = 10, linewidth = 2,
    )
    scatterlines!(ax, xs, lazy_times;
        label = "Lazy + normal_order()",
        color = :coral, marker = :utriangle, markersize = 10, linewidth = 2,
    )

    # Annotate ratios
    for i in 1:n
        ymax = max(eager_times[i], lazy_times[i])
        text!(
            ax, i, ymax;
            text = "$(round(ratios[i]; digits = 2))×",
            align = (:center, :bottom),
            fontsize = 9,
            color = ratios[i] > 1.1 ? :red : ratios[i] < 0.95 ? :green : :gray50,
        )
    end

    axislegend(ax; position = :lt)

    svgpath = joinpath(@__DIR__, "..", "docs", "src", "assets", "eager_vs_lazy.svg")
    mkpath(dirname(svgpath))
    save(svgpath, fig)
    println("\nPlot saved to docs/src/assets/eager_vs_lazy.svg")
catch e
    if e isa ArgumentError || e isa LoadError
        println("\nCairoMakie not available — skipping plot. Install with:")
        println("  ] add CairoMakie")
    else
        println("\nPlot failed: ", e)
        rethrow(e)
    end
end
