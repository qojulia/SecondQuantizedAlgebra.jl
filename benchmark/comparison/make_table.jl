"""
    benchmark/comparison/make_table.jl

Merge the per-package `results/<package>.json` files (see BENCHMARKS.md) into
the markdown results table and the scaling figure used by
`docs/src/comparison.md`. Packages whose JSON is missing are skipped, so the
table can be regenerated from partial results.

Usage:
    julia --project=benchmark benchmark/comparison/make_table.jl
"""

import JSON
import CairoMakie as Makie
using BenchmarkTools: prettytime

const RESULTS_DIR = joinpath(@__DIR__, "results")
const DOCS_ASSETS = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets"))

# Fixed package order and display names; SQA is the ratio baseline.
const PACKAGES = [
    ("sqa", "SQA"),
    ("quantumalgebra", "QuantumAlgebra.jl"),
    ("sympy", "SymPy"),
    ("openfermion", "OpenFermion"),
    ("sneg", "sneg"),
]

# Canonical keys in display order (see BENCHMARKS.md).
const ROWS = [
    ("jc_build", "Jaynes–Cummings: build H"),
    ("jc_H2", "Jaynes–Cummings: H²"),
    ("jc_heisenberg", "Heisenberg equation: [H, a]"),
    ("jc_nested_d2", "Nested [H, [H, … σ⁻]] depth 2"),
    ("jc_nested_d4", "Nested [H, [H, … σ⁻]] depth 4"),
    ("jc_nested_d6", "Nested [H, [H, … σ⁻]] depth 6"),
    ("fock_reorder_n4", "Normal-order (a·a†)⁴"),
    ("fock_reorder_n8", "Normal-order (a·a†)⁸"),
    ("bose_hubbard_build", "Bose–Hubbard (M = 8): build H"),
    ("bose_hubbard_H2", "Bose–Hubbard (M = 8): H²"),
    ("tavis_cummings_build", "Tavis–Cummings Σᵢ: build H"),
    ("tavis_cummings_comm", "Tavis–Cummings: [H, σ⁺ⱼσ⁻ⱼ]"),
    ("meanfield_average", "Mean-field ⟨H⟩"),
]

function load_results()
    data = Dict{String, Any}()
    for (name, _) in PACKAGES
        path = joinpath(RESULTS_DIR, "$name.json")
        isfile(path) || continue
        data[name] = JSON.parsefile(path)
    end
    isempty(data) && error("no result files found in $RESULTS_DIR")
    return data
end

time_ns(data, pkg, key) =
    haskey(data, pkg) && haskey(data[pkg]["results"], key) ?
    Float64(data[pkg]["results"][key]["time_ns"]) : nothing

is_capped(data, pkg, key) = haskey(data, pkg) && haskey(get(data[pkg], "capped", Dict()), key)

function cell(data, pkg, key, baseline)
    is_capped(data, pkg, key) && return "capped (>3 s)"
    t = time_ns(data, pkg, key)
    t === nothing && return "n/a"
    pkg == "sqa" && return prettytime(t)
    baseline === nothing && return prettytime(t)
    r = round(t / baseline; sigdigits = 2)
    return prettytime(t) * " ($(isinteger(r) ? string(Int(r)) : string(r))×)"
end

function markdown_table(data)
    io = IOBuffer()
    present = [(name, disp) for (name, disp) in PACKAGES if haskey(data, name)]
    println(io, "| Benchmark | ", join((disp for (_, disp) in present), " | "), " |")
    println(io, "|---|", repeat("---:|", length(present)))
    for (key, disp) in ROWS
        baseline = time_ns(data, "sqa", key)
        cells = [cell(data, name, key, baseline) for (name, _) in present]
        println(io, "| ", disp, " | ", join(cells, " | "), " |")
    end
    return String(take!(io))
end

function versions_line(data)
    parts = String[]
    for (name, disp) in PACKAGES
        haskey(data, name) || continue
        m = data[name]["meta"]
        push!(parts, "$disp $(m["version"]) ($(m["language"]) $(m["runtime_version"]))")
    end
    return join(parts, ", ")
end

# ── Scaling figure ────────────────────────────────────────────────────────────
# Validated categorical palette (light mode), fixed package order; marker
# shapes are the secondary encoding.
const COLORS = Dict(
    "sqa" => "#2a78d6",
    "quantumalgebra" => "#008300",
    "sympy" => "#e87ba4",
    "openfermion" => "#eda100",
    "sneg" => "#1baf7a",
)
const MARKERS = Dict(
    "sqa" => :circle,
    "quantumalgebra" => :rect,
    "sympy" => :utriangle,
    "openfermion" => :diamond,
    "sneg" => :cross,
)

function compact_time(ns)
    ns < 1.0e3 && return "$(round(Int, ns)) ns"
    ns < 1.0e6 && return "$(round(Int, ns / 1.0e3)) μs"
    ns < 1.0e9 && return "$(round(Int, ns / 1.0e6)) ms"
    return "$(round(ns / 1.0e9; digits = 1)) s"
end

function scaling_figure(data, path)
    panels = [
        ("Nested commutator [H, [H, … σ⁻]]", "depth", [(d, "jc_nested_d$d") for d in 1:6]),
        ("Operator power Hⁿ (Jaynes–Cummings)", "n", [(p, "jc_pow_n$p") for p in 2:6]),
        ("Normal ordering (a·a†)ⁿ", "n", [(n, "fock_reorder_n$n") for n in 2:2:10]),
        (
            "Bose–Hubbard chain: build H (M modes)", "M",
            [(2, "bh_chain_M2"), (4, "bh_chain_M4"), (8, "bose_hubbard_build"), (16, "bh_chain_M16")],
        ),
    ]
    fig = Makie.Figure(; size = (1150, 800), backgroundcolor = :white)
    for (idx, (title, xlabel, points)) in enumerate(panels)
        row = (idx - 1) ÷ 2 + 1
        col = (idx - 1) % 2 + 1
        ax = Makie.Axis(
            fig[row, col];
            title, xlabel, ylabel = col == 1 ? "time per operation" : "",
            yscale = log10,
            xticks = first.(points),
            ytickformat = ts -> [compact_time(t) for t in ts],
            xgridvisible = false,
            ygridcolor = (:black, 0.08),
            titlealign = :left,
        )
        for (name, disp) in PACKAGES
            xs = Float64[]; ys = Float64[]
            for (x, key) in points
                t = time_ns(data, name, key)
                t === nothing && continue
                push!(xs, x); push!(ys, t)
            end
            isempty(xs) && continue
            Makie.scatterlines!(
                ax, xs, ys;
                label = disp, color = COLORS[name], marker = MARKERS[name],
                linewidth = 2, markersize = 12,
            )
        end
    end
    present = [(name, disp) for (name, disp) in PACKAGES if haskey(data, name)]
    elements = [
        [
                Makie.LineElement(; color = COLORS[name], linewidth = 2),
                Makie.MarkerElement(; color = COLORS[name], marker = MARKERS[name], markersize = 12),
            ] for (name, _) in present
    ]
    fig[1:2, 3] = Makie.Legend(fig, elements, [disp for (_, disp) in present]; framevisible = false)
    Makie.save(path, fig)
    return path
end

function main()
    data = load_results()
    table = markdown_table(data)
    footer = versions_line(data)

    mkpath(RESULTS_DIR)
    tablepath = joinpath(RESULTS_DIR, "comparison_table.md")
    open(tablepath, "w") do io
        println(io, "<!-- Generated by benchmark/comparison/make_table.jl; do not edit by hand. -->")
        println(io, "<!-- $footer -->\n")
        println(io, table)
    end
    println(table)
    println("Versions: ", footer)

    mkpath(DOCS_ASSETS)
    figpath = scaling_figure(data, joinpath(DOCS_ASSETS, "comparison_scaling.png"))
    println(stderr, "wrote $tablepath and $figpath")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
