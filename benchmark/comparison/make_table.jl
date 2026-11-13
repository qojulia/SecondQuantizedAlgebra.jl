"""
    benchmark/comparison/make_table.jl

Render one validated comparison run into tables and figures. Use `run.py report`
to check provenance and coverage before rendering; publication replaces the
marked documentation results block from the same dataset.
"""

import JSON
import CairoMakie as Makie
using BenchmarkTools: prettytime, median

const RESULTS_DIR = get(ENV, "SQA_COMPARE_OUTPUT", joinpath(@__DIR__, "results"))
const REPORT_DIR = get(ENV, "SQA_COMPARE_REPORT", RESULTS_DIR)
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
    ("meanfield_average", "Interface construction: ⟨H⟩ (different representations)"),
]

function load_results()
    manifest_path = joinpath(RESULTS_DIR, "run.json")
    isfile(manifest_path) || error("Use run.py to create a validated run archive; legacy results cannot be republished")
    manifest = JSON.parsefile(manifest_path)
    all(==("passed"), values(manifest["jobs"])) || error("Run has incomplete jobs")
    data = Dict{String, Any}()
    for name in manifest["packages"]
        path = joinpath(RESULTS_DIR, "$name.json")
        isfile(path) || error("Missing package results: $name")
        result = JSON.parsefile(path)
        result["meta"]["run_id"] == manifest["run_id"] || error("Mixed run ids")
        result["meta"]["contract"] == "canonical-symbolic-v2" || error("Wrong contract")
        result["meta"]["contract_sha256"] == manifest["contract_sha256"] || error("Mixed contracts")
        data[name] = result
    end
    return data
end

time_ns(data, pkg, key) =
    haskey(data, pkg) && haskey(data[pkg]["results"], key) ?
    Float64(data[pkg]["results"][key]["time_ns"]) : nothing

is_capped(data, pkg, key) = haskey(data, pkg) && haskey(get(data[pkg], "capped", Dict()), key)

function cell(data, pkg, key, baseline)
    is_capped(data, pkg, key) && return "capped (>30 s)"
    t = time_ns(data, pkg, key)
    t === nothing && return replace(get(get(data[pkg], "statuses", Dict()), key, "not_run"), "_" => " ")
    (pkg == "sqa" || key == "meanfield_average") && return prettytime(t)
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
        ("Nested commutator [H, [H, … σ⁻]]", "depth", [(d, "jc_nested_d$d") for d in 1:7]),
        ("Operator power Hⁿ (Jaynes–Cummings)", "n", [(p, "jc_pow_n$p") for p in 2:7]),
        ("Normal ordering (a·a†)ⁿ", "n", [(n, "fock_reorder_n$n") for n in 2:2:14]),
        (
            "Bose–Hubbard chain: build H (M modes)", "M",
            [(2, "bh_chain_M2"), (4, "bh_chain_M4"), (8, "bose_hubbard_build"), (16, "bh_chain_M16"), (32, "bh_chain_M32")],
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
                push!(xs, x); push!(ys, t === nothing ? NaN : t)
            end
            any(isfinite, ys) || continue
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

function overview_figure(data, path)
    rows = filter(row -> row[1] != "meanfield_average", ROWS)
    fig = Makie.Figure(; size = (1100, 600), backgroundcolor = :white)
    ax = Makie.Axis(
        fig[1, 1]; xscale = log10,
        xlabel = "relative time (package / SQA); below 1 means package is faster",
        yticks = (1:length(rows), last.(rows)), yreversed = true
    )
    Makie.vlines!(ax, [1.0]; color = :gray, linestyle = :dash)
    for (name, disp) in PACKAGES
        name == "sqa" && continue
        for (i, (key, _)) in enumerate(rows)
            baseline = time_ns(data, "sqa", key)
            value = time_ns(data, name, key)
            (baseline === nothing || value === nothing) && continue
            Makie.scatter!(
                ax, [value / baseline], [i]; label = disp,
                color = COLORS[name], marker = MARKERS[name], markersize = 10
            )
        end
    end
    Makie.axislegend(ax; unique = true, position = :rb)
    Makie.save(path, fig)
    return path
end

function main()
    data = load_results()
    table = markdown_table(data)
    footer = versions_line(data)

    mkpath(REPORT_DIR)
    tablepath = joinpath(REPORT_DIR, "comparison_table.md")
    open(tablepath, "w") do io
        println(io, "<!-- Generated by benchmark/comparison/make_table.jl; do not edit by hand. -->")
        println(io, "<!-- $footer -->\n")
        println(io, table)
    end
    println(table)
    println("Versions: ", footer)

    figpath = scaling_figure(data, joinpath(REPORT_DIR, "comparison_scaling.png"))
    scaling_figure(data, joinpath(REPORT_DIR, "comparison_scaling.svg"))
    scaling_figure(data, joinpath(REPORT_DIR, "comparison_scaling.pdf"))
    overview_figure(data, joinpath(REPORT_DIR, "comparison_overview.svg"))
    overview_figure(data, joinpath(REPORT_DIR, "comparison_overview.png"))
    if get(ENV, "SQA_COMPARE_PUBLISH", "false") == "true"
        manifest = JSON.parsefile(joinpath(RESULTS_DIR, "run.json"))
        length(manifest["keys"]) == 32 || error("Cannot publish a smoke run")
        docs = joinpath(@__DIR__, "..", "..", "docs", "src", "comparison.md")
        source = read(docs, String)
        generated = "<!-- comparison-results:start -->\n" *
            "Run `$(manifest["run_id"])`; contract `$(manifest["contract"])`; minimum warm runtime.\n\n" *
            footer * ".\n\n" * table * "\n<!-- comparison-results:end -->"
        occursin("<!-- comparison-results:start -->", source) || error("Missing results markers")
        source = replace(source, r"<!-- comparison-results:start -->.*?<!-- comparison-results:end -->"s => generated)
        write(docs, source)
        for name in manifest["packages"]
            summary = JSON.parsefile(joinpath(RESULTS_DIR, "$name.json"))
            summary["meta"]["raw_samples_archive"] = "reproduction.tar.gz"
            for result in values(summary["results"])
                if haskey(result, "samples_ns")
                    samples = pop!(result, "samples_ns")
                    result["sample_count"] = length(samples)
                    result["median_ns"] = median(samples)
                end
            end
            open(joinpath(@__DIR__, "results", "$name.json"), "w") do io
                JSON.print(io, summary, 2)
            end
        end
        cp(tablepath, joinpath(@__DIR__, "results", "comparison_table.md"); force = true)
        cp(figpath, joinpath(DOCS_ASSETS, "comparison_scaling.png"); force = true)
        cp(joinpath(REPORT_DIR, "comparison_overview.svg"), joinpath(DOCS_ASSETS, "comparison_overview.svg"); force = true)
    end
    println(stderr, "wrote $tablepath and $figpath")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
