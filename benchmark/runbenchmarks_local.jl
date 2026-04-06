"""
    benchmark/runbenchmarks_local.jl

Run benchmarks, save results to `benchmark/data/`, and print a changelog
comparing all previous runs on the same CPU.

Usage:
    julia --project=benchmark benchmark/runbenchmarks_local.jl
"""

using BenchmarkTools
using SecondQuantizedAlgebra
import Dates
using Dates: now, DateTime, @dateformat_str

const SUITE = BenchmarkGroup()

include("commutator.jl")
include("simplify_and_normal_order.jl")
include("indexing.jl")

benchmark_commutator!(SUITE)
benchmark_simplify_and_normal_order!(SUITE)
benchmark_indexing!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
medians = median(results)

# ── Save results ─────────────────────────────────────────────────────────────

datadir = joinpath(@__DIR__, "data")
mkpath(datadir)

cpu_name = replace(Sys.cpu_info()[1].model, r"[^a-zA-Z0-9]" => "_")
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
filename = joinpath(datadir, "benchmark_$(cpu_name)_$(timestamp).json")

BenchmarkTools.save(filename, medians)

println("\nResults saved to: ", filename)

# ── Prune old results (keep last 5 per CPU) ──────────────────────────────────

let
    all_jsons = sort(
        filter(readdir(datadir; join = true)) do f
            endswith(f, ".json") && contains(basename(f), cpu_name)
        end
    )
    if length(all_jsons) > 5
        for old in all_jsons[1:(end - 5)]
            rm(old)
            println("Pruned old result: ", basename(old))
        end
    end
end

# ── Changelog comparison ─────────────────────────────────────────────────────

function load_result(path::String)
    data = BenchmarkTools.load(path)
    result = first(data)
    m = match(r"(\d{4}-\d{2}-\d{2}_\d{6})\.json$", basename(path))
    ts = if m !== nothing
        DateTime(m[1], dateformat"yyyy-mm-dd_HHMMSS")
    else
        DateTime(0)
    end
    return (; timestamp = ts, results = result, path = path)
end

function leaf_ids(bg::BenchmarkGroup; prefix::String = "")
    ids = String[]
    for (k, v) in bg
        full = isempty(prefix) ? string(k) : string(prefix, "/", k)
        if v isa BenchmarkGroup
            append!(ids, leaf_ids(v; prefix = full))
        else
            push!(ids, full)
        end
    end
    return sort(ids)
end

function get_leaf(bg::BenchmarkGroup, id::String)
    parts = split(id, "/")
    node = bg
    for p in parts
        haskey(node, p) || return nothing
        node = node[p]
    end
    return node
end

function fmt_time(t)
    t === nothing && return "—"
    return BenchmarkTools.prettytime(time(t))
end

function fmt_change(old, new)
    (old === nothing || new === nothing) && return ""
    t_old = time(old)
    t_new = time(new)
    t_old == 0 && return ""
    pct = (t_new / t_old - 1) * 100
    if abs(pct) < 1
        return "  ~same"
    elseif pct < 0
        return " $(round(pct; digits = 1))%"
    else
        return " +$(round(pct; digits = 1))%"
    end
end

function rpad_or_trim(s::AbstractString, n::Int)
    return length(s) >= n ? s[1:min(n, end)] : rpad(s, n)
end

function print_changelog(datadir::String, cpu_name::String)
    json_files = filter(readdir(datadir; join = true)) do f
        endswith(f, ".json") && contains(basename(f), cpu_name)
    end

    length(json_files) < 2 && (println("\nOnly one run recorded — changelog available after the next run."); return)

    entries = sort([load_result(f) for f in json_files]; by = e -> e.timestamp)
    entries = entries[max(1, end - 2):end]  # keep last 3 runs

    all_ids = String[]
    for e in entries
        append!(all_ids, leaf_ids(e.results))
    end
    unique!(sort!(all_ids))

    n = length(entries)

    println("\n", "="^70)
    println("CHANGELOG — ", Sys.cpu_info()[1].model, " (", n, " runs)")
    println("="^70)

    name_w = max(30, maximum(length, all_ids; init = 0) + 2)
    col_w = 22

    # Header
    print(rpad_or_trim("Benchmark", name_w))
    for (i, e) in enumerate(entries)
        label = Dates.format(e.timestamp, "yyyy-mm-dd HH:MM")
        print(" | ", rpad_or_trim(label, col_w))
        if i > 1
            print(rpad_or_trim("  Δ", 10))
        end
    end
    println()
    println("—"^(name_w + n * (col_w + 3) + (n - 1) * 10))

    # Rows
    for id in all_ids
        print(rpad_or_trim(id, name_w))
        prev = nothing
        for (i, e) in enumerate(entries)
            val = get_leaf(e.results, id)
            print(" | ", rpad_or_trim(fmt_time(val), col_w))
            if i > 1
                print(rpad_or_trim(fmt_change(prev, val), 10))
            end
            prev = val
        end
        println()
    end

    println()
    return println("Legend: Δ = change from previous run. Negative = faster.")
end

print_changelog(datadir, cpu_name)
