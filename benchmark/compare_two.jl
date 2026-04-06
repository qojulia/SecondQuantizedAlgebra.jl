"""
    benchmark/compare_two.jl

Compare two benchmark data files side by side.

Usage:
    julia --project=benchmark benchmark/compare_two.jl <baseline.json> <target.json>
    julia --project=benchmark benchmark/compare_two.jl          # interactive pick from data/
"""

using BenchmarkTools
import Dates
using Dates: DateTime, @dateformat_str

const DATADIR = joinpath(@__DIR__, "data")

# ── Helpers ──────────────────────────────────────────────────────────────────

function load_result(path::String)
    data = BenchmarkTools.load(path)
    result = first(data)
    m = match(r"(\d{4}-\d{2}-\d{2}_\d{6})\.json$", basename(path))
    ts = m !== nothing ? DateTime(m[1], dateformat"yyyy-mm-dd_HHMMSS") : DateTime(0)
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

function fmt_allocs(t)
    t === nothing && return "—"
    return string(allocs(t))
end

function fmt_memory(t)
    t === nothing && return "—"
    return BenchmarkTools.prettymemory(memory(t))
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

# ── File selection ───────────────────────────────────────────────────────────

function pick_files()
    jsons = sort(filter(f -> endswith(f, ".json"), readdir(DATADIR; join = true)))
    if length(jsons) < 2
        error("Need at least 2 benchmark files in $DATADIR")
    end

    println("Available benchmark files:")
    for (i, f) in enumerate(jsons)
        println("  [$i] $(basename(f))")
    end

    print("\nBaseline index [1]: ")
    b = readline()
    baseline_idx = isempty(b) ? 1 : parse(Int, b)

    print("Target index [$(length(jsons))]: ")
    t = readline()
    target_idx = isempty(t) ? length(jsons) : parse(Int, t)

    return jsons[baseline_idx], jsons[target_idx]
end

# ── Main comparison ──────────────────────────────────────────────────────────

function compare(baseline_path::String, target_path::String)
    base = load_result(baseline_path)
    tgt = load_result(target_path)

    all_ids = unique!(sort!(vcat(leaf_ids(base.results), leaf_ids(tgt.results))))

    base_label = Dates.format(base.timestamp, "yyyy-mm-dd HH:MM")
    tgt_label = Dates.format(tgt.timestamp, "yyyy-mm-dd HH:MM")

    name_w = max(45, maximum(length, all_ids; init = 0) + 2)
    time_w = 16
    change_w = 10
    alloc_w = 10
    mem_w = 12

    # Header
    println()
    println("="^(name_w + 2 * time_w + change_w + 2 * alloc_w + change_w + 2 * mem_w + 20))
    println("Baseline: ", basename(baseline_path), "  (", base_label, ")")
    println("Target:   ", basename(target_path), "  (", tgt_label, ")")
    println("="^(name_w + 2 * time_w + change_w + 2 * alloc_w + change_w + 2 * mem_w + 20))
    println()

    print(rpad_or_trim("Benchmark", name_w))
    print(" │ ", rpad_or_trim("Base time", time_w))
    print(" │ ", rpad_or_trim("Target time", time_w))
    print(" │ ", rpad_or_trim("Δ time", change_w))
    print(" │ ", rpad_or_trim("Base alloc", alloc_w))
    print(" │ ", rpad_or_trim("Tgt alloc", alloc_w))
    print(" │ ", rpad_or_trim("Δ alloc", change_w))
    print(" │ ", rpad_or_trim("Base mem", mem_w))
    print(" │ ", rpad_or_trim("Tgt mem", mem_w))
    println()
    println("─"^(name_w + 2 * time_w + change_w + 2 * alloc_w + change_w + 2 * mem_w + 20))

    n_faster = 0
    n_slower = 0
    n_same = 0

    for id in all_ids
        b = get_leaf(base.results, id)
        t = get_leaf(tgt.results, id)

        change_str = fmt_change(b, t)

        # Alloc change
        alloc_change = if b !== nothing && t !== nothing
            ab, at = allocs(b), allocs(t)
            ab == 0 ? "" : begin
                pct = (at / ab - 1) * 100
                abs(pct) < 1 ? "  ~same" : pct < 0 ? " $(round(Int, pct))%" : " +$(round(Int, pct))%"
            end
        else
            ""
        end

        # Count regressions/improvements
        if b !== nothing && t !== nothing
            pct = (time(t) / time(b) - 1) * 100
            if pct < -1
                n_faster += 1
            elseif pct > 1
                n_slower += 1
            else
                n_same += 1
            end
        end

        print(rpad_or_trim(id, name_w))
        print(" │ ", rpad_or_trim(fmt_time(b), time_w))
        print(" │ ", rpad_or_trim(fmt_time(t), time_w))
        print(" │ ", rpad_or_trim(change_str, change_w))
        print(" │ ", rpad_or_trim(fmt_allocs(b), alloc_w))
        print(" │ ", rpad_or_trim(fmt_allocs(t), alloc_w))
        print(" │ ", rpad_or_trim(alloc_change, change_w))
        print(" │ ", rpad_or_trim(fmt_memory(b), mem_w))
        print(" │ ", rpad_or_trim(fmt_memory(t), mem_w))
        println()
    end

    println()
    println("Summary: $n_faster faster, $n_slower slower, $n_same unchanged")
    println("Legend: Δ = change from baseline → target. Negative = faster/fewer.")
    println()
end

# ── Entry point ──────────────────────────────────────────────────────────────

if length(ARGS) >= 2
    compare(ARGS[1], ARGS[2])
elseif length(ARGS) == 1
    error("Provide two file paths, or none for interactive mode.")
else
    baseline, target = pick_files()
    compare(baseline, target)
end
