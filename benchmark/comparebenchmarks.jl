using BenchmarkTools: BenchmarkGroup, TrialEstimate, load

function benchmark_times!(times::Dict{String, Float64}, group::BenchmarkGroup, prefix::String)
    for (key, value) in group
        name = isempty(prefix) ? string(key) : string(prefix, "/", key)
        if value isa BenchmarkGroup
            benchmark_times!(times, value, name)
        else
            times[name] = (value::TrialEstimate).time
        end
    end
    return times
end

function compare_benchmarks(base_file::String, head_file::String)
    base = benchmark_times!(Dict{String, Float64}(), only(load(base_file)), "")
    head = benchmark_times!(Dict{String, Float64}(), only(load(head_file)), "")
    common = sort!(collect(intersect(keys(base), keys(head))))
    isempty(common) && error("No matching benchmarks in the base and PR measurements")
    failures = String[]
    for name in common
        ratio = head[name] / base[name]
        isfinite(ratio) || error("Invalid benchmark ratio for $name")
        if ratio > 1.3
            println(name, ": ", round(ratio; digits = 3), "x (PR/base)")
        end
        ratio > 1.7 && push!(failures, name)
    end
    println("Compared ", length(common), " benchmarks on the same runner")
    isempty(failures) || error("Benchmarks exceeded the 1.70x failure threshold: " * join(failures, ", "))
    return nothing
end

length(ARGS) == 2 || error("Usage: comparebenchmarks.jl BASE.json HEAD.json")
compare_benchmarks(ARGS[1], ARGS[2])
