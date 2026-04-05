using BenchmarkTools
using SecondQuantizedAlgebra

const SUITE = BenchmarkGroup()

include("commutator.jl")
include("simplify_and_normal_order.jl")
include("indexing.jl")

benchmark_commutator!(SUITE)
benchmark_simplify_and_normal_order!(SUITE)
benchmark_indexing!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
