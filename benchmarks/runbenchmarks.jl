using BenchmarkTools
using SecondQuantizedAlgebra

const SUITE = BenchmarkGroup()

include("SW.jl")
include("numeric_conversion.jl")

# benchmark_numeric_conversion!(SUITE) # need a more complex example (this one is too fast)
benchmark_Schrieffer_Wolf!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
