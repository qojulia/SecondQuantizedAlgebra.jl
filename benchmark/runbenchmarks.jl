using BenchmarkTools: BenchmarkTools

include("suite.jl")
const SUITE = benchmark_suite()

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(minimum(results))

BenchmarkTools.save("benchmarks_output.json", minimum(results))
