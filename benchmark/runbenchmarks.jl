using BenchmarkTools
using SecondQuantizedAlgebra
using Latexify: latexify
using Symbolics: @variables

const SUITE = BenchmarkGroup()

include("commutator.jl")
include("simplify_and_normal_order.jl")
include("indexing.jl")
include("accumulation.jl")
include("unitary.jl")

benchmark_commutator!(SUITE)
benchmark_simplify_and_normal_order!(SUITE)
benchmark_indexing!(SUITE)
benchmark_accumulation!(SUITE)
benchmark_unitary!(SUITE)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(minimum(results))

BenchmarkTools.save("benchmarks_output.json", minimum(results))
