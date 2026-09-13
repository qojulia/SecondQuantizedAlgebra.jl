using BenchmarkTools
using SecondQuantizedAlgebra
using Symbolics: @variables, derivative

const SUITE = BenchmarkGroup()
include("operator_functions.jl")
benchmark_operator_functions!(SUITE)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 2000
BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
minimums = minimum(results)
display(minimums)

output = isempty(ARGS) ? "qexpr_research_output.json" : only(ARGS)
BenchmarkTools.save(output, minimums)
