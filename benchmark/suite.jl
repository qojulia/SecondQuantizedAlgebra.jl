using BenchmarkTools: BenchmarkGroup, @benchmarkable
using SecondQuantizedAlgebra: SecondQuantizedAlgebra, FockSpace, NLevelSpace,
    ProductSpace, Destroy, Transition, SpinSpace, Spin, PauliSpace, Pauli,
    PhaseSpace, Position, Momentum, Index, IndexedVariable, IndexedOperator,
    DoubleIndexedVariable, Σ, ⊗, commutator, simplify, normal_order,
    expand_completeness, Displace, Rotation, Squeeze, transform, conjugate,
    generators, gauge_term, substitute, expim, taylor, dissipator,
    hamiltonian_action, to_numeric, QuantumOpticsBackend
using Symbolics: @variables, derivative

include("commutator.jl")
include("simplify_and_normal_order.jl")
include("indexing.jl")
include("accumulation.jl")
include("unitary.jl")
include("operator_functions.jl")
include("superoperators.jl")

function validate_suite(group::BenchmarkGroup, path::Vector{String} = String[])
    isempty(group) && error("Empty benchmark group: $(join(path, '/'))")
    for (key, value) in group
        key isa String && !isempty(key) && !occursin('/', key) ||
            error("Benchmark keys must be nonempty strings without '/': $key")
        child = [path; key]
        if value isa BenchmarkGroup
            validate_suite(value, child)
        else
            length(child) == 3 || error("Expected topic/case/legend: $(join(child, '/'))")
        end
    end
    return group
end

function benchmark_suite()
    suite = BenchmarkGroup()
    benchmark_commutator!(suite)
    benchmark_simplify_and_normal_order!(suite)
    benchmark_indexing!(suite)
    benchmark_accumulation!(suite)
    benchmark_unitary!(suite)
    benchmark_operator_functions!(suite)
    benchmark_superoperators!(suite)
    return validate_suite(suite)
end
