module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, simplify
using Symbolics: Symbolics, Num, expand
using TermInterface: TermInterface

const CNum = Complex{Num}

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using Combinatorics: with_replacement_combinations
using Latexify: Latexify, latexify, @latexrecipe

include("types.jl")
include("hilbertspace.jl")
include("index_types.jl")
include("fock.jl")
include("nlevel.jl")
include("pauli.jl")
include("spin.jl")
include("phase_space.jl")
include("cluster.jl")
include("qmul.jl")
include("qadd.jl")
include("index.jl")
include("interface.jl")
include("simplify.jl")
include("commutator.jl")
include("normal_order.jl")
include("average.jl")
include("operators.jl")
include("numeric.jl")
include("printing.jl")
include("latexify_recipes.jl")

"""
    @qnumbers

Convenience macro for the construction of operators.

Examples
========
```julia
h = FockSpace(:fock)
@qnumbers a::Destroy(h)

h = FockSpace(:one) ⊗ FockSpace(:two)
@qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
```
"""
macro qnumbers(qs...)
    ex = Expr(:block)
    qnames = []
    for q in qs
        @assert q isa Expr && q.head == :(::)
        name = q.args[1]
        @assert name isa Symbol
        push!(qnames, name)
        f = q.args[2]
        @assert f isa Expr && f.head == :call
        op_type = f.args[1]
        op_args = f.args[2:end]
        name_quoted = Expr(:quote, name)
        # Insert name as second argument: Op(hilbert, name, extra_args...)
        construction = Expr(:call, esc(op_type), esc(op_args[1]), name_quoted, map(esc, op_args[2:end])...)
        push!(ex.args, :($(esc(name)) = $(construction)))
    end
    push!(ex.args, Expr(:tuple, map(esc, qnames)...))
    return ex
end


export FockSpace, ProductSpace,
    NLevelSpace, Transition,
    PauliSpace, Pauli,
    SpinSpace, Spin,
    PhaseSpace, Position, Momentum,
    ClusterSpace, cluster_expand, has_cluster,
    Index, has_index, IndexedOperator,
    IndexedVariable, DoubleIndexedVariable,
    Σ, ∑, change_index, expand_sums, get_indices,
    ⊗, tensor,
    Destroy, Create,
    @qnumbers,
    NormalOrder,
    average, undo_average,
    acts_on, is_average,
    has_sum_metadata, get_sum_indices, get_sum_non_equal,
    fundamental_operators, find_operators, unique_ops,
    normal_order, simplify, expand, commutator,
    to_numeric, numeric_average,
    transition_superscript

end
