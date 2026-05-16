module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, simplify, substitute
using Symbolics: Symbolics, Num, expand, @variables
using TermInterface: TermInterface

const CNum = Complex{Num}

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor, expect

using Combinatorics: with_replacement_combinations
using Latexify: Latexify, latexify, @latexrecipe
using PrecompileTools: @setup_workload, @compile_workload

include("types.jl")
include("operators/hilbertspace.jl")
include("expressions/index_types.jl")

include("operators/fock.jl")
include("operators/nlevel.jl")
include("operators/pauli.jl")
include("operators/spin.jl")
include("operators/phase_space.jl")
include("operators/operators.jl")
include("expressions/cnum.jl")
include("expressions/qterm.jl")
include("expressions/qadd.jl")

include("algebra/passes.jl")
include("algebra/pipelines.jl")

include("expressions/index.jl")

include("algebra/algebra.jl")
include("algebra/weyl.jl")

include("average.jl")
include("numeric.jl")

include("printing/printing.jl")
include("printing/latexify_recipes.jl")

"""
    @qnumbers ops...

Convenience macro for constructing named quantum operators.

Each argument has the form `name::OperatorType(hilbert_space, args...)`. The macro
creates the operator via `OperatorType(hilbert_space, :name, args...)` and binds it
to `name` in the calling scope.

# Examples
```julia
h = FockSpace(:fock)
@qnumbers a::Destroy(h)               # creates Destroy(h, :a) bound to `a`

h = FockSpace(:one) ⊗ FockSpace(:two)
@qnumbers a::Destroy(h, 1) b::Destroy(h, 2)   # operators on subspaces 1 and 2

h = NLevelSpace(:atom, 2)
@qnumbers σ::Transition(h, 1, 2)  # creates Transition(h, :σ, 1, 2) bound to `σ`
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
    Index, has_index, IndexedOperator,
    IndexedVariable, DoubleIndexedVariable,
    Σ, ∑, change_index, get_indices,
    ⊗, tensor, Destroy, Create,
    @qnumbers, @variables,
    average, undo_average,
    acts_on, is_average,
    fundamental_operators, find_operators, unique_ops,
    prefactor, operators,
    substitute,
    normal_order, normal_to_symmetric, symmetric_to_normal,
    simplify, expand, expand_completeness, assume_distinct_index, commutator, anticommutator,
    to_numeric, numeric_average


# Public API that is intentionally NOT exported — accessed as
# `SecondQuantizedAlgebra.symbol`. The `public` keyword is Julia ≥ 1.11; on
# 1.10 the contract is documented in `docs/src/API.md`.
macro public(ex)
    return if VERSION >= v"1.11.0-DEV.469"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("something informative")
        esc(Expr(:public, args...))
    else
        nothing
    end
end

@public HilbertSpace, QField, QSym,
    QAdd, QTerm, QTermDict, has_sum_metadata, get_sum_indices, get_sum_non_equal,
    transition_superscript, constraint_pairs

include("precompile.jl")

end
