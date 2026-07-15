module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, simplify, substitute, add_worker
using Symbolics: Symbolics, Num, expand, @variables, build_function, symbolic_to_float
using TermInterface: TermInterface

# `CNum` (the coefficient type `Coeff`) is defined in `expressions/cnum.jl`.

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor, expect

using Combinatorics: with_replacement_combinations
using FunctionWrappers: FunctionWrapper
using Latexify: Latexify, latexify, @latexrecipe
using PrecompileTools: @setup_workload, @compile_workload
using SciMLPublic: @public
import MutableArithmetics as MA

include("types.jl")
include("operators/hilbertspace.jl")
include("expressions/intern.jl")
include("expressions/index_types.jl")
include("operators/op.jl")

include("operators/fock.jl")
include("operators/nlevel.jl")
include("operators/pauli.jl")
include("operators/spin.jl")
include("operators/phase_space.jl")
include("operators/operators.jl")
include("expressions/monomial.jl")
include("expressions/cnum.jl")
include("expressions/qterm.jl")
include("expressions/qadd.jl")

include("algebra/passes.jl")
include("algebra/pipelines.jl")

include("expressions/index.jl")

include("algebra/algebra.jl")
include("algebra/mutable_arithmetics.jl")
include("algebra/weyl.jl")

include("average.jl")
include("numeric.jl")

include("printing/printing.jl")
include("printing/latexify_recipes.jl")

"""
    @qnumbers ops...

Convenience macro for constructing named quantum operators.

Each argument has the form `name::OperatorType(hilbert_space, args...)`. The macro
calls `OperatorType(hilbert_space, :name, args...)` and binds the result to `name`
in the calling scope. Multiple operators can be declared in one call.

# Examples
```jldoctest
julia> h = FockSpace(:fock);

julia> @qnumbers a::Destroy(h)
(a,)
```

See also [`Destroy`](@ref), [`Transition`](@ref), [`CollectiveTransition`](@ref),
[`Pauli`](@ref), [`Spin`](@ref).
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
    NLevelSpace, Transition, CollectiveNLevelSpace, CollectiveTransition,
    PauliSpace, Pauli,
    SpinSpace, Spin,
    PhaseSpace, Position, Momentum,
    Index, has_index, index_slot, index_range, index_name, index_sym, IndexedOperator,
    IndexedVariable, DoubleIndexedVariable,
    Σ, ∑, change_index, get_indices,
    ⊗, tensor, Destroy, Create,
    @qnumbers, @variables,
    average, undo_average, make_time_dependent,
    acts_on, is_average, is_indexed_sum,
    fundamental_operators, find_operators, unique_ops,
    prefactor, operators,
    substitute,
    normal_order, normal_to_symmetric, symmetric_to_normal,
    simplify, expand, expand_completeness, assume_distinct_index, commutator, anticommutator,
    to_numeric, numeric_average,
    qadjoint, qconj, dagger, inner_adjoint,
    Op, operator_name, is_destroy, is_create, is_transition, is_collective_transition,
    is_pauli, is_spin, is_position, is_momentum, optype


# Public API that is intentionally NOT exported — accessed as
# `SecondQuantizedAlgebra.symbol`.
@public HilbertSpace, QField, QSym, OpKind,
    OP_DESTROY, OP_CREATE, OP_TRANSITION, OP_PAULI, OP_SPIN, OP_POSITION, OP_MOMENTUM,
    OP_COLLECTIVE_TRANSITION,
    QAdd, QTerm, QTermDict, Coeff, CNum,
    has_sum_metadata, get_sum_indices, get_sum_non_equal,
    transition_superscript, constraint_pairs,
    to_num, order_key, term_order_key, qadd_order_key

include("precompile.jl")

end
