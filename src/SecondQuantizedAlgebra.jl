module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, simplify, substitute, add_worker
using Symbolics: Symbolics, Num, expand, taylor, @variables, build_function, get_variables
using TermInterface: TermInterface

import QuantumInterface
import QuantumInterface: expect, basis, dagger
using QuantumInterface: AbstractOperator, StateVector, Basis, CompositeBasis
import TensorCore: ⊗, tensor

using Combinatorics: combinations, with_replacement_combinations
using Latexify: Latexify, latexify, @latexrecipe
using LaTeXStrings: LaTeXString
using LinearAlgebra: LinearAlgebra, istriu, istril
using PrecompileTools: @setup_workload, @compile_workload
using Primes: isprime
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
include("expressions/gaussian.jl")
include("expressions/exact.jl")
include("expressions/monomial.jl")
include("expressions/coeff.jl")
include("expressions/phase.jl")
include("expressions/recognize.jl")
include("expressions/lowering.jl")
include("expressions/phase_forms.jl")
include("expressions/rational_zero.jl")
include("expressions/reduce.jl")
include("expressions/qterm.jl")
include("expressions/qadd.jl")

include("algebra/passes.jl")
include("algebra/pipelines.jl")

include("expressions/index.jl")

include("algebra/algebra.jl")
include("expressions/coefficient_terms.jl")
include("algebra/mutable_arithmetics.jl")
include("expressions/qexpr.jl")
include("expressions/qexpr_structural.jl")
include("expressions/qexpr_taylor.jl")
include("algebra/weyl.jl")
include("algebra/unitary.jl")
include("algebra/unitary_qexpr.jl")
include("algebra/unitary_affine.jl")
include("algebra/unitary_affine_provenance.jl")
include("algebra/unitary_constructors.jl")
include("algebra/unitary_displacement_frames.jl")
include("algebra/unitary_bogoliubov.jl")
include("algebra/unitary_generators.jl")
include("algebra/harmonics.jl")

include("superoperators/superterm.jl")
include("superoperators/superoperator.jl")
include("superoperators/scope.jl")
include("superoperators/channel.jl")
include("superoperators/lindbladian.jl")
include("superoperators/structural.jl")
include("superoperators/gksl.jl")
include("superoperators/superexpr.jl")
include("superoperators/frames.jl")
include("superoperators/harmonics.jl")

include("average.jl")

include("numeric/backend.jl")
include("numeric/coeff.jl")
include("numeric/core.jl")
include("numeric/indexed.jl")
include("numeric/api.jl")
include("numeric/qexpr.jl")
include("numeric/superoperator.jl")

include("printing/printing.jl")
include("printing/latexify_recipes.jl")
include("printing/qexpr.jl")
include("printing/superoperator.jl")

"""
    @qnumbers ops...

Define operators and bind each to a variable of the same name.

Each argument reads `name::Constructor(space, args...)` and stands for
`name = Constructor(space, :name, args...)`, so the operator name is written once.
Several operators can be declared in one call.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2)
(a, σ₁₂)

julia> a' * σ
a' * σ₁₂
```

See also [`Destroy`](@ref), [`Transition`](@ref), [`Pauli`](@ref), [`Spin`](@ref),
[`Position`](@ref), [`Momentum`](@ref), [`CollectiveTransition`](@ref).
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
    fundamental_operators, find_operators, unique_up_to_adjoint,
    get_prefactor, get_operators, get_variables,
    substitute,
    normal_order, normal_to_symmetric, symmetric_to_normal,
    UnitaryTransform, RotatingFrame, Displace, DisplacementFrame, Rotation, Squeeze, Bogoliubov,
    transform, conjugate, gauge_term, generators, displacement,
    simplify, expand, taylor, expim, expand_completeness, assume_distinct_index, commutator, anticommutator,
    to_numeric, numeric_average,
    NumericBackend, QuantumOpticsBackend, QuantumToolboxBackend,
    numeric_operator, numeric_basis, numeric_subbasis, numeric_embed,
    numeric_identity, numeric_num_subsystems,
    numeric_assemble, numeric_assemble_td, numeric_materialize, EagerOperator, numeric_expect,
    numeric_backend,
    qadjoint, qconj, dagger, inner_adjoint,
    Op, operator_name, operator_index, is_destroy, is_create, is_transition, is_collective_transition,
    is_pauli, is_spin, is_position, is_momentum, optype,
    Lindbladian, hamiltonian_action, dissipator, JumpBasis, kossakowski

@public HilbertSpace, QField, QSym, OpKind,
    NumericContext, expect,
    OP_DESTROY, OP_CREATE, OP_TRANSITION, OP_PAULI, OP_SPIN, OP_POSITION, OP_MOMENTUM,
    OP_COLLECTIVE_TRANSITION,
    QAdd, QTerm, QTermDict, QExpr, Coeff, CNum, PhaseTerm,
    has_sum_metadata, get_sum_indices, get_sum_non_equal, get_sum_body, indexed_sum,
    set_acts_on, rename, n_levels, ground_state, transition_levels, operator_axis,
    CoefficientTerm, coefficient_terms, has_coefficient_terms, coefficient,
    transition_superscript, constraint_pairs,
    to_num, order_key, term_order_key, qadd_order_key,
    exponential_form, trigonometric_form, phase_terms, harmonic, harmonics, secular, map_coefficients,
    SuperOperator, SuperTerm, SuperExpr, Channel, SuperBuilder, LindbladianBuilder,
    spre, spost, sprepost, liouvillian, compose, superterm, left_operator, right_operator,
    hamiltonian, channels, jump, partner, scope, formal_terms, poly, accumulate!, build,
    superterm_order_key, superoperator_order_key, channel_order_key,
    independent_basis, support_basis, GKSLFormError,
    to_numeric_master, numeric_superoperator, numeric_matrix, numeric_scaled_td

include("precompile.jl")

end
