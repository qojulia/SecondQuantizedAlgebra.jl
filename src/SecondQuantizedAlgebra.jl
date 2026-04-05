module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, simplify
using Symbolics: Symbolics, Num, expand
using TermInterface: TermInterface

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using Latexify: Latexify, latexify, @latexrecipe

include("types.jl")
include("hilbertspace.jl")
include("index.jl")
include("fock.jl")
include("nlevel.jl")
include("pauli.jl")
include("spin.jl")
include("phase_space.jl")
include("cluster.jl")
include("indexed_variables.jl")
include("qmul.jl")
include("qadd.jl")
include("qsum.jl")
include("macros.jl")
include("interface.jl")
include("simplify.jl")
include("commutator.jl")
include("change_index.jl")
include("expand_sums.jl")
include("normal_order.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")

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
    normal_order, simplify, expand, commutator,
    to_numeric, numeric_average,
    transition_superscript

end
