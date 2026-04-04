module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics
using TermInterface: TermInterface

using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using Latexify: Latexify, latexify, @latexrecipe

include("types.jl")
include("hilbertspace.jl")
include("fock.jl")
include("nlevel.jl")
include("pauli.jl")
include("spin.jl")
include("phase_space.jl")
include("qmul.jl")
include("qadd.jl")
include("macros.jl")
include("interface.jl")
include("simplify.jl")
include("normal_order.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    NLevelSpace, Transition,
    PauliSpace, Pauli,
    SpinSpace, Spin,
    PhaseSpace, Position, Momentum,
    ⊗, tensor,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    OrderingConvention, NormalOrder,
    normal_order, simplify,
    to_numeric, numeric_average,
    transition_superscript

end
