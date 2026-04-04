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
include("qmul.jl")
include("qadd.jl")
include("macros.jl")
include("interface.jl")
include("normal_order.jl")
include("simplify.jl")
include("printing.jl")
include("latexify_recipes.jl")
include("numeric.jl")

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    normal_order, simplify,
    to_numeric, numeric_average

end
