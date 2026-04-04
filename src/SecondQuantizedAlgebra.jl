module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils, BasicSymbolic, arguments, iscall, operation, substitute
using Symbolics: Symbolics
using TermInterface: TermInterface

using Combinatorics: combinations

using SciMLBase: SciMLBase
using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using LaTeXStrings: LaTeXStrings, @L_str, latexstring
using Latexify: Latexify, latexify, @latexrecipe
using MacroTools: MacroTools

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

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor,
    Destroy, Create,
    QMul, QAdd,
    @qnumbers,
    normal_order, simplify

end
