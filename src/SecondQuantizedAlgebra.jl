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

export QField, QSym, QTerm,
    HilbertSpace, FockSpace, ProductSpace,
    ⊗, tensor

end
