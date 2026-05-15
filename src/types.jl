"""
    QField

Abstract supertype for all second-quantized operator expressions.

Subtypes:
- [`QSym`](@ref): leaf operators (e.g. [`Destroy`](@ref), [`Create`](@ref), [`Transition`](@ref))
- [`QAdd`](@ref): compound expressions (a sum of [`QTerm`](@ref) products)

Supports arithmetic (`+`, `-`, `*`, `^`, `/`), `adjoint`, and comparison via `==`/`isequal`.
All arithmetic eagerly applies normal ordering and returns [`QAdd`](@ref).
"""
abstract type QField end

"""
    QSym <: QField

Abstract type for fundamental (leaf) operators in the expression tree.

Every `QSym` carries three fields identifying its site:
- `name::Symbol` — display name
- `space_index::Int` — position in [`ProductSpace`](@ref) (1 for single spaces)
- `index::Index` — symbolic summation index ([`Index`](@ref)), or `NO_INDEX`

Concrete subtypes: [`Destroy`](@ref), [`Create`](@ref), [`Transition`](@ref),
[`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref), [`Momentum`](@ref).
"""
abstract type QSym <: QField end

Base.one(::T) where {T <: QField} = one(T)
Base.one(::Type{<:QField}) = 1
Base.isone(::QField) = false
Base.zero(::T) where {T <: QField} = zero(T)
Base.zero(::Type{<:QField}) = 0
