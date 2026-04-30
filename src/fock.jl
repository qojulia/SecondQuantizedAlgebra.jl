"""
    Destroy <: QSym

Bosonic annihilation operator ``a`` on a [`FockSpace`](@ref).

Satisfies the canonical commutation relation ``[a, a^\\dagger] = 1``.
The adjoint `a'` returns the corresponding [`Create`](@ref) operator.

# Construction
```julia
h = FockSpace(:cavity)
a = Destroy(h, :a)              # single-space
h2 = FockSpace(:a) ⊗ FockSpace(:b)
a = Destroy(h2, :a, 1)          # on first subspace of ProductSpace
```
Or via the [`@qnumbers`](@ref) macro:
```julia
@qnumbers a::Destroy(h)
```
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Destroy(name::Symbol, si::Int) = Destroy(name, si, NO_INDEX)

"""
    Create <: QSym

Bosonic creation operator ``a^\\dagger`` on a [`FockSpace`](@ref).

The adjoint of [`Destroy`](@ref). Satisfies `[a, a'] = 1` under [`NormalOrder`](@ref).
Constructed implicitly via `adjoint(::Destroy)`, or directly with the same signatures
as [`Destroy`](@ref).
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Create(name::Symbol, si::Int) = Create(name, si, NO_INDEX)

# Construction from Hilbert spaces
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Create(name, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one FockSpace.
Destroy(h::ProductSpace, name::Symbol) = Destroy(h, name, _unique_subspace_index(h, FockSpace))
Create(h::ProductSpace, name::Symbol) = Create(h, name, _unique_subspace_index(h, FockSpace))

"""
    IndexedOperator(op::QSym, i::Index) -> QSym

Attach a symbolic summation index to an operator.

Returns a new operator of the same type with `index = i`. To clear the index, use
`IndexedOperator(op, NO_INDEX)`.

# Examples
```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
i = Index(h, :i, 5, h)
a_i = IndexedOperator(a, i)   # symbolic: a_i
```
"""
function IndexedOperator end
IndexedOperator(op::Destroy, i::Index) = Destroy(op.name, op.space_index, i)
IndexedOperator(op::Create, i::Index) = Create(op.name, op.space_index, i)

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.index, h))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.index, h))))

# Canonical ordering
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within operator product sequences.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1
