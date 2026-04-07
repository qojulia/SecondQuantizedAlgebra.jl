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
    copy_index::Int
    index::Index
end
Destroy(name::Symbol, si::Int, ci::Int) = Destroy(name, si, ci, NO_INDEX)
Destroy(name::Symbol, si::Int) = Destroy(name, si, 1, NO_INDEX)

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
    copy_index::Int
    index::Index
end
Create(name::Symbol, si::Int, ci::Int) = Create(name, si, ci, NO_INDEX)
Create(name::Symbol, si::Int) = Create(name, si, 1, NO_INDEX)

# Construction from Hilbert spaces
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Create(name, idx)
end

"""
    IndexedOperator(op::QSym, i::Index) -> QSym
    IndexedOperator(op::QSym, k::Int) -> QSym

Attach a symbolic or concrete index to an operator.

- `IndexedOperator(op, i::Index)`: set the symbolic summation index to `i`,
  preserving the operator's `copy_index`.
- `IndexedOperator(op, k::Int)`: set `copy_index = k` and clear the symbolic
  index (sets `index = NO_INDEX`).

Returns a new operator of the same type with the updated index fields.

# Examples
```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
i = Index(h, :i, 5, h)
a_i = IndexedOperator(a, i)   # symbolic: a_i
a_2 = IndexedOperator(a, 2)   # concrete: copy 2
```
"""
function IndexedOperator end
IndexedOperator(op::Destroy, i::Index) = Destroy(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Create, i::Index) = Create(op.name, op.space_index, op.copy_index, i)

# NumberedOperator: concrete integer index sets copy_index, clears symbolic index
IndexedOperator(op::Destroy, k::Int) = Destroy(op.name, op.space_index, k, NO_INDEX)
IndexedOperator(op::Create, k::Int) = Create(op.name, op.space_index, k, NO_INDEX)

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.copy_index, op.index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.copy_index, op.index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))

# Canonical ordering
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within operator product sequences.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1
