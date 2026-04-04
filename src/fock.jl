"""
    Destroy <: QSym

Bosonic annihilation operator. Stores only name and space index.
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
end

"""
    Create <: QSym

Bosonic creation operator. Stores only name and space index.
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
end

# Construction from Hilbert spaces (validation, then discard)
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"
    @assert h.spaces[idx] isa FockSpace "Space at index $idx is not a FockSpace"
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    @assert 1 <= idx <= length(h.spaces) "Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"
    @assert h.spaces[idx] isa FockSpace "Space at index $idx is not a FockSpace"
    return Create(name, idx)
end

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, h)))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, h)))

# Canonical ordering: Create (0) before Destroy (1), then by space_index, then name
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within `QMul.args_nc`.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1

"""
    canonical_lt(a::QSym, b::QSym)

Canonical ordering comparator: sort by space_index only.
Operators on the same space are NOT reordered — their relative order
is preserved so that `normal_order()` can detect and apply commutation relations.
"""
function canonical_lt(a::QSym, b::QSym)
    return a.space_index < b.space_index
end
