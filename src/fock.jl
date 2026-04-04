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
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
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

