"""
    Destroy <: QSym

Bosonic annihilation operator.
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Destroy(name::Symbol, space_index::Int) = Destroy(name, space_index, 1)

"""
    Create <: QSym

Bosonic creation operator.
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
Create(name::Symbol, space_index::Int) = Create(name, space_index, 1)

# Construction from Hilbert spaces (validation, then discard)
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

# Adjoint
Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.copy_index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.copy_index)

# Equality
Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

# Hashing
Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.copy_index, h))))

# Canonical ordering: Create (0) before Destroy (1), then by space_index, then name
"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within `QMul.args_nc`.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1
