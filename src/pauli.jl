"""
    PauliSpace <: HilbertSpace

Hilbert space for two-level Pauli operators.
"""
struct PauliSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

"""
    Pauli <: QSym

Pauli operator (σx, σy, σz) on a [`PauliSpace`](@ref).
Axis: 1=x, 2=y, 3=z.
"""
struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
    index::Index
    function Pauli(name::Symbol, axis::Int, si::Int, ci::Int, idx::Index)
        1 <= axis <= 3 || throw(ArgumentError("Pauli axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, si, ci, idx)
    end
end
Pauli(name::Symbol, axis::Int, si::Int, ci::Int) = Pauli(name, axis, si, ci, NO_INDEX)
Pauli(name::Symbol, axis::Int, si::Int) = Pauli(name, axis, si, 1, NO_INDEX)

# Construction from Hilbert spaces
Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    _unwrap_space(h.spaces[idx]) isa PauliSpace || throw(ArgumentError("Space at index $idx is not a PauliSpace"))
    return Pauli(name, axis, idx)
end

# IndexedOperator convenience
IndexedOperator(op::Pauli, i::Index) = Pauli(op.name, op.axis, op.space_index, op.copy_index, i)

# Adjoint — Hermitian
Base.adjoint(op::Pauli) = op

# Equality
Base.isequal(a::Pauli, b::Pauli) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Pauli, b::Pauli) = isequal(a, b)

# Hashing
Base.hash(a::Pauli, h::UInt) = hash(:Pauli, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.copy_index, hash(a.index, h))))))

# Ladder (not applicable)
ladder(::Pauli) = 0
