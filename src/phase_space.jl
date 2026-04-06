"""
    PhaseSpace <: HilbertSpace

Hilbert space for position and momentum (quadrature) operators.
"""
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

"""
    Position <: QSym

Position (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Position <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Position(name::Symbol, si::Int, ci::Int) = Position(name, si, ci, NO_INDEX)
Position(name::Symbol, si::Int) = Position(name, si, 1, NO_INDEX)

"""
    Momentum <: QSym

Momentum (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Momentum <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
    index::Index
end
Momentum(name::Symbol, si::Int, ci::Int) = Momentum(name, si, ci, NO_INDEX)
Momentum(name::Symbol, si::Int) = Momentum(name, si, 1, NO_INDEX)

# Construction from Hilbert spaces
Position(h::PhaseSpace, name::Symbol) = Position(name, 1)
Momentum(h::PhaseSpace, name::Symbol) = Momentum(name, 1)

function Position(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Position(name, idx)
end
function Momentum(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    _unwrap_space(h.spaces[idx]) isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Momentum(name, idx)
end

# IndexedOperator convenience
IndexedOperator(op::Position, i::Index) = Position(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Momentum, i::Index) = Momentum(op.name, op.space_index, op.copy_index, i)
IndexedOperator(op::Position, k::Int) = Position(op.name, op.space_index, k, NO_INDEX)
IndexedOperator(op::Momentum, k::Int) = Momentum(op.name, op.space_index, k, NO_INDEX)

# Adjoint — Hermitian
Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

# Equality
Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)

# Hashing
Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))

# Ladder
ladder(::Position) = 0
ladder(::Momentum) = 1
