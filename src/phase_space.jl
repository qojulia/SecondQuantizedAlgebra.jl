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
end

"""
    Momentum <: QSym

Momentum (quadrature) operator on a [`PhaseSpace`](@ref).
"""
struct Momentum <: QSym
    name::Symbol
    space_index::Int
end

# Construction from Hilbert spaces
Position(h::PhaseSpace, name::Symbol) = Position(name, 1)
Momentum(h::PhaseSpace, name::Symbol) = Momentum(name, 1)

function Position(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Position(name, idx)
end
function Momentum(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa PhaseSpace || throw(ArgumentError("Space at index $idx is not a PhaseSpace"))
    return Momentum(name, idx)
end

# Adjoint — Hermitian
Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

# Equality
Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)

# Hashing
Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, h)))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, h)))

# Ladder (not applicable — Hermitian operators, no creation/annihilation distinction)
ladder(::Position) = 0
ladder(::Momentum) = 1
