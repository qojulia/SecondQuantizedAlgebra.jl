"""
    NLevelSpace <: HilbertSpace

Hilbert space for N-level systems (atoms, qubits, etc.).
"""
struct NLevelSpace <: HilbertSpace
    name::Symbol
    n::Int
    ground_state::Int
    function NLevelSpace(name::Symbol, n::Int, ground_state::Int)
        1 <= ground_state <= n || throw(ArgumentError("Ground state $ground_state out of range 1:$n"))
        return new(name, n, ground_state)
    end
end
Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.name == b.name && a.n == b.n && a.ground_state == b.ground_state
Base.hash(a::NLevelSpace, h::UInt) = hash(:NLevelSpace, hash(a.name, hash(a.n, hash(a.ground_state, h))))

"""
    Transition <: QSym

Transition operator |i⟩⟨j| on an [`NLevelSpace`](@ref).
"""
struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
end

# Construction from Hilbert spaces
function Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return Transition(name, i, j, 1)
end
function Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return Transition(name, i, j, idx)
end

# Adjoint: |i⟩⟨j|† = |j⟩⟨i|
Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index)

# Equality
Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

# Hashing
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, h)))))

# Ladder (not applicable to Transition)
ladder(::Transition) = 0
