"""
    NLevelSpace(name::Symbol, n::Int)
    NLevelSpace(name::Symbol, n::Int, ground_state::Int)
    NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space for an N-level system (atoms, qubits, qudits, etc.).

Supports [`Transition`](@ref) operators ``|i\\rangle\\langle j|`` with the composition
rule ``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``.

The `ground_state` (default `1`) is used by [`normal_order`](@ref) and [`simplify`](@ref)
with a Hilbert space argument to eliminate ground-state projectors via the completeness
relation ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g} |k\\rangle\\langle k|``.

Levels can be integers (default `1:n`) or symbolic names:

# Examples
```julia
NLevelSpace(:atom, 3)              # levels 1, 2, 3 with ground state 1
NLevelSpace(:atom, 3, 2)           # ground state = level 2
NLevelSpace(:atom, (:g, :e, :a))   # symbolic level names
```
"""
struct NLevelSpace <: HilbertSpace
    name::Symbol
    n::Int
    ground_state::Int
    levels::Vector{Symbol}
    function NLevelSpace(name::Symbol, n::Int, ground_state::Int, levels::Vector{Symbol})
        1 <= ground_state <= n || throw(ArgumentError("Ground state $ground_state out of range 1:$n"))
        isempty(levels) || length(levels) == n || throw(ArgumentError("levels length $(length(levels)) != n=$n"))
        return new(name, n, ground_state, levels)
    end
end
const _EMPTY_LEVELS = Symbol[]
NLevelSpace(name::Symbol, n::Int, ground_state::Int) = NLevelSpace(name, n, ground_state, _EMPTY_LEVELS)
NLevelSpace(name::Symbol, n::Int) = NLevelSpace(name, n, 1, _EMPTY_LEVELS)
function NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})
    return NLevelSpace(name, length(levels), 1, collect(levels))
end

"""
    _level_index(h::NLevelSpace, s::Symbol) -> Int

Look up the integer index for a symbolic level name.
"""
function _level_index(h::NLevelSpace, s::Symbol)
    isempty(h.levels) && throw(ArgumentError("NLevelSpace $(h.name) has no symbolic levels"))
    idx = findfirst(==(s), h.levels)
    idx === nothing && throw(ArgumentError("Level :$s not found in $(h.levels)"))
    return idx
end

Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.name == b.name && a.n == b.n && a.ground_state == b.ground_state && a.levels == b.levels
Base.hash(a::NLevelSpace, h::UInt) = hash(:NLevelSpace, hash(a.name, hash(a.n, hash(a.ground_state, hash(a.levels, h)))))

"""
    Transition <: QSym

Transition operator ``|i\\rangle\\langle j|`` on an [`NLevelSpace`](@ref).

Satisfies the composition rule ``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``.
The adjoint is ``|i\\rangle\\langle j|^\\dagger = |j\\rangle\\langle i|``.

# Construction
```julia
h = NLevelSpace(:atom, 3)
σ12 = Transition(h, :σ, 1, 2)          # |1⟩⟨2| with integer levels

h_sym = NLevelSpace(:atom, (:g, :e))
σge = Transition(h_sym, :σ, :g, :e)    # |g⟩⟨e| with symbolic levels

hp = FockSpace(:c) ⊗ NLevelSpace(:a, 2)
σ = Transition(hp, :σ, 1, 2, 2)        # on 2nd subspace of ProductSpace
```
Or via the [`@qnumbers`](@ref) macro:
```julia
@qnumbers σ::Transition(h, :σ, 1, 2)
```
"""
struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    copy_index::Int
    index::Index
end
Transition(name::Symbol, i::Int, j::Int, si::Int, ci::Int) = Transition(name, i, j, si, ci, NO_INDEX)
Transition(name::Symbol, i::Int, j::Int, si::Int) = Transition(name, i, j, si, 1, NO_INDEX)

# Construction from Hilbert spaces
function Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return Transition(name, i, j, 1)
end
function Transition(h::NLevelSpace, name::Symbol, i::Symbol, j::Symbol)
    return Transition(name, _level_index(h, i), _level_index(h, j), 1)
end
function Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = _unwrap_space(h.spaces[idx])
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return Transition(name, i, j, idx)
end
function Transition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = _unwrap_space(h.spaces[idx])
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    return Transition(name, _level_index(space, i), _level_index(space, j), idx)
end

# IndexedOperator convenience
IndexedOperator(op::Transition, i::Index) = Transition(op.name, op.i, op.j, op.space_index, op.copy_index, i)
IndexedOperator(op::Transition, k::Int) = Transition(op.name, op.i, op.j, op.space_index, k, NO_INDEX)

# Adjoint: |i⟩⟨j|† = |j⟩⟨i|
Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index, op.copy_index, op.index)

# Equality
Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.copy_index == b.copy_index && a.index == b.index
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

# Hashing
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.copy_index, hash(a.index, h)))))))

# Ladder (not applicable to Transition)
ladder(::Transition) = 0
