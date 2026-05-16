"""
    PhaseSpace(name::Symbol)

Hilbert space for position and momentum (quadrature) operators.

Supports [`Position`](@ref) and [`Momentum`](@ref) operators with the canonical
commutation relation ``[p, x] = -i`` (eagerly applied by `*`: position is
ordered left of momentum).

# Examples
```julia
h = PhaseSpace(:osc)
@qnumbers x::Position(h) p::Momentum(h)
```
"""
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

"""
    Position <: QSym

Position (quadrature) operator ``x`` on a [`PhaseSpace`](@ref).

Hermitian (`x' == x`). Related to Fock operators by ``x = (a + a^\\dagger)/\\sqrt{2}``.

# Construction
```julia
h = PhaseSpace(:osc)
x = Position(h, :x)             # single space
hp = PhaseSpace(:a) ⊗ PhaseSpace(:b)
x1 = Position(hp, :x, 1)        # first subspace
```
"""
struct Position <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Position(name::Symbol, si::Int) = Position(name, si, NO_INDEX)

"""
    Momentum <: QSym

Momentum (quadrature) operator ``p`` on a [`PhaseSpace`](@ref).

Hermitian (`p' == p`). Related to Fock operators by ``p = i(a^\\dagger - a)/\\sqrt{2}``.

# Construction
```julia
h = PhaseSpace(:osc)
p = Momentum(h, :p)             # single space
hp = PhaseSpace(:a) ⊗ PhaseSpace(:b)
p2 = Momentum(hp, :p, 2)        # second subspace
```
"""
struct Momentum <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Momentum(name::Symbol, si::Int) = Momentum(name, si, NO_INDEX)

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

# Auto-detect subspace when the ProductSpace contains exactly one PhaseSpace.
Position(h::ProductSpace, name::Symbol) = Position(h, name, _unique_subspace_index(h, PhaseSpace))
Momentum(h::ProductSpace, name::Symbol) = Momentum(h, name, _unique_subspace_index(h, PhaseSpace))

IndexedOperator(op::Position, i::Index) = Position(op.name, op.space_index, i)
IndexedOperator(op::Momentum, i::Index) = Momentum(op.name, op.space_index, i)

Base.adjoint(op::Position) = op
Base.adjoint(op::Momentum) = op

Base.isequal(a::Position, b::Position) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.isequal(a::Momentum, b::Momentum) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::Position, b::Position) = isequal(a, b)
Base.:(==)(a::Momentum, b::Momentum) = isequal(a, b)

Base.hash(a::Position, h::UInt) = hash(:Position, hash(a.name, hash(a.space_index, hash(a.index, h))))
Base.hash(a::Momentum, h::UInt) = hash(:Momentum, hash(a.name, hash(a.space_index, hash(a.index, h))))

ladder(::Position) = 0
ladder(::Momentum) = 1

# --- Operator hooks ---

function _site_compare(a::Position, b::Position, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    return a.index == b.index ? Equal : Undetermined
end
function _site_compare(a::Momentum, b::Momentum, ne::Vector{NonEqualPair})::SiteCmp
    return _site_compare(
        Position(a.name, a.space_index, a.index),
        Position(b.name, b.space_index, b.index), ne
    )
end

# Cross-type same-site test ignores name: position x and momentum p on the same
# Hilbert subspace and index are conjugate variables on one site.
function _site_compare(a::Position, b::Momentum, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    return a.index == b.index ? Equal : Undetermined
end
function _site_compare(a::Momentum, b::Position, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    return a.index == b.index ? Equal : Undetermined
end

# Same-site P·X = X·P - i (commute residual is -i times identity).
_can_commute(a::Position, b::Momentum) = true
_can_commute(a::Momentum, b::Position) = false
_can_commute(a::Position, b::Position) = true
_can_commute(a::Momentum, b::Momentum) = true

_commute_pair(a::Momentum, b::Position) = (b, a, _to_cnum(-im), _EMPTY_OPS)   # P·X = X·P - i·I

_reduce_pair(::Position, ::Momentum) = nothing
_reduce_pair(::Momentum, ::Position) = nothing
_reduce_pair(::Position, ::Position) = nothing
_reduce_pair(::Momentum, ::Momentum) = nothing
