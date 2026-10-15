"""
    PhaseSpace(name::Symbol) <: HilbertSpace

Hilbert space of one mode in quadrature form, home of [`Position`](@ref) and
[`Momentum`](@ref).

A `PhaseSpace` holds a single canonical pair. The quadrature names are only
labels: a second `Position` under another name is the same observable. Several
modes therefore need several spaces combined with [`⊗`](@ref).

# Examples

```jldoctest
julia> h = PhaseSpace(:osc)
ℋ(osc)

julia> x = Position(h, :x); y = Position(h, :y); p = Momentum(h, :p);

julia> commutator(x, p), commutator(y, p)
(im, im)
```

See also [`Position`](@ref), [`Momentum`](@ref), [`FockSpace`](@ref).
"""
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
PhaseSpace(name::AbstractString) = name_must_be_symbol(name)
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

"""
    Position(h::PhaseSpace, name::Symbol)
    Position(h::ProductSpace, name::Symbol, [subsystem::Int])

Position quadrature of a mode.

The operator ``x = (a + a^\\dagger)/\\sqrt{2}`` is Hermitian and canonically
conjugate to [`Momentum`](@ref), ``[x, p] = i``. Products are ordered with position
to the left of momentum. On a [`ProductSpace`](@ref), `subsystem` selects the
[`PhaseSpace`](@ref) factor and may be omitted when there is only one.

# Examples

```jldoctest
julia> h = PhaseSpace(:osc);

julia> @qnumbers x::Position(h) p::Momentum(h);

julia> commutator(x, p)
im
```

See also [`Momentum`](@ref), [`PhaseSpace`](@ref).
"""
Position(name::Symbol, si::Int, idx::Index) = Op(OP_POSITION, name_id(name), si, idx, 0, 0, 0, 0)
Position(name::Symbol, si::Int) = Position(name, si, NO_INDEX)

"""
    Momentum(h::PhaseSpace, name::Symbol)
    Momentum(h::ProductSpace, name::Symbol, [subsystem::Int])

Momentum quadrature of a mode.

The operator ``p = i(a^\\dagger - a)/\\sqrt{2}`` is Hermitian and canonically
conjugate to [`Position`](@ref). Arguments are as for [`Position`](@ref).

# Examples

```jldoctest
julia> h = PhaseSpace(:osc);

julia> @qnumbers x::Position(h) p::Momentum(h);

julia> p' == p
true

julia> p * x
-im + x * p
```

See also [`Position`](@ref), [`PhaseSpace`](@ref).
"""
Momentum(name::Symbol, si::Int, idx::Index) = Op(OP_MOMENTUM, name_id(name), si, idx, 0, 0, 0, 0)
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

Position(h::ProductSpace, name::Symbol) = Position(h, name, unique_subspace_index(h, PhaseSpace))
Momentum(h::ProductSpace, name::Symbol) = Momentum(h, name, unique_subspace_index(h, PhaseSpace))

Position(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
Momentum(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
