"""
    PhaseSpace(name::Symbol) <: HilbertSpace

Hilbert space for position and momentum quadratures. Hosts [`Position`](@ref)
and [`Momentum`](@ref) operators satisfying ``[p, x] = -i``; arithmetic
canonicalizes products to place position left of momentum.

# Examples

```jldoctest
julia> PhaseSpace(:osc)
ℋ(osc)
```

See also [`Position`](@ref), [`Momentum`](@ref), [`FockSpace`](@ref).
"""
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PhaseSpace, b::PhaseSpace) = a.name == b.name
Base.hash(a::PhaseSpace, h::UInt) = hash(:PhaseSpace, hash(a.name, h))

"""
    Position(h::PhaseSpace, name::Symbol) -> Op

Position (quadrature) operator ``x`` on a [`PhaseSpace`](@ref). Hermitian
(`x' == x`) and related to Fock operators by
``x = (a + a^\\dagger)/\\sqrt{2}``. Canonical pair with [`Momentum`](@ref):
``[p, x] = -i``. Returns an [`Op`](@ref) tagged `OP_POSITION`.

# Examples

```jldoctest
julia> h = PhaseSpace(:osc);

julia> @qnumbers x::Position(h) p::Momentum(h);

julia> p * x
-im + x * p
```

See also [`Momentum`](@ref), [`PhaseSpace`](@ref).
"""
Position(name::Union{Symbol, Int32}, si::Int, idx::Index) = Op(OP_POSITION, _name_id(name), si, idx, 0, 0, 0, 0)
Position(name::Union{Symbol, Int32}, si::Int) = Position(name, si, NO_INDEX)

"""
    Momentum(h::PhaseSpace, name::Symbol) -> Op

Momentum (quadrature) operator ``p`` on a [`PhaseSpace`](@ref). Hermitian
(`p' == p`) and related to Fock operators by
``p = i(a^\\dagger - a)/\\sqrt{2}``. Canonical pair with [`Position`](@ref):
``[p, x] = -i``. Returns an [`Op`](@ref) tagged `OP_MOMENTUM`.

See also [`Position`](@ref), [`PhaseSpace`](@ref).
"""
Momentum(name::Union{Symbol, Int32}, si::Int, idx::Index) = Op(OP_MOMENTUM, _name_id(name), si, idx, 0, 0, 0, 0)
Momentum(name::Union{Symbol, Int32}, si::Int) = Momentum(name, si, NO_INDEX)

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
