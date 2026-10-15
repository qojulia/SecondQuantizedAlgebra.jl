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
Position(name::Symbol, si::Int) = Op(OP_POSITION, name_id(name), si, NO_INDEX, 0, 0, 0, 0)

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
Momentum(name::Symbol, si::Int) = Op(OP_MOMENTUM, name_id(name), si, NO_INDEX, 0, 0, 0, 0)

Position(h::PhaseSpace, name::Symbol) = Position(name, 1)
Momentum(h::PhaseSpace, name::Symbol) = Momentum(name, 1)

Position(h::ProductSpace, name::Symbol, idx::Int = unique_subspace_index(h, PhaseSpace)) =
    on_factor(Position, h, PhaseSpace, idx, name)
Momentum(h::ProductSpace, name::Symbol, idx::Int = unique_subspace_index(h, PhaseSpace)) =
    on_factor(Momentum, h, PhaseSpace, idx, name)

Position(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
Momentum(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
