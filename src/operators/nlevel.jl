"""
    NLevelSpace(name::Symbol, n::Int)
    NLevelSpace(name::Symbol, n::Int, ground_state::Int)
    NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space of a system with `n` discrete levels, home of the [`Transition`](@ref)
operators.

Levels are numbered `1:n`, or named by a tuple of symbols in the order of the basis.
The `ground_state` (default `1`; the first name for named levels) marks the
projector ``|g\\rangle\\langle g|`` that [`expand_completeness`](@ref) can eliminate
through ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|``.

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, (:g, :e))
ℋ(atom)

julia> σ = Transition(h, :σ, :g, :e);

julia> σ' * σ
σ₂₂

julia> h3 = NLevelSpace(:atom, 3, 2);

julia> expand_completeness(Transition(h3, :σ, 2, 2))
1 - σ₁₁ - σ₃₃
```

See also [`Transition`](@ref), [`expand_completeness`](@ref), [`PauliSpace`](@ref).
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
const EMPTY_LEVELS = Symbol[]
NLevelSpace(name::Symbol, n::Int, ground_state::Int) = NLevelSpace(name, n, ground_state, EMPTY_LEVELS)
NLevelSpace(name::Symbol, n::Int) = NLevelSpace(name, n, 1, EMPTY_LEVELS)
function NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})
    return NLevelSpace(name, length(levels), 1, collect(Symbol, levels))
end
NLevelSpace(name::AbstractString, args...) = name_must_be_symbol(name)

Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.name == b.name && a.n == b.n && a.ground_state == b.ground_state && a.levels == b.levels
Base.hash(a::NLevelSpace, h::UInt) = hash(:NLevelSpace, hash(a.name, hash(a.n, hash(a.ground_state, hash(a.levels, h)))))

"""
    CollectiveNLevelSpace(name::Symbol, n::Int)
    CollectiveNLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space of an ensemble of identical `n`-level systems, home of the
[`CollectiveTransition`](@ref) operators.

The number of systems is not part of the space; it is fixed only by the basis chosen for
numerics. Levels are numbered or named as for [`NLevelSpace`](@ref). There is no ground
state to eliminate, since the collective projectors add up to the number of systems rather
than to one.

# Examples
```jldoctest
julia> h = CollectiveNLevelSpace(:atoms, (:g, :e))
ℋ(atoms)

julia> S = CollectiveTransition(h, :S, :e, :g);

julia> commutator(S, S')
-S₁₁ + S₂₂
```

See also [`CollectiveTransition`](@ref), [`NLevelSpace`](@ref).
"""
struct CollectiveNLevelSpace <: HilbertSpace
    name::Symbol
    n::Int
    levels::Vector{Symbol}
    function CollectiveNLevelSpace(name::Symbol, n::Int, levels::Vector{Symbol})
        n >= 1 || throw(ArgumentError("Number of levels must be positive, got $n"))
        isempty(levels) || length(levels) == n || throw(ArgumentError("levels length $(length(levels)) != n=$n"))
        return new(name, n, levels)
    end
end
CollectiveNLevelSpace(name::Symbol, n::Int) = CollectiveNLevelSpace(name, n, EMPTY_LEVELS)
CollectiveNLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}}) =
    CollectiveNLevelSpace(name, length(levels), collect(Symbol, levels))

Base.:(==)(a::CollectiveNLevelSpace, b::CollectiveNLevelSpace) =
    a.name == b.name && a.n == b.n && a.levels == b.levels
Base.hash(a::CollectiveNLevelSpace, h::UInt) =
    hash(:CollectiveNLevelSpace, hash(a.name, hash(a.n, hash(a.levels, h))))

const LevelSpace = Union{NLevelSpace, CollectiveNLevelSpace}

"""
    level_index(h, level::Int) -> Int
    level_index(h, level::Symbol) -> Int

Number of `level` in the [`NLevelSpace`](@ref) or [`CollectiveNLevelSpace`](@ref) `h`. An
integer is checked to lie in `1:n`, a symbol is looked up among the level names.
"""
function level_index(h::LevelSpace, i::Int)
    1 <= i <= h.n || throw(ArgumentError("Level $i out of range 1:$(h.n)"))
    return i
end
function level_index(h::LevelSpace, s::Symbol)
    isempty(h.levels) && throw(ArgumentError("$(nameof(typeof(h))) $(h.name) has no symbolic levels"))
    idx = position_of(==(s), h.levels)
    iszero(idx) && throw(ArgumentError("Level :$s not found in $(h.levels)"))
    return idx
end

"""
    Transition(h::NLevelSpace, name::Symbol, i, j)
    Transition(h::ProductSpace, name::Symbol, i, j, [subsystem::Int])

Transition operator from level `j` to level `i`.

The operator is ``\\sigma^{ij} = |i\\rangle\\langle j|``. The levels are integers or
the symbolic names of the [`NLevelSpace`](@ref); either way it prints with level
numbers.
Transitions on the same system compose,
``\\sigma^{ij}\\sigma^{kl} = \\delta_{jk}\\sigma^{il}``, and the adjoint swaps the
levels. The ground-state projector ``\\sigma^{gg}`` is kept as it is; call
[`expand_completeness`](@ref) to rewrite it in terms of the other projectors.
On a [`ProductSpace`](@ref), `subsystem` selects the factor and may be omitted
when only one [`NLevelSpace`](@ref) is present.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> σ = Transition(h, :σ, 1, 2);

julia> σ * σ', σ' * σ, σ * σ
(σ₁₁, σ₂₂, 0)
```

See also [`NLevelSpace`](@ref), [`expand_completeness`](@ref), [`@qnumbers`](@ref).
"""
Transition(name::Symbol, i::Integer, j::Integer, si::Integer, idx::Index, ground_state::Integer, n_levels::Integer) =
    Op(OP_TRANSITION, name_id(name), si, idx, i, j, ground_state, n_levels)

Transition(h::NLevelSpace, name::Symbol, i::L, j::L) where {L <: Union{Int, Symbol}} =
    Transition(name, level_index(h, i), level_index(h, j), 1, NO_INDEX, h.ground_state, h.n)
Transition(h::ProductSpace, name::Symbol, i::L, j::L, idx::Int = unique_subspace_index(h, NLevelSpace)) where {L <: Union{Int, Symbol}} =
    on_factor(Transition, h, NLevelSpace, idx, name, i, j)

Transition(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)

"""
    CollectiveTransition(h::CollectiveNLevelSpace, name::Symbol, i, j)
    CollectiveTransition(h::ProductSpace, name::Symbol, i, j, [subsystem::Int])

Collective transition from level `j` to level `i`, summed over all systems of the ensemble.

The operator is ``S^{ij} = \\sum_k |i\\rangle_k\\langle j|``. Collective transitions obey
``[S^{ij}, S^{kl}] = \\delta_{jk}S^{il} - \\delta_{li}S^{kj}`` and are normal ordered like
ladder operators. Unlike a single [`Transition`](@ref), a product of two of them does not
compose, since it includes pairs of different systems. On a [`ProductSpace`](@ref),
`subsystem` selects the [`CollectiveNLevelSpace`](@ref) factor and may be omitted when there
is only one.

# Examples
```jldoctest
julia> h = CollectiveNLevelSpace(:atoms, (:g, :e));

julia> @qnumbers S::CollectiveTransition(h, :e, :g);

julia> S' * S
S₁₁ - S₂₂ + S₂₁ * S₁₂
```

See also [`CollectiveNLevelSpace`](@ref), [`Transition`](@ref).
"""
CollectiveTransition(name::Symbol, i::Integer, j::Integer, si::Integer) =
    Op(OP_COLLECTIVE_TRANSITION, name_id(name), si, NO_INDEX, i, j, 0, 0)

CollectiveTransition(h::CollectiveNLevelSpace, name::Symbol, i::L, j::L) where {L <: Union{Int, Symbol}} =
    CollectiveTransition(name, level_index(h, i), level_index(h, j), 1)
CollectiveTransition(h::ProductSpace, name::Symbol, i::L, j::L, idx::Int = unique_subspace_index(h, CollectiveNLevelSpace)) where {L <: Union{Int, Symbol}} =
    on_factor(CollectiveTransition, h, CollectiveNLevelSpace, idx, name, i, j)

CollectiveTransition(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
