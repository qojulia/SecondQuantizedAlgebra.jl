"""
    NLevelSpace(name::Symbol, n::Int)
    NLevelSpace(name::Symbol, n::Int, ground_state::Int)
    NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space for an N-level system (atoms, qubits, qudits). Hosts
[`Transition`](@ref) operators
``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``.

The `ground_state` (default `1`) selects the projector ``|g\\rangle\\langle g|``
that the canonical basis omits; arithmetic keeps it atomic and
[`expand_completeness`](@ref) materialises
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|`` on
demand. Levels can be integers (default `1:n`) or symbolic names.

# Examples

```jldoctest
julia> NLevelSpace(:atom, 3)
ℋ(atom)

julia> NLevelSpace(:atom, (:g, :e))
ℋ(atom)
```

See also [`Transition`](@ref), [`expand_completeness`](@ref).
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
    return NLevelSpace(name, length(levels), 1, collect(Symbol, levels))
end
NLevelSpace(name::AbstractString, args...) = _name_must_be_symbol(name)

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
    CollectiveNLevelSpace(name::Symbol, n::Int)
    CollectiveNLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space for collective transitions of identical N-level systems. It hosts
[`CollectiveTransition`](@ref) operators and deliberately carries no ground-state
projector convention: collectively, the diagonal sum is the particle-number
operator rather than the identity.
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
CollectiveNLevelSpace(name::Symbol, n::Int) = CollectiveNLevelSpace(name, n, _EMPTY_LEVELS)
CollectiveNLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}}) =
    CollectiveNLevelSpace(name, length(levels), collect(Symbol, levels))

Base.:(==)(a::CollectiveNLevelSpace, b::CollectiveNLevelSpace) =
    a.name == b.name && a.n == b.n && a.levels == b.levels
Base.hash(a::CollectiveNLevelSpace, h::UInt) =
    hash(:CollectiveNLevelSpace, hash(a.name, hash(a.n, hash(a.levels, h))))

function _level_index(h::CollectiveNLevelSpace, s::Symbol)
    isempty(h.levels) && throw(ArgumentError("CollectiveNLevelSpace $(h.name) has no symbolic levels"))
    idx = findfirst(==(s), h.levels)
    idx === nothing && throw(ArgumentError("Level :$s not found in $(h.levels)"))
    return idx
end

"""
    Transition(h::NLevelSpace, name::Symbol, i, j) -> Op

Transition operator ``|i\\rangle\\langle j|`` on an [`NLevelSpace`](@ref).
Satisfies the composition rule
``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``,
with adjoint ``|i\\rangle\\langle j|^\\dagger = |j\\rangle\\langle i|``. Returns an
[`Op`](@ref) tagged `OP_TRANSITION`, packing `i`, `j`, the ground state, and the
number of levels into `l1`, `l2`, `g`, `nlev`.

The carried ground-state level and level count let the eager arithmetic keep
``\\sigma^{gg}`` atomic in canonical form (use [`expand_completeness`](@ref) to
materialise ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|``
when needed).

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, 3);

julia> σ = Transition(h, :σ, 1, 2);

julia> σ * σ'
σ₁₁

julia> σ' * σ
σ₂₂
```

See also [`NLevelSpace`](@ref), [`expand_completeness`](@ref), [`@qnumbers`](@ref).
"""
Transition(name::Symbol, i::Integer, j::Integer, si::Integer, idx::Index, ground_state::Integer, n_levels::Integer) =
    Op(OP_TRANSITION, _name_id(name), si, idx, i, j, ground_state, n_levels)

function Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return Transition(name, i, j, 1, NO_INDEX, h.ground_state, h.n)
end
function Transition(h::NLevelSpace, name::Symbol, i::Symbol, j::Symbol)
    return Transition(name, _level_index(h, i), _level_index(h, j), 1, NO_INDEX, h.ground_state, h.n)
end
function Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return Transition(name, i, j, idx, NO_INDEX, space.ground_state, space.n)
end
function Transition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    return Transition(name, _level_index(space, i), _level_index(space, j), idx, NO_INDEX, space.ground_state, space.n)
end

# Auto-detect subspace when the ProductSpace contains exactly one NLevelSpace.
Transition(h::ProductSpace, name::Symbol, i::Int, j::Int) =
    Transition(h, name, i, j, _unique_subspace_index(h, NLevelSpace))
Transition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol) =
    Transition(h, name, i, j, _unique_subspace_index(h, NLevelSpace))

Transition(::HilbertSpace, name::AbstractString, args...) = _name_must_be_symbol(name)

"""
    CollectiveTransition(h::CollectiveNLevelSpace, name::Symbol, i, j) -> Op

Collective transition ``S^{ij} = \\sum_k |i\\rangle_k\\langle j|``. These operators
obey ``[S^{ij},S^{kl}] = \\delta_{jk}S^{il}-\\delta_{li}S^{kj}`` and do not compose
like single-site [`Transition`](@ref) operators.
"""
CollectiveTransition(name::Symbol, i::Integer, j::Integer, si::Integer) =
    Op(OP_COLLECTIVE_TRANSITION, _name_id(name), si, NO_INDEX, i, j, 0, 0)

function CollectiveTransition(h::CollectiveNLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return CollectiveTransition(name, i, j, 1)
end
CollectiveTransition(h::CollectiveNLevelSpace, name::Symbol, i::Symbol, j::Symbol) =
    CollectiveTransition(name, _level_index(h, i), _level_index(h, j), 1)

function CollectiveTransition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa CollectiveNLevelSpace || throw(ArgumentError("Space at index $idx is not a CollectiveNLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return CollectiveTransition(name, i, j, idx)
end
function CollectiveTransition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa CollectiveNLevelSpace || throw(ArgumentError("Space at index $idx is not a CollectiveNLevelSpace"))
    return CollectiveTransition(name, _level_index(space, i), _level_index(space, j), idx)
end

CollectiveTransition(h::ProductSpace, name::Symbol, i::Int, j::Int) =
    CollectiveTransition(h, name, i, j, _unique_subspace_index(h, CollectiveNLevelSpace))
CollectiveTransition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol) =
    CollectiveTransition(h, name, i, j, _unique_subspace_index(h, CollectiveNLevelSpace))

CollectiveTransition(::HilbertSpace, name::AbstractString, args...) = _name_must_be_symbol(name)  
