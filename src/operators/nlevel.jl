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
Satisfies the composition rule
``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``,
with adjoint ``|i\\rangle\\langle j|^\\dagger = |j\\rangle\\langle i|``.

Each `Transition` carries the ground-state level and number of levels of its
host [`NLevelSpace`](@ref) so the eager arithmetic can keep
``\\sigma^{gg}`` atomic in canonical form (use
[`expand_completeness`](@ref) to materialise the identity
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|`` when
needed).

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
struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    index::Index
    ground_state::Int
    n_levels::Int
end

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

IndexedOperator(op::Transition, i::Index) = Transition(op.name, op.i, op.j, op.space_index, i, op.ground_state, op.n_levels)

Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index, op.index, op.ground_state, op.n_levels)

Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.index == b.index && a.ground_state == b.ground_state && a.n_levels == b.n_levels
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.index, hash(a.ground_state, hash(a.n_levels, h))))))))

ladder(::Transition) = 0

# --- Operator hooks ---

function _site_compare(a::Transition, b::Transition, ne::Vector{NonEqualPair})
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

# Transitions on the same site never commute freely (composition fires).
_can_commute(a::Transition, b::Transition) = false

# Composition: σⁱʲ · σᵏˡ = δⱼₖ σⁱˡ.
function _reduce_pair(a::Transition, b::Transition)
    a.name == b.name || return (NoReduction, a, _CNUM_ZERO)
    a.space_index == b.space_index || return (NoReduction, a, _CNUM_ZERO)
    a.index == b.index || return (NoReduction, a, _CNUM_ZERO)
    if a.j == b.i
        new = Transition(a.name, a.i, b.j, a.space_index, a.index, a.ground_state, a.n_levels)
        return (OpReduction, new, _CNUM_ONE)
    else
        return (ScalarReduction, a, _CNUM_ZERO)    # δ_{j,k} = 0
    end
end

_commute_pair(a::Transition, b::Transition) = (b, a, _CNUM_ZERO, _EMPTY_OPS)

# Ground-state expansion: σᵍᵍ -> 1 - Σ_{k≠g} σᵏᵏ.
function _ground_state_expand(op::Transition)
    op.i == op.ground_state && op.j == op.ground_state || return (false, 0, 0, 0)
    return (true, op.ground_state, op.n_levels, op.space_index)
end
