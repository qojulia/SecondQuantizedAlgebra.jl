"""
    NLevelSpace(name::Symbol, n::Int)
    NLevelSpace(name::Symbol, n::Int, ground_state::Int)
    NLevelSpace(name::Symbol, levels::Tuple{Vararg{Symbol}})

Hilbert space for an N-level system (atoms, qubits, qudits, etc.).

Supports [`Transition`](@ref) operators ``|i\\rangle\\langle j|`` with the composition
rule ``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``.

The `ground_state` (default `1`) defines the canonical basis as
``\\{|i\\rangle\\langle j| : (i,j) \\neq (g,g)\\} \\cup \\{1\\}``: ground-state projectors
are eliminated via completeness ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g} |k\\rangle\\langle k|``.
Under [`NormalOrder`](@ref) this rewrite fires eagerly during multiplication; under
[`LazyOrder`](@ref) the user opts in via [`simplify`](@ref)`(expr, h)` or
[`normal_order`](@ref)`(expr, h)`. Each [`Transition`](@ref) carries `ground_state`
and `n_levels` directly so the algebra stays Hilbert-space-decoupled.

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

Each `Transition` carries the ground-state level (`ground_state`) and number of levels
(`n_levels`) of its host [`NLevelSpace`](@ref). This lets the [`NormalOrder`](@ref) eager
arithmetic apply the completeness relation ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|``
without consulting the Hilbert space — keeping the algebra Hilbert-space-decoupled
while preserving canonical form.

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
@qnumbers σ::Transition(h, 1, 2)
```
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

# Construction from Hilbert spaces
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

# IndexedOperator convenience
IndexedOperator(op::Transition, i::Index) = Transition(op.name, op.i, op.j, op.space_index, i, op.ground_state, op.n_levels)

# Adjoint: |i⟩⟨j|† = |j⟩⟨i|
Base.adjoint(op::Transition) = Transition(op.name, op.j, op.i, op.space_index, op.index, op.ground_state, op.n_levels)

# Equality
Base.isequal(a::Transition, b::Transition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.index == b.index && a.ground_state == b.ground_state && a.n_levels == b.n_levels
Base.:(==)(a::Transition, b::Transition) = isequal(a, b)

# Hashing
Base.hash(a::Transition, h::UInt) = hash(:Transition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.index, hash(a.ground_state, hash(a.n_levels, h))))))))

# Ladder (not applicable to Transition)
ladder(::Transition) = 0

"""
    CollectiveTransition <: QSym

Collective transition operator ``S^{ij} = \\sum_{k=1}^{N} \\sigma_k^{ij}`` on an
[`NLevelSpace`](@ref). Lives in the *symmetric (bosonic) subspace* of `N`
indistinguishable atoms — the natural representation for Dicke-style physics
where atoms are not individually addressable and the dynamics preserves
permutation symmetry.

# When to use this vs indexed `Σ`

Two regimes — pick the one that matches your physics, do **not** mix:

| | Distinguishable atoms | Indistinguishable atoms |
|---|---|---|
| Representation | [`Σ`](@ref) of [`IndexedOperator`](@ref) on a [`ProductSpace`](@ref) | `CollectiveTransition` on a single [`NLevelSpace`](@ref) |
| Physical setting | Atom-resolved coupling, disorder, individual addressing | Dicke physics, single-mode driving, no disorder |
| Hilbert space | Full ``n^N`` (one atom per subspace) | Symmetric subspace, dimension ``\\binom{N+n-1}{n-1}`` |
| Numeric scaling | Exponential in `N` | Polynomial in `N` — tractable for thousands |

For the same physics where both are valid (e.g. Tavis–Cummings with identical
atoms), the indexed form is preferred; `CollectiveTransition` is for problems
where the symmetric-subspace assumption is essential, typically because the full
Hilbert space is intractable.

# Algebra

The operators satisfy the ``\\mathfrak{su}(N)`` Lie algebra
``[S^{ij}, S^{kl}] = \\delta_{jk} S^{il} - \\delta_{li} S^{kj}``. Under
[`NormalOrder`](@ref) this commutator fires eagerly during multiplication: the
operator with the larger `(i, j)` (lex-descending) ends up on the left, with
the appropriate remainder terms appended.

The single-atom completeness ``\\sigma^{gg} = 1 - \\sum_{k \\neq g}\\sigma^{kk}``
does **not** apply here — for collective operators ``\\sum_j S^{jj} = N \\cdot I``
(note the ``N``), and that rewrite is intentionally left for numeric-conversion
time rather than the symbolic algebra.

# Construction

```julia
h = NLevelSpace(:atom, 3)
S(i, j) = CollectiveTransition(h, :S, i, j)
S(2, 1) * S(1, 3)         # ordered — stays as written
S(1, 2) * S(2, 3)         # unordered — eagerly normal-ordered: S(1, 3) + S(2, 3) * S(1, 2)
```

In a [`ProductSpace`](@ref) supply the subspace index:
```julia
hp = FockSpace(:c) ⊗ NLevelSpace(:a, 3)
S = CollectiveTransition(hp, :S, 1, 2, 2)
```

# Numeric conversion

Numeric simulation requires a `QuantumOpticsBase.ManyBodyBasis` built over the
single-atom basis with `bosonstates`:
```julia
b1 = NLevelBasis(3)
b  = ManyBodyBasis(b1, bosonstates(b1, N))    # symmetric subspace for N atoms
to_numeric(S(1, 2), b)                        # many-body operator
```
"""
struct CollectiveTransition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    index::Index   # always NO_INDEX; carried for site-sort compatibility
end

# Construction from Hilbert spaces
function CollectiveTransition(h::NLevelSpace, name::Symbol, i::Int, j::Int)
    1 <= i <= h.n || throw(ArgumentError("Level i=$i out of range 1:$(h.n)"))
    1 <= j <= h.n || throw(ArgumentError("Level j=$j out of range 1:$(h.n)"))
    return CollectiveTransition(name, i, j, 1, NO_INDEX)
end
function CollectiveTransition(h::NLevelSpace, name::Symbol, i::Symbol, j::Symbol)
    return CollectiveTransition(name, _level_index(h, i), _level_index(h, j), 1, NO_INDEX)
end
function CollectiveTransition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    1 <= i <= space.n || throw(ArgumentError("Level i=$i out of range 1:$(space.n)"))
    1 <= j <= space.n || throw(ArgumentError("Level j=$j out of range 1:$(space.n)"))
    return CollectiveTransition(name, i, j, idx, NO_INDEX)
end
function CollectiveTransition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    space = h.spaces[idx]
    space isa NLevelSpace || throw(ArgumentError("Space at index $idx is not an NLevelSpace"))
    return CollectiveTransition(name, _level_index(space, i), _level_index(space, j), idx, NO_INDEX)
end

# Auto-detect subspace when the ProductSpace contains exactly one NLevelSpace.
CollectiveTransition(h::ProductSpace, name::Symbol, i::Int, j::Int) =
    CollectiveTransition(h, name, i, j, _unique_subspace_index(h, NLevelSpace))
CollectiveTransition(h::ProductSpace, name::Symbol, i::Symbol, j::Symbol) =
    CollectiveTransition(h, name, i, j, _unique_subspace_index(h, NLevelSpace))

# IndexedOperator does not apply — CollectiveTransition is already a sum over atoms.
function IndexedOperator(::CollectiveTransition, ::Index)
    throw(ArgumentError("CollectiveTransition is already a sum over atoms; indexing it is not meaningful. Use Transition for individually-addressed atoms."))
end

# Adjoint: (S^{ij})† = S^{ji}
Base.adjoint(op::CollectiveTransition) = CollectiveTransition(op.name, op.j, op.i, op.space_index, op.index)

# Equality and hashing
Base.isequal(a::CollectiveTransition, b::CollectiveTransition) = a.name == b.name && a.i == b.i && a.j == b.j && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::CollectiveTransition, b::CollectiveTransition) = isequal(a, b)
Base.hash(a::CollectiveTransition, h::UInt) = hash(:CollectiveTransition, hash(a.name, hash(a.i, hash(a.j, hash(a.space_index, hash(a.index, h))))))

# Ladder (not applicable)
ladder(::CollectiveTransition) = 0
