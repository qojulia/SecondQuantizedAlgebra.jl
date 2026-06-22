"""
    Index(h::HilbertSpace, name::Symbol, range, space)
    Index(h::HilbertSpace, name::Symbol, range, space_index::Int)

Symbolic summation index for site- or mode-resolved operator expressions.

Represents labels such as `i`, `j`, or `k` that attach to operators and
parameters. Pass a subspace type or integer position to specify which factor of
a [`ProductSpace`](@ref) the index ranges over. Use with [`IndexedOperator`](@ref)
to build objects like ``a_i`` and with [`Σ`](@ref) to build sums
``\\sum_i a_i^\\dagger a_i``.

# Examples
```jldoctest
julia> h = FockSpace(:site) ⊗ FockSpace(:cavity);

julia> i = Index(h, :i, 10, FockSpace(:site));

julia> j = Index(h, :j, 10, 1); # same subspace via integer position

julia> i == j
false
```

See also [`has_index`](@ref), [`IndexedOperator`](@ref), [`Σ`](@ref),
[`change_index`](@ref).
"""
struct Index
    name::Symbol
    range::Num
    space_index::Int
    sym::Num
end

const NO_INDEX = Index(:_, Num(0), 0, Num(0))

"""
Metadata key carrying the concrete slot integer on a per-slot index sym minted by
`(i::Index)(k)`. Sym equality/hashing ignore metadata, so a slot-stamped sym still
dedup-equals its name-only counterpart.
"""
struct IndexSlot end

"""
    has_index(idx::Index) -> Bool

Return `true` when `idx` is an actual symbolic index.

Returns `false` for the sentinel `NO_INDEX`, used for operators without an
attached index.
"""
has_index(idx::Index) = idx.space_index != 0

# A pair of indices `(α, β)` meaning `α ≠ β`.
const NonEqualPair = Tuple{Index, Index}

# Shared sentinels for empty vectors on hot paths. Never mutated.
# `_EMPTY_OPS` lives in `operators/op.jl` (needs `Op`, defined after `Index`).
const _EMPTY_NE = NonEqualPair[]
const _EMPTY_INDICES = Index[]

function Base.:(==)(a::Index, b::Index)
    a === b && return true
    return a.name == b.name && a.space_index == b.space_index && isequal(a.range, b.range)
end
function Base.hash(a::Index, h::UInt)
    return hash(:Index, hash(a.name, hash(a.space_index, h)))
end
Base.isless(a::Index, b::Index) = (a.space_index, a.name) < (b.space_index, b.name)

# Ordering key for an index; mirrors isless's (space_index, name) (range excluded).
_index_key(idx::Index) = (idx.space_index, idx.name)

function Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, space::HilbertSpace)
    si = _find_space_index(h, space)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, Num(range), si, Num(sym_var))
end
function Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, si::Int)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, Num(range), si, Num(sym_var))
end

"""
    (i::Index)(k::Integer) -> Index

Return a fresh per-slot `Index` for position `k`, named `Symbol(i.name, "_", k)`
and inheriting `i.range` and `i.space_index`.

Use as a slot template after `evaluate` unrolls a sum over `i`: `i(3)` is the
concrete index for the third atom, matching the indices that `evaluate` mints
internally so the resulting operators dedup-equal `evaluate`'s output. Naming
seeds from `i.name`, so the user's vocabulary is preserved (SQA's naming policy).
"""
function (i::Index)(k::Integer)
    name = Symbol(i.name, "_", k)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    sym_var = SymbolicUtils.setmetadata(sym_var, IndexSlot, k)
    return Index(name, i.range, i.space_index, Num(sym_var))
end

"""
    index_slot(x) -> Union{Int, Nothing}

Concrete slot integer of a per-slot index sym minted by an [`Index`](@ref) call `(i::Index)(k)`,
or `nothing` for any other symbol. Lets consumers recover the position `k` without
parsing the symbol's name.
"""
index_slot(x) = SymbolicUtils.getmetadata(SymbolicUtils.unwrap(x), IndexSlot, nothing)

_find_space_index(::HilbertSpace, ::HilbertSpace) = 1
function _find_space_index(h::ProductSpace, space::HilbertSpace)
    for (i, s) in enumerate(h.spaces)
        s == space && return i
    end
    throw(ArgumentError("Space $space not found in $h"))
end
