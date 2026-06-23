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
    name_id::Int32       # interned name (see intern.jl); 0 == anonymous concrete site / NO_INDEX
    range_id::Int32      # interned range Num; 0 == no range
    space_index::Int32
    slot::Int32          # 0 == abstract; k>0 == concrete site k
end

const NO_INDEX = Index(Int32(0), Int32(0), Int32(0), Int32(0))

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

# Equality compares the interned ids (pure integers). `slot` is excluded, exactly
# as the old equality excluded the per-slot `sym`: concrete sites are distinguished
# out-of-band (numeric reads `slot` via `index_slot`), never via `Index ==`.
function Base.:(==)(a::Index, b::Index)
    a === b && return true
    return a.name_id == b.name_id && a.space_index == b.space_index && a.range_id == b.range_id
end
function Base.hash(a::Index, h::UInt)
    return hash(:Index, hash(a.name_id, hash(a.space_index, h)))
end
# Order by lexicographic name rank (not raw id, which is insertion-order), so
# canonical form stays alphabetical and deterministic across sessions.
Base.isless(a::Index, b::Index) = (a.space_index, _name_rank(a.name_id)) < (b.space_index, _name_rank(b.name_id))

# Ordering key for an index; mirrors isless's (space_index, name rank) (range excluded).
_index_key(idx::Index) = (idx.space_index, _name_rank(idx.name_id))

function Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, space::HilbertSpace)
    si = _find_space_index(h, space)
    return Index(_intern_name(name), _intern_range(range), Int32(si), Int32(0))
end
function Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, si::Int)
    return Index(_intern_name(name), _intern_range(range), Int32(si), Int32(0))
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
    name = Symbol(_name_from_id(i.name_id), "_", k)
    return Index(_intern_name(name), i.range_id, i.space_index, Int32(k))
end

"""
    index_slot(x) -> Union{Int, Nothing}
    index_slot(idx::Index) -> Union{Int, Nothing}

Concrete slot integer of a per-slot index, or `nothing` for an abstract index.
On an [`Index`](@ref) reads the `slot` field directly; on a symbol reads the
`IndexSlot` metadata minted by `(i::Index)(k)`, so consumers can recover the
position `k` without parsing the symbol's name.
"""
index_slot(idx::Index)::Union{Int, Nothing} = idx.slot == 0 ? nothing : Int(idx.slot)
index_slot(x) = SymbolicUtils.getmetadata(SymbolicUtils.unwrap(x), IndexSlot, nothing)

"""
    index_range(idx::Index) -> Num

The summation range of `idx` (the user's `Num`, recovered from the intern table).
Returns `Num(0)` for the sentinel `NO_INDEX`.
"""
index_range(idx::Index)::Num = _range_from_id(idx.range_id)

"""
    index_name(idx::Index) -> Symbol

The display name of `idx` (recovered from the intern table).
"""
index_name(idx::Index)::Symbol = _name_from_id(idx.name_id)

"""
    index_sym(idx::Index) -> Num

The symbolic variable for `idx`, reconstructed from its interned name (plus the
`IndexSlot` metadata for a per-slot index). SymbolicUtils hashconsing makes the
reconstruction identical (`===`) to the originally minted symbol, so substitution
and `get_variables` are unaffected. An anonymous concrete site (`name_id == 0`,
`slot == k`) reconstructs to the integer `Num(k)`, matching `to_numeric`'s
resolved-site convention.
"""
function index_sym(idx::Index)::Num
    idx.name_id == 0 && return Num(Int(idx.slot))
    base = _base_sym_from_id(idx.name_id)            # cached name-only Sym (hot path: cached read)
    idx.slot == 0 && return base
    return Num(SymbolicUtils.setmetadata(SymbolicUtils.unwrap(base), IndexSlot, Int(idx.slot)))
end

_find_space_index(::HilbertSpace, ::HilbertSpace) = 1
function _find_space_index(h::ProductSpace, space::HilbertSpace)
    for (i, s) in enumerate(h.spaces)
        s == space && return i
    end
    throw(ArgumentError("Space $space not found in $h"))
end
