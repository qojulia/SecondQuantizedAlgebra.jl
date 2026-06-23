# Global intern tables backing the isbits `Op`/`Index` (issue #137).
#
# Operator/index names and index ranges are stored as `Int32` ids so that `Op`
# and `Index` become `isbits` (dense inline storage, integer hashing, no GC
# scanning of the operator vector). The id to object direction is a never-
# shrinking, never-relocating `Vector` (lock-free `@inbounds` reads on the hot
# path); the object to id direction is a `Dict` mutated only at construction
# (cold path) under `_INTERN_LOCK`.
#
# `_NAME_RANK` gives each name id its lexicographic position so canonical
# ordering stays alphabetical and deterministic across sessions even though the
# ids themselves are insertion-order. Ids are NEVER serialized (not portable
# across sessions); `Index`/`Op` are not `Serialization`-stable.

const _INTERN_LOCK = ReentrantLock()

# Name interning.
const _NAME_BY_ID = Symbol[]              # id (1-based) -> Symbol
const _NAME_TO_ID = Dict{Symbol, Int32}() # Symbol -> id
const _NAME_RANK = Int32[]                # id -> lexicographic position

# id 0 is the reserved sentinel (NO_INDEX / no name).
_name_from_id(id::Int32)::Symbol = id == 0 ? :_ : @inbounds _NAME_BY_ID[id]
_name_rank(id::Int32)::Int32 = id == 0 ? Int32(0) : @inbounds _NAME_RANK[id]

function _recompute_name_ranks!()
    n = length(_NAME_BY_ID)
    resize!(_NAME_RANK, n)
    order = sortperm(_NAME_BY_ID)         # order[r] == id whose name sorts at rank r
    for (r, id) in enumerate(order)
        _NAME_RANK[id] = Int32(r)
    end
    return nothing
end

function _intern_name(s::Symbol)::Int32
    return lock(_INTERN_LOCK) do
        id = get(_NAME_TO_ID, s, Int32(0))
        id != 0 && return id
        push!(_NAME_BY_ID, s)
        id = Int32(length(_NAME_BY_ID))
        _NAME_TO_ID[s] = id
        _recompute_name_ranks!()          # O(n log n) over distinct names; cold path
        return id
    end
end

# Constructor name argument: a `Symbol` is interned; an already-interned `Int32`
# id is passed through (internal rebuilds forward ids without re-interning).
_name_id(s::Symbol)::Int32 = _intern_name(s)
_name_id(id::Int32)::Int32 = id

# Range interning.
const _RANGE_BY_ID = Num[]               # id (1-based) -> Num
const _RANGE_TO_ID = Dict{Num, Int32}()  # Num -> id (dedup by isequal)

_range_from_id(id::Int32)::Num = id == 0 ? Num(0) : @inbounds _RANGE_BY_ID[id]

function _intern_range(r::Num)::Int32
    return lock(_INTERN_LOCK) do
        id = get(_RANGE_TO_ID, r, Int32(0))
        id != 0 && return id
        push!(_RANGE_BY_ID, r)
        id = Int32(length(_RANGE_BY_ID))
        _RANGE_TO_ID[r] = id
        return id
    end
end
_intern_range(r::Integer) = _intern_range(Num(r))
