# Global intern tables backing the isbits `Op`/`Index` (issue #137).
#
# Names/ranges are stored as `Int32` ids so `Op`/`Index` are `isbits`. The
# id->object tables are append-only `Vector`s (bounds-checked reads on the hot
# path); the object->id tables are `Dict`s mutated only at construction under
# `_INTERN_LOCK`. `_NAME_RANK` holds each id's lexicographic position so
# canonical ordering stays name-sorted and deterministic across sessions.
#
# Ids are insertion-order: NEVER serialized (`Index`/`Op` are not
# `Serialization`-stable; reads bounds-check so a stale id throws, not corrupts).
# Construction is the only writer and is NOT thread-safe; populate the tables
# from one thread before any concurrent canonicalization.

const _INTERN_LOCK = ReentrantLock()

# Name interning.
const _NAME_BY_ID = Symbol[]              # id (1-based) -> Symbol
const _NAME_TO_ID = Dict{Symbol, Int32}() # Symbol -> id
const _NAME_RANK = Int32[]                # id -> lexicographic position
const _SYM_BY_ID = Union{Nothing, Num}[]  # id -> cached Num(Sym(name)); lazy (see _base_sym_from_id)

# id 0 is the reserved sentinel (NO_INDEX / no name).
_name_from_id(id::Int32)::Symbol = id == 0 ? :_ : _NAME_BY_ID[id]
_name_rank(id::Int32)::Int32 = id == 0 ? Int32(0) : _NAME_RANK[id]

function _base_sym_from_id(id::Int32)::Num
    if 1 <= id <= length(_SYM_BY_ID)
        cached = _SYM_BY_ID[id]
        cached === nothing || return cached
    end
    return _build_base_sym(id)
end

@noinline function _build_base_sym(id::Int32)::Num
    return lock(_INTERN_LOCK) do
        while length(_SYM_BY_ID) < length(_NAME_BY_ID)  # defensive lockstep
            push!(_SYM_BY_ID, nothing)
        end
        cached = _SYM_BY_ID[id]
        cached === nothing || return cached
        s = Num(SymbolicUtils.Sym{SymbolicUtils.SymReal}(_name_from_id(id); type = Int))
        _SYM_BY_ID[id] = s
        return s
    end
end

# get-or-push with a lock-free fast path for the common already-interned case;
# `on_new(id)` runs under the lock for a freshly assigned id. `F` is a sink
# type-parameter so the callback specializes; explicit lock/try (not a `do`
# closure) keeps the fast path allocation-free.
function _intern!(by_id::Vector{T}, to_id::Dict{T, Int32}, x::T, on_new::F)::Int32 where {T, F}
    id = get(to_id, x, Int32(0))
    id != 0 && return id
    lock(_INTERN_LOCK)
    try
        id = get(to_id, x, Int32(0))
        id != 0 && return id
        push!(by_id, x)
        id = Int32(length(by_id))
        to_id[x] = id
        on_new(id)
        return id
    finally
        unlock(_INTERN_LOCK)
    end
end

# Freshly interned name: add its lazy `index_sym` slot, then splice its
# lexicographic rank in place (shift later ranks by one) instead of re-sorting.
function _new_name!(id::Int32)
    push!(_SYM_BY_ID, nothing)
    s = _NAME_BY_ID[id]
    r = Int32(1)
    for j in 1:(id - 1)
        _NAME_BY_ID[j] < s ? (r += Int32(1)) : (_NAME_RANK[j] += Int32(1))
    end
    push!(_NAME_RANK, r)
    return nothing
end

_intern_name(s::Symbol)::Int32 = _intern!(_NAME_BY_ID, _NAME_TO_ID, s, _new_name!)
# Constructor name argument: always a `Symbol`, interned to its id. Internal
# rebuilds forward an existing `name_id` by constructing `Op` directly.
_name_id(s::Symbol)::Int32 = _intern_name(s)

# Range interning.
const _RANGE_BY_ID = Num[]               # id (1-based) -> Num
const _RANGE_TO_ID = Dict{Num, Int32}()  # Num -> id (dedup by isequal)

_range_from_id(id::Int32)::Num = id == 0 ? Num(0) : _RANGE_BY_ID[id]

_no_init(::Int32) = nothing
_intern_range(r::Num)::Int32 = _intern!(_RANGE_BY_ID, _RANGE_TO_ID, r, _no_init)
_intern_range(r::Integer) = _intern_range(Num(r))
