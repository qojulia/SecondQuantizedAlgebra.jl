const INTERN_LOCK = ReentrantLock()

const NAME_BY_ID = Symbol[]
const NAME_TO_ID = Dict{Symbol, Int32}()
const NAME_RANK = Int32[]
const SYM_BY_ID = Num[]
const UNBUILT_SYM = Num(SymbolicUtils.Sym{SymbolicUtils.SymReal}(:__unbuilt_index_sym__; type = Int))

name_from_id(id::Int32) = id == 0 ? :_ : NAME_BY_ID[id]
name_rank(id::Int32) = id == 0 ? Int32(0) : NAME_RANK[id]

function base_sym_from_id(id::Int32)
    if 1 <= id <= length(SYM_BY_ID)
        cached = SYM_BY_ID[id]
        cached === UNBUILT_SYM || return cached
    end
    return build_base_sym(id)
end

@noinline function build_base_sym(id::Int32)
    return lock(INTERN_LOCK) do
        while length(SYM_BY_ID) < length(NAME_BY_ID)
            push!(SYM_BY_ID, UNBUILT_SYM)
        end
        cached = SYM_BY_ID[id]
        cached === UNBUILT_SYM || return cached
        s = Num(SymbolicUtils.Sym{SymbolicUtils.SymReal}(name_from_id(id); type = Int))
        SYM_BY_ID[id] = s
        return s
    end
end

function intern!(by_id::Vector{T}, to_id::Dict{T, Int32}, x::T, on_new::F) where {T, F}
    id = get(to_id, x, Int32(0))
    id != 0 && return id
    lock(INTERN_LOCK)
    try
        id = get(to_id, x, Int32(0))
        id != 0 && return id
        push!(by_id, x)
        id = Int32(length(by_id))
        to_id[x] = id
        on_new(id)
        return id
    finally
        unlock(INTERN_LOCK)
    end
end

function new_name!(id::Int32)
    push!(SYM_BY_ID, UNBUILT_SYM)
    s = NAME_BY_ID[id]
    r = Int32(1)
    for j in 1:(id - 1)
        NAME_BY_ID[j] < s ? (r += Int32(1)) : (NAME_RANK[j] += Int32(1))
    end
    push!(NAME_RANK, r)
    return id
end

intern_name(s::Symbol) = intern!(NAME_BY_ID, NAME_TO_ID, s, new_name!)
name_id(s::Symbol) = intern_name(s)

const RANGE_BY_ID = Num[]
const RANGE_TO_ID = Dict{Num, Int32}()

range_from_id(id::Int32) = id == 0 ? Num(0) : RANGE_BY_ID[id]

no_init(id::Int32) = id
intern_range(r::Num) = intern!(RANGE_BY_ID, RANGE_TO_ID, r, no_init)
intern_range(r::Integer) = intern_range(Num(r))
