# In-place accumulator for additive reductions. Folding many `QAdd`s with
# `Base.:+` is O(n²) because every `+` copies the whole backing dict; the
# builder threads one shared dict through the whole reduction instead. See
# devdocs for the rationale and why the product path is intentionally left as
# repeated `*`.
struct QAddBuilder
    args::QTermDict
    indices::Vector{Index}
end
QAddBuilder() = QAddBuilder(QTermDict(), Index[])

# In-place dedup-append union of `src` into `dest` (the mutating analog of
# `merge_unique`); keeps the builder's identity stable.
function merge_indices!(dest::Vector{Index}, src::Vector{Index})
    for x in src
        x ∉ dest && push!(dest, x)
    end
    return dest
end

# One accumulate primitive per argument shape, reusing the exact in-place ops of
# the corresponding `Base.:+` method so results match a `+`-chain byte-for-byte.
function accumulate!(b::QAddBuilder, x::QAdd)
    for (term, c) in x.arguments
        addto_key!(b.args, copy_key(term), c)
    end
    isempty(x.indices) || merge_indices!(b.indices, x.indices)
    return b
end
accumulate!(b::QAddBuilder, x::QSym) = (addto!(b.args, Op[x], CNUM_ONE); b)
function accumulate!(b::QAddBuilder, x::Coefficient)
    x isa Number && iszero(x) && return b
    addto!(b.args, EMPTY_OPS, to_cnum(x))
    return b
end

# Materialize once: the single `prune_dead_ne` (constructor invariant) and
# `drop_unused_indices` are equivalent to doing them per `+` step because
# mid-chain index drops only trim the index vector, never delete terms. Returns
# the shared `ZERO_QADD` for the empty case without mutating it.
function build(b::QAddBuilder)
    isempty(b.args) && return zero_qadd()
    return QAdd(b.args, drop_unused_indices(b.args, b.indices))
end

function MA.operate(
        ::typeof(sum), a::AbstractArray{<:QAdd}; init = zero(QAdd),
    )
    isempty(a) && return init
    if init isa QAdd || init isa QSym || init isa Coefficient
        b = QAddBuilder()
        accumulate!(b, init)
        for x in a
            accumulate!(b, x)
        end
        return build(b)
    end
    return mapreduce(identity, MA.add!!, a; init)
end
MA.promote_operation(::typeof(sum), ::Type{<:AbstractArray{<:QAdd}}) = QAdd

# Ergonomic entry points: only `AbstractArray{<:QAdd}` (element is a `QAdd`, so
# the result is a `QAdd`). A bracketed comprehension `sum([… for …])` is the fast
# path; a bare generator stays on Base's generic `+`-fold.
function Base.sum(a::AbstractArray{<:QAdd})
    isempty(a) && return zero(QAdd)
    b = QAddBuilder()
    for x in a
        accumulate!(b, x)
    end
    return build(b)
end
function Base.sum(f, a::AbstractArray{<:QAdd})
    isempty(a) && return zero(QAdd)
    b = QAddBuilder()
    for x in a
        accumulate!(b, f(x))
    end
    return build(b)
end
# Narrowed to the concrete `Vector` (what comprehensions produce) rather than
# `AbstractArray` to avoid a method ambiguity with `StaticArrays.reduce`.
Base.reduce(::typeof(+), a::Vector{<:QAdd}) = sum(a)
