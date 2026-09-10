"""
    QAdd <: QField

The sole compound expression type: a sum of eagerly-ordered operator products.

All arithmetic on [`QSym`](@ref) operators produces a `QAdd`. Iterating over a
`QAdd` yields `Pair{QTerm, CNum}` entries; read `term.ops` for the operator
sequence and `term.ne` for the scoped index constraints.

See also [`QTerm`](@ref), [`get_prefactor`](@ref), [`get_operators`](@ref), [`Σ`](@ref),
[`constraint_pairs`](@ref).
"""
struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
    function QAdd(arguments::QTermDict, indices::Vector{Index})
        return new(prune_dead_ne(arguments, indices), indices)
    end
end

# Strip NE pairs that reference no op index, coefficient index, or scope index.
function prune_dead_ne(args::QTermDict, indices::Vector{Index})
    needs_rebuild = false
    for (term, c) in args
        for p in term.ne
            if !pair_referenced(p, term.ops, c, indices)
                needs_rebuild = true
                break
            end
        end
        needs_rebuild && break
    end
    needs_rebuild || return args
    out = QTermDict()
    sizehint!(out, length(args))
    for (term, c) in args
        kept = NonEqualPair[]
        for p in term.ne
            pair_referenced(p, term.ops, c, indices) && push_ne_unique!(kept, p)
        end
        kept_ne = isempty(kept) ? EMPTY_NE : insertion_sort!(kept, isless)
        addto_key!(out, QTerm(copy(term.ops), kept_ne), c)
    end
    return out
end

@inline function pair_referenced(
        p::NonEqualPair, ops::Vector{Op}, c::CNum, scope::Vector{Index},
    )
    α, β = p
    α in scope && return true
    β in scope && return true
    depends_on_index_term(c, ops, α) && return true
    depends_on_index_term(c, ops, β) && return true
    return false
end

"""
    constraint_pairs(q::QAdd) -> Vector{Tuple{Index, Index}}

Return the deduplicated union of every term's `non_equal` pairs in `q`. This is
an introspection helper only; it does not define the expression semantics.

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h); j = Index(h, :j, 3, h);

julia> q = Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i, [j]);

julia> length(SecondQuantizedAlgebra.constraint_pairs(q))
1
```
"""
function constraint_pairs(q::QAdd)
    return constraint_pairs(q.arguments)
end

function single_qadd(c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair} = EMPTY_NE)
    iszero_cnum(c) && return QAdd(QTermDict(), EMPTY_INDICES)
    d = QTermDict()
    addto!(d, ops, c, ne)
    return QAdd(d, EMPTY_INDICES)
end

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

"""
    Base.one(q::QField) -> QAdd

Multiplicative identity as a unit `QAdd`.
"""
Base.one(::Type{<:QField}) = single_qadd(CNUM_ONE, Op[])
Base.one(q::QField) = one(typeof(q))

"""
    Base.isone(q::QField) -> Bool

True iff `q` is structurally the multiplicative identity.
"""
Base.isone(::QSym) = false
function Base.isone(q::QAdd)
    length(q.arguments) == 1 || return false
    (term, c) = first(q.arguments)
    isempty(term.ops) || return false
    re, im = realimag(c)
    iszero_num(im) || return false
    v = SymbolicUtils.unwrap(re)
    return (v isa Number && isone(v)) || isequal(re, NUM_ONE)
end

# `indices` is a multiset of bound sum indices (`Σ_iΣ_j ≡ Σ_jΣ_i`); compare/hash it
# order-insensitively but multiplicity-faithfully. A plain subset test is asymmetric
# and the XOR hash collapses repeated indices, so both use per-element multiplicity.
function same_index_set(a::Vector{Index}, b::Vector{Index})
    length(a) == length(b) || return false
    for idx in a
        count(x -> isequal(x, idx), a) == count(x -> isequal(x, idx), b) || return false
    end
    return true
end

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    return same_index_set(a.indices, b.indices)
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
function Base.hash(q::QAdd, h::UInt)
    # `+` (not `⊻`) so repeated indices contribute, matching `same_index_set`.
    ih = foldl((acc, idx) -> acc + hash(idx), q.indices; init = zero(UInt))
    return hash(:QAdd, hash(q.arguments, hash(ih, hash(length(q.indices), h))))
end

function Base.adjoint(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        rev = Op[adjoint(o) for o in Iterators.reverse(t.ops)]
        canonicalize!(out, rev, conj_cnum(c), t.ne)
    end
    return QAdd(out, copy(q.indices))
end

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

ne_sort_key(ne::Vector{NonEqualPair}) = Tuple(ne)

# Type-stable display order. Compares operator sequences field-by-field (each
# `full_op_key` is a fixed-length concrete tuple), then the constraint set. This
# replaces `isless` on the variable-length `Vararg` tuples the `_*sort_key`
# helpers build, which inference cannot resolve concretely (a print-path cost).
function ne_less(a::Vector{NonEqualPair}, b::Vector{NonEqualPair})
    la, lb = length(a), length(b)
    @inbounds for i in 1:min(la, lb)
        a[i] != b[i] && return isless(a[i], b[i])
    end
    return la < lb
end

function sorted_args_less(x::Pair{QTerm, CNum}, y::Pair{QTerm, CNum})
    xo, yo = x.first.ops, y.first.ops
    lx, ly = length(xo), length(yo)
    lx != ly && return lx < ly
    @inbounds for i in 1:lx
        kx, ky = full_op_key(xo[i]), full_op_key(yo[i])
        kx != ky && return isless(kx, ky)
    end
    return ne_less(x.first.ne, y.first.ne)
end

"""
    sorted_arguments(q::QAdd) -> Vector{QAdd}

Return each term of `q` as a single-entry [`QAdd`](@ref), in deterministic sort
order.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> length(SecondQuantizedAlgebra.sorted_arguments(a + a'))
2
```
"""
function sorted_arguments(q::QAdd)
    isempty(q.arguments) && return QAdd[]
    pairs = collect(q.arguments)
    insertion_sort!(pairs, sorted_args_less)
    return QAdd[single_qadd(c, term.ops, term.ne) for (term, c) in pairs]
end

# Total: without the packed level/axis fields two axes of one spin triple (or two levels of
# one transition set) tie, and the display order falls back to dict iteration order.
full_op_key(op::QSym) = (
    sort_key(op)..., type_order(op), name_rank(op.name_id),
    Int(op.l1), Int(op.l2), Int(op.g), Int(op.nlev),
)

"""
    term_order_key(t::QTerm) -> Tuple

Total ordering key for an operator product, built from `order_key`: length, then
operator-by-operator, then the canonical-sorted non-equal constraints.
"""
term_order_key(t::QTerm) = (length(t.ops), map(order_key, t.ops), ne_sort_key(t.ne))

"""
    qadd_order_key(q::QAdd) -> Vector

Total, reproducible ordering key for a sum: its term/coefficient pairs in sorted order, so
two `QAdd`s compare with `<` on their keys and tie exactly when they are `isequal`. The
coefficient contributes its printed form (a reproducible tiebreak, not a numeric order).
"""
function qadd_order_key(q::QAdd)
    pairs = [(term_order_key(t), coeff_key(c)) for (t, c) in q.arguments]
    insertion_sort!(pairs, isless)
    return pairs
end

Base.isless(a::QAdd, b::QAdd) = isless(qadd_order_key(a), qadd_order_key(b))
Base.isless(a::Type{<:QField}, b::Type{<:QField}) = isless(nameof(a), nameof(b))

function coeff_key(c::CNum)
    re, im = realimag(c)
    return (string(re), string(im))
end

"""
    Base.getindex(q::QAdd, key::AbstractVector{<:QSym}) -> Complex{Num}

Look up the prefactor for a given operator sequence. Returns zero if absent.
Throws if more than one constrained term shares the same operator sequence.
"""
function Base.getindex(q::QAdd, key::AbstractVector{<:QSym})
    found = nothing
    for (term, c) in q.arguments
        term.ops == key || continue
        found === nothing || throw(
            ArgumentError(
                "operator sequence has multiple constrained terms; iterate `q.arguments` " *
                    "or `sorted_arguments(q)` to inspect them explicitly"
            )
        )
        found = c
    end
    found === nothing && return to_num(CNUM_ZERO)
    return to_num(found)
end

"""
    Base.haskey(q::QAdd, key::AbstractVector{<:QSym}) -> Bool

Return `true` if some stored term has operator sequence `key`, ignoring
constraint scope. Pair with [`Base.getindex`](@ref) for the unique-prefactor
lookup; iterate `q.arguments` for full scope-aware access.
"""
function Base.haskey(q::QAdd, key::AbstractVector{<:QSym})
    for term in keys(q.arguments)
        term.ops == key && return true
    end
    return false
end

"""
    get_prefactor(s::QAdd) -> Complex{Num}

Return the `Complex{Num}` prefactor of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term
expressions, iterate over the `QAdd` directly.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> get_prefactor(2 * a' * a)
2
```

See also [`get_operators`](@ref), [`sorted_arguments`](@ref).
"""
function get_prefactor(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("get_prefactor requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return to_num(first(values(s.arguments)))
end

"""
    get_operators(s::QAdd) -> Vector{Op}

Return the unique operators occurring in [`QAdd`](@ref) `s`.

The result is sorted by the canonical [`order_key`](@ref). An operator and
its adjoint are distinct entries; use [`unique_up_to_adjoint`](@ref) when adjoint pairs
should be treated as one degree of freedom. Repeated factors are returned once;
read `term.ops` directly when multiplicity must be preserved.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> length(get_operators(a' * a))
2
```

See also [`get_prefactor`](@ref), [`unique_up_to_adjoint`](@ref), [`sorted_arguments`](@ref).
"""
function get_operators(s::QAdd)
    ops = Op[]
    for term in keys(s.arguments)
        append!(ops, term.ops)
    end
    unique!(ops)
    sort!(ops; by = order_key)
    return ops
end

"""
    get_variables(s::QAdd) -> Vector{SymbolicUtils.BasicSymbolic}

Return the unique scalar symbolic variables occurring in the prefactors of
[`QAdd`](@ref) `s`. Operator leaves are not included. Variables in indexed
coefficients are returned according to `Symbolics.get_variables`; use
[`get_indices`](@ref) to inspect the expression's structured index scope.

The returned variables follow the representation used by
`Symbolics.get_variables`, so they can be passed directly to Symbolics
functions such as `build_function`. Pass a `varlist` as the second argument
to restrict the result, or use `Symbolics.get_variables!` to fill a caller-owned
buffer.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables ω g;

julia> length(get_variables(ω * a' * a + g * a))
2
```
"""
function Symbolics.get_variables(s::QAdd; kwargs...)::Vector{SymbolicUtils.BasicSymbolic}
    vars = SymbolicUtils.BasicSymbolic[]
    for c in values(s.arguments)
        append!(vars, Symbolics.get_variables(real(c); kwargs...))
        append!(vars, Symbolics.get_variables(imag(c); kwargs...))
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables(
        s::QAdd,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )::Vector{SymbolicUtils.BasicSymbolic}
    vars = SymbolicUtils.BasicSymbolic[]
    for c in values(s.arguments)
        append!(
            vars,
            Symbolics.get_variables(
                real(c), varlist; is_atomic = is_atomic, kwargs...
            ),
        )
        append!(
            vars,
            Symbolics.get_variables(
                imag(c), varlist; is_atomic = is_atomic, kwargs...
            ),
        )
    end
    unique!(vars)
    sort!(vars; by = string)
    return vars
end

function Symbolics.get_variables!(buffer, s::QAdd; kwargs...)
    for c in values(s.arguments)
        Symbolics.get_variables!(buffer, real(c); kwargs...)
        Symbolics.get_variables!(buffer, imag(c); kwargs...)
    end
    return buffer
end

function Symbolics.get_variables!(
        buffer,
        s::QAdd,
        varlist;
        is_atomic = SymbolicUtils.default_is_atomic,
        kwargs...,
    )
    for c in values(s.arguments)
        Symbolics.get_variables!(
            buffer, real(c), varlist; is_atomic = is_atomic, kwargs...
        )
        Symbolics.get_variables!(
            buffer, imag(c), varlist; is_atomic = is_atomic, kwargs...
        )
    end
    return buffer
end
