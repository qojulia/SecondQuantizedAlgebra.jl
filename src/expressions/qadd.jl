"""
    QAdd <: QField

The sole compound expression type: a sum of eagerly-ordered operator products.

All arithmetic on [`QSym`](@ref) operators produces a `QAdd`. Iterating over a
`QAdd` yields `Pair{QTerm, CNum}` entries; read `term.ops` for the operator
sequence and `term.ne` for the scoped index constraints.

See also [`QTerm`](@ref), [`prefactor`](@ref), [`operators`](@ref), [`Σ`](@ref),
[`constraint_pairs`](@ref).
"""
struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
    function QAdd(arguments::QTermDict, indices::Vector{Index})
        return new(_prune_dead_ne(arguments, indices), indices)
    end
end

# Strip NE pairs that reference no op index, coefficient index, or scope index.
function _prune_dead_ne(args::QTermDict, indices::Vector{Index})
    needs_rebuild = false
    for (term, c) in args
        for p in term.ne
            if !_pair_referenced(p, term.ops, c, indices)
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
            _pair_referenced(p, term.ops, c, indices) && _push_ne_unique!(kept, p)
        end
        kept_ne = isempty(kept) ? _EMPTY_NE : sort!(kept)
        _addto_key!(out, QTerm(copy(term.ops), kept_ne), c)
    end
    return out
end

@inline function _pair_referenced(
        p::NonEqualPair, ops::Vector{Op}, c::CNum, scope::Vector{Index},
    )
    α, β = p
    α in scope && return true
    β in scope && return true
    _depends_on_index_term(c, ops, α) && return true
    _depends_on_index_term(c, ops, β) && return true
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
    return _constraint_pairs(q.arguments)
end

function _single_qadd(c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair} = _EMPTY_NE)
    _iszero_cnum(c) && return QAdd(QTermDict(), _EMPTY_INDICES)
    d = QTermDict()
    _addto!(d, ops, c, ne)
    return QAdd(d, _EMPTY_INDICES)
end

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

"""
    Base.one(q::QField) -> QAdd

Multiplicative identity as a unit `QAdd`.
"""
Base.one(::Type{<:QField}) = _single_qadd(_CNUM_ONE, Op[])
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
    _iszero_num(imag(c)) || return false
    v = SymbolicUtils.unwrap(real(c))
    return (v isa Number && isone(v)) || isequal(real(c), _NUM_ONE)
end

# `indices` is a multiset of bound sum indices (`Σ_iΣ_j ≡ Σ_jΣ_i`); compare/hash it
# order-insensitively but multiplicity-faithfully. A plain subset test is asymmetric
# and the XOR hash collapses repeated indices, so both use per-element multiplicity.
function _same_index_set(a::Vector{Index}, b::Vector{Index})
    length(a) == length(b) || return false
    for idx in a
        count(x -> isequal(x, idx), a) == count(x -> isequal(x, idx), b) || return false
    end
    return true
end

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    return _same_index_set(a.indices, b.indices)
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
function Base.hash(q::QAdd, h::UInt)
    # `+` (not `⊻`) so repeated indices contribute, matching `_same_index_set`.
    ih = foldl((acc, idx) -> acc + hash(idx), q.indices; init = zero(UInt))
    return hash(:QAdd, hash(q.arguments, hash(ih, hash(length(q.indices), h))))
end

function Base.adjoint(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        rev = Op[adjoint(o) for o in Iterators.reverse(t.ops)]
        _canonicalize!(out, rev, _conj_cnum(c), t.ne)
    end
    return QAdd(out, copy(q.indices))
end

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

_ne_sort_key(ne::Vector{NonEqualPair}) = Tuple(ne)

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
    pairs = sort!(
        collect(q.arguments);
        by = p -> (_term_sort_key(p.first.ops), _ne_sort_key(p.first.ne)),
    )
    return QAdd[_single_qadd(c, term.ops, term.ne) for (term, c) in pairs]
end

_term_sort_key(ops::Vector{Op}) = (length(ops), map(_full_op_key, ops)...)
_full_op_key(op::QSym) = (_sort_key(op)..., _type_order(op), _name_rank(op.name_id))

"""
    term_order_key(t::QTerm) -> Tuple

Total ordering key for an operator product, built from `order_key`: length, then
operator-by-operator, then the canonical-sorted non-equal constraints.
"""
term_order_key(t::QTerm) = (length(t.ops), map(order_key, t.ops), _ne_sort_key(t.ne))

"""
    qadd_order_key(q::QAdd) -> Vector

Total, reproducible ordering key for a sum: its term/coefficient pairs in sorted order, so
two `QAdd`s compare with `<` on their keys and tie exactly when they are `isequal`. The
coefficient contributes its printed form (a reproducible tiebreak, not a numeric order).
"""
function qadd_order_key(q::QAdd)
    pairs = [(term_order_key(t), _coeff_key(c)) for (t, c) in q.arguments]
    sort!(pairs)
    return pairs
end

Base.isless(a::QAdd, b::QAdd) = isless(qadd_order_key(a), qadd_order_key(b))
Base.isless(a::Type{<:QField}, b::Type{<:QField}) = isless(nameof(a), nameof(b))

_coeff_key(c::CNum) = (string(real(c)), string(imag(c)))

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
    found === nothing && return to_num(_CNUM_ZERO)
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
    prefactor(s::QAdd) -> Complex{Num}

Return the `Complex{Num}` prefactor of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term
expressions, iterate over the `QAdd` directly.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> prefactor(2 * a' * a)
2
```

See also [`operators`](@ref), [`sorted_arguments`](@ref).
"""
function prefactor(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("prefactor requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return to_num(first(values(s.arguments)))
end

"""
    operators(s::QAdd) -> Vector{Op}

Return the ordered operator sequence of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term
expressions, iterate over the `QAdd` directly and read `term.ops` from each
[`QTerm`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> length(operators(a' * a))
2
```

See also [`prefactor`](@ref), [`sorted_arguments`](@ref).
"""
function operators(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("operators requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return first(keys(s.arguments)).ops
end
