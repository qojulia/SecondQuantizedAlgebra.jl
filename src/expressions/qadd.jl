"""
    QAdd <: QField

The sole compound expression type — a sum of eagerly-ordered operator products.

All arithmetic on [`QSym`](@ref) operators returns `QAdd`. Internally stores a
dictionary whose keys are exact term identities `(ops, non_equal)` and whose
values are `Complex{Num}` prefactors. Like-term collection applies only to
terms with identical operator strings *and* identical scoped constraints.

# Fields
- `arguments::QTermDict` — `Dict{QTerm, CNum}` of exact term identity → prefactor
- `indices::Vector{Index}` — summation indices (empty for a regular sum)

# Iteration
Iterating over a `QAdd` yields `Pair{QTerm, CNum}` entries; access
`term.ops` and `term.ne` on the key directly.

See also [`QTerm`](@ref), [`prefactor`](@ref), [`operators`](@ref), [`Σ`](@ref),
[`constraint_pairs`](@ref).
"""
struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
end

"""
    constraint_pairs(q::QAdd) -> Vector{Tuple{Index, Index}}

Return the deduplicated union of every term's `non_equal` pairs in `q`. This is
an introspection helper only; it does not define the expression semantics.
"""
function constraint_pairs(q::QAdd)
    return _constraint_pairs(q.arguments)
end

# Convenience: single-term QAdd
function _single_qadd(c::CNum, ops::Vector{QSym}, ne::Vector{NonEqualPair} = _EMPTY_NE)
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[])
    d = QTermDict()
    _addto!(d, ops, c, ne)
    return QAdd(d, Index[])
end

function _qadd_from_oterms(
        terms::Vector{OrderedTerm}, indices::Vector{Index},
        ne::Vector{NonEqualPair} = _EMPTY_NE
    )
    d = QTermDict()
    for t in terms
        for ot in _apply_ordering(t.prefactor, t.ops, get_ordering())
            _addto!(d, ot.ops, ot.prefactor, ne)
        end
    end
    return QAdd(d, indices)
end

function _ordered_qadd(c::CNum, ops::Vector{QSym}, ne::Vector{NonEqualPair} = _EMPTY_NE)
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[])
    d = QTermDict()
    for t in _apply_ordering(c, ops, get_ordering())
        _addto!(d, t.ops, t.prefactor, ne)
    end
    return QAdd(d, Index[])
end

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    a.indices == b.indices || return false
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
function Base.hash(q::QAdd, h::UInt)
    return hash(:QAdd, hash(q.arguments, hash(q.indices, h)))
end

function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (term, c) in q.arguments
        adj_ops = QSym[adjoint(op) for op in reverse(term.ops)]
        _site_sort!(adj_ops)
        _addto!(d, adj_ops, conj(c), term.ne)
    end
    return QAdd(d, copy(q.indices))
end

# --- Iteration: yield `term::QTerm => coeff::CNum` directly ---

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

# --- Sorted term access for printing and comparison ---

_ne_sort_key(ne::Vector{NonEqualPair}) = Tuple(ne)

"""
    sorted_arguments(q::QAdd) -> Vector{QAdd}

Return each term of `q` as a single-entry [`QAdd`](@ref), in deterministic sort
order.
"""
function sorted_arguments(q::QAdd)
    isempty(q.arguments) && return QAdd[]
    pairs = sort!(
        collect(q.arguments);
        by = p -> (_term_sort_key(p.first.ops), _ne_sort_key(p.first.ne)),
    )
    return QAdd[_single_qadd(c, term.ops, term.ne) for (term, c) in pairs]
end

_term_sort_key(ops::Vector{QSym}) = (length(ops), map(_full_op_key, ops)...)
_full_op_key(op::QSym) = (_sort_key(op)..., _type_order(op), op.name)
_type_order(::Destroy) = 0
_type_order(::Create) = 1
_type_order(::Transition) = 2
_type_order(::Pauli) = 3
_type_order(::Spin) = 4
_type_order(::Position) = 5
_type_order(::Momentum) = 6

"""
    Base.getindex(q::QAdd, key::AbstractVector{<:QSym}) -> CNum

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
    found === nothing && return _CNUM_ZERO
    return found
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

# --- QAdd accessor helpers ---

"""
    prefactor(s::QAdd) -> CNum

Return the `Complex{Num}` prefactor of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term
expressions, iterate over the `QAdd` directly.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
prefactor(2 * a' * a)   # 2 + 0im
```

See also [`operators`](@ref), [`sorted_arguments`](@ref).
"""
function prefactor(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("prefactor requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return first(values(s.arguments))
end

"""
    operators(s::QAdd) -> Vector{QSym}

Return the ordered operator sequence of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term
expressions, iterate over the `QAdd` directly and read `term.ops` from each
[`QTerm`](@ref).

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
operators(a' * a)   # [Create(:a, 1, NO_INDEX, ...), Destroy(:a, 1, NO_INDEX)]
```

See also [`prefactor`](@ref), [`sorted_arguments`](@ref).
"""
function operators(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("operators requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return first(keys(s.arguments)).ops
end
