"""
    QAdd <: QField

The sole compound expression type — a sum of canonical monomials.

All arithmetic on [`QSym`](@ref) operators returns `QAdd`. Internally stores a
dictionary whose keys are exact term identities and whose values are `Complex{Num}`
prefactors.
"""
struct QAdd <: QField
    arguments::QTermDict
    indices::Vector{Index}
end

function _indices_union(d::QTermDict, extra::Vector{Index} = _empty_bound())
    out = _bound_indices(d)
    isempty(extra) && return out
    return _merge_unique(out, _canonical_bound(extra))
end

_qadd(d::QTermDict, indices::Vector{Index} = _empty_bound()) = QAdd(d, _indices_union(d, indices))

"""
    constraint_pairs(q::QAdd) -> Vector{Tuple{Index, Index}}

Return the deduplicated union of every term's `non_equal` pairs in `q`.
"""
function constraint_pairs(q::QAdd)
    return _constraint_pairs(q.arguments)
end

function _single_qadd(
        c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _empty_ne()
    )
    _iszero_cnum(c) && return _zero_qadd()
    d = QTermDict()
    _addto!(d, ops, c, ne)
    return _qadd(d)
end

function _single_qadd(
        c::CNum, ops::Vector{QSym},
        bound::Vector{Index},
        ne::Vector{NonEqualPair} = _empty_ne()
    )
    _iszero_cnum(c) && return _zero_qadd()
    d = QTermDict()
    _addto!(d, ops, c, bound, ne)
    return _qadd(d, bound)
end

function _qadd_from_oterms(
        terms::Vector{OrderedTerm}, indices::Vector{Index},
        ne::Vector{NonEqualPair} = _empty_ne()
    )
    d = QTermDict()
    for t in terms
        for ot in _apply_ordering(t.prefactor, t.ops, get_ordering())
            _addto!(d, ot.ops, ot.prefactor, indices, ne)
        end
    end
    return _qadd(d, indices)
end

function _ordered_qadd(
        c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    return _ordered_qadd(c, ops, _empty_bound(), ne, trace_ops)
end

function _ordered_qadd(
        c::CNum, ops::Vector{QSym},
        bound::Vector{Index},
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    _iszero_cnum(c) && return _zero_qadd()
    d = QTermDict()
    _accumulate_normalized!(d, c, ops, bound, ne, get_ordering())
    return _qadd(d, bound)
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
        adj_ops = QSym[adjoint(op) for op in reverse(_flatten_chains(term.chains))]
        _accumulate_normalized!(d, conj(c), adj_ops, term.bound, term.ne, get_ordering())
    end
    return _qadd(d, copy(q.indices))
end

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{QTerm, CNum}

_ne_sort_key(ne::Vector{NonEqualPair}) = Tuple(ne)
_bound_sort_key(bound::Vector{Index}) = Tuple(bound)

"""
    sorted_arguments(q::QAdd) -> Vector{QAdd}

Return each term of `q` as a single-entry [`QAdd`](@ref), in deterministic sort
order.
"""
function sorted_arguments(q::QAdd)
    isempty(q.arguments) && return QAdd[]
    pairs = sort!(
        collect(q.arguments);
        by = p -> (_term_sort_key(p.first.ops), _bound_sort_key(p.first.bound), _ne_sort_key(p.first.ne)),
    )
    return QAdd[_qadd(QTermDict(_copy_key(term) => c), copy(term.bound)) for (term, c) in pairs]
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
Throws if more than one scoped term shares the same operator sequence.
"""
function Base.getindex(q::QAdd, key::AbstractVector{<:QSym})
    found = nothing
    for (term, c) in q.arguments
        term.ops == key || continue
        found === nothing || throw(
            ArgumentError(
                "operator sequence has multiple scoped terms; iterate `q.arguments` " *
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

Return `true` if some stored term has operator sequence `key`, ignoring scope.
"""
function Base.haskey(q::QAdd, key::AbstractVector{<:QSym})
    for term in keys(q.arguments)
        term.ops == key && return true
    end
    return false
end

"""
    prefactor(s::QAdd) -> CNum

Return the prefactor of a single-term [`QAdd`](@ref).
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
"""
function operators(s::QAdd)
    length(s.arguments) == 1 || throw(
        ArgumentError("operators requires a single-term expression, got $(length(s.arguments)) terms")
    )
    return first(keys(s.arguments)).ops
end
