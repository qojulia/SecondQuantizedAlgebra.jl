# A pair of indices `(α, β)` meaning `α ≠ β`.
const NonEqualPair = Tuple{Index, Index}

# Shared sentinel for "no constraints". Never mutated.
const _EMPTY_NE = NonEqualPair[]

"""
    QTermKey

Internal identity of a term in a `QAdd`: an ordered operator product plus the
pairwise index inequality constraints that scope that product.
"""
struct QTermKey
    ops::Vector{QSym}
    ne::Vector{NonEqualPair}
end

Base.isequal(a::QTermKey, b::QTermKey) = isequal(a.ops, b.ops) && isequal(a.ne, b.ne)
Base.:(==)(a::QTermKey, b::QTermKey) = isequal(a, b)
function Base.hash(term::QTermKey, h::UInt)
    return hash(:QTermKey, hash(term.ops, hash(term.ne, h)))
end

"""
    QTermDict

Alias for `Dict{QTermKey, CNum}` — the storage type backing every [`QAdd`](@ref).
Each entry maps an exact term identity `(ops, ne)` to its prefactor. Iteration
yields `Pair{QTermKey, CNum}`; callers reach `term.ops` and `term.ne` directly.
"""
const QTermDict = Dict{QTermKey, CNum}

@inline _copy_ne(ne::Vector{NonEqualPair}) = isempty(ne) ? _EMPTY_NE : copy(ne)

@inline _canonical_pair(p::NonEqualPair) = isless(p[2], p[1]) ? (p[2], p[1]) : p

function _push_ne_unique!(out::Vector{NonEqualPair}, p::NonEqualPair)
    cp = _canonical_pair(p)
    cp ∉ out && push!(out, cp)
    return out
end

"""
    _canonical_ne(ne) -> Vector{NonEqualPair}

Canonical form for a constraint set: each pair oriented via [`_canonical_pair`](@ref),
duplicates removed, then sorted. The output is the form stored in every
[`QTermKey`](@ref)`.ne` — two constraint sets compare equal as dict keys iff their
canonical forms are identical. Returns the shared `_EMPTY_NE` sentinel for empty input.
"""
function _canonical_ne(ne::Vector{NonEqualPair})
    isempty(ne) && return _EMPTY_NE
    out = NonEqualPair[]
    sizehint!(out, length(ne))
    for p in ne
        _push_ne_unique!(out, p)
    end
    sort!(out)
    return isempty(out) ? _EMPTY_NE : out
end

@inline function _ne_contains(ne::Vector{NonEqualPair}, α::Index, β::Index)
    return _canonical_pair((α, β)) in ne
end

function _merge_ne(a::Vector{NonEqualPair}, b::Vector{NonEqualPair})
    isempty(b) && return a
    out = a === _EMPTY_NE ? NonEqualPair[] : copy(a)
    for p in b
        _push_ne_unique!(out, p)
    end
    return isempty(out) ? _EMPTY_NE : out
end

function _merge_ne_pair(ne::Vector{NonEqualPair}, α::Index, β::Index)
    return _merge_ne(ne, NonEqualPair[(α, β)])
end

function _copy_key(term::QTermKey)
    return QTermKey(copy(term.ops), _copy_ne(term.ne))
end

"""
    _term_key(ops, ne = _EMPTY_NE) -> QTermKey

Storage-key constructor: copies `ops` (callers can mutate freely afterwards) and
canonicalizes `ne`. Every [`QTermDict`](@ref) insertion goes through this so that
structural equality of `(ops, ne)` always implies equal hash keys.
"""
function _term_key(ops::Vector{QSym}, ne::Vector{NonEqualPair} = _EMPTY_NE)
    return QTermKey(copy(ops), _canonical_ne(ne))
end

function _push_key_unique!(out::Vector{QTermKey}, term::QTermKey)
    term in out || push!(out, term)
    return out
end

"""
    _addto!(d, ops, c[, ne]) -> d

Insert (or merge into) the constrained term `(ops, ne)` in `d`. Like-term
collection applies only when both the operator sequence and the constraint set
match exactly.
"""
function _addto!(
        d::QTermDict, ops::Vector{QSym}, c::CNum,
        ne::Vector{NonEqualPair} = _EMPTY_NE
    )
    return _addto_key!(d, _term_key(ops, ne), c)
end

"""
    _addto_key!(d, term, c) -> d

Like-term collection at the dict level. If `term` is absent, store `c` (unless
zero). If present, add into the existing prefactor and delete the entry when
the result is zero. Zero coefficients are never stored, so the dict never needs
a separate cleanup pass.
"""
function _addto_key!(d::QTermDict, term::QTermKey, c::CNum)
    existing = get(d, term, nothing)
    if existing === nothing
        _iszero_cnum(c) && return d
        d[term] = c
        return d
    end
    new_c = _add_cnum(existing, c)
    if _iszero_cnum(new_c)
        delete!(d, term)
    else
        d[term] = new_c
    end
    return d
end

function _copy_args(d::QTermDict)
    out = QTermDict()
    sizehint!(out, length(d))
    for (term, c) in d
        out[_copy_key(term)] = c
    end
    return out
end

function _constraint_pairs(args::QTermDict)
    out = NonEqualPair[]
    for term in keys(args), p in term.ne
        _push_ne_unique!(out, p)
    end
    return out
end

function _merge_unique(a::Vector{T}, b::Vector{T}) where {T}
    isempty(a) && return b
    isempty(b) && return a
    result = copy(a)
    for x in b
        x ∉ result && push!(result, x)
    end
    return result
end
