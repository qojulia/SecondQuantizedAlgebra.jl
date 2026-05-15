# A pair of indices `(α, β)` meaning `α ≠ β`.
const NonEqualPair = Tuple{Index, Index}

# Legacy internal sentinels. New code should prefer `_empty_ne()` / `_empty_bound()`
# so an accidental mutation of one empty vector cannot leak into later terms.
const _EMPTY_NE = NonEqualPair[]
const _EMPTY_BOUND = Index[]
_empty_ne() = NonEqualPair[]
_empty_bound() = Index[]

struct _TermScope
    bound::Vector{Index}
    ne::Vector{NonEqualPair}
end

struct _SiteFactor
    index::Index
    ops::Vector{QSym}
end

struct _SpaceChain
    space_index::Int
    factors::Vector{_SiteFactor}
end

struct _SiteBlock
    space_index::Int
    index::Index
    ops::Vector{QSym}
end

struct _Monomial
    scope::_TermScope
    chains::Vector{_SpaceChain}
    ops::Vector{QSym}
end

"""
    QTerm

A single canonical monomial of a [`QAdd`](@ref): explicit scope plus canonical
per-space factor chains.
"""
struct QTerm
    monomial::_Monomial
end

QTerm(ops::Vector{QSym}, ne::Vector{NonEqualPair}) = _term_key(ops, _empty_bound(), ne, ops)

@inline _copy_ne(ne::Vector{NonEqualPair}) = isempty(ne) ? _empty_ne() : copy(ne)
@inline _copy_bound(bound::Vector{Index}) = isempty(bound) ? _empty_bound() : copy(bound)

@inline _canonical_pair(p::NonEqualPair) = isless(p[2], p[1]) ? (p[2], p[1]) : p

function _push_ne_unique!(out::Vector{NonEqualPair}, p::NonEqualPair)
    cp = _canonical_pair(p)
    cp ∉ out && push!(out, cp)
    return out
end

function _push_index_unique!(out::Vector{Index}, idx::Index)
    idx ∉ out && push!(out, idx)
    return out
end

function _canonical_ne(ne::Vector{NonEqualPair})
    isempty(ne) && return _empty_ne()
    out = NonEqualPair[]
    sizehint!(out, length(ne))
    for p in ne
        _push_ne_unique!(out, p)
    end
    sort!(out)
    return isempty(out) ? _empty_ne() : out
end

function _canonical_bound(bound::Vector{Index})
    isempty(bound) && return _empty_bound()
    out = Index[]
    sizehint!(out, length(bound))
    for idx in bound
        _push_index_unique!(out, idx)
    end
    return isempty(out) ? _empty_bound() : out
end

@inline function _ne_contains(ne::Vector{NonEqualPair}, α::Index, β::Index)
    return _canonical_pair((α, β)) in ne
end

function _merge_ne(a::Vector{NonEqualPair}, b::Vector{NonEqualPair})
    isempty(b) && return a
    out = copy(a)
    for p in b
        _push_ne_unique!(out, p)
    end
    return isempty(out) ? _empty_ne() : out
end

function _merge_ne_pair(ne::Vector{NonEqualPair}, α::Index, β::Index)
    return _merge_ne(ne, NonEqualPair[(α, β)])
end

function _merge_bound(a::Vector{Index}, b::Vector{Index})
    isempty(b) && return a
    out = copy(a)
    for idx in b
        _push_index_unique!(out, idx)
    end
    return isempty(out) ? _empty_bound() : out
end

function _merge_bound_idx(bound::Vector{Index}, idx::Index)
    return _merge_bound(bound, Index[idx])
end

function _drop_bound_with(bound::Vector{Index}, idx::Index)
    isempty(bound) && return _empty_bound()
    out = Index[]
    for x in bound
        x == idx && continue
        push!(out, x)
    end
    return isempty(out) ? _empty_bound() : out
end

function _substitute_ne(ne::Vector{NonEqualPair}, from::Index, to::Index)
    isempty(ne) && return _empty_ne()
    out = NonEqualPair[]
    for (α, β) in ne
        new_α = α == from ? to : α
        new_β = β == from ? to : β
        new_α == new_β && continue
        _push_ne_unique!(out, (new_α, new_β))
    end
    return isempty(out) ? _empty_ne() : out
end

function _substitute_bound(bound::Vector{Index}, from::Index, to::Index)
    isempty(bound) && return _empty_bound()
    out = Index[]
    for idx in bound
        _push_index_unique!(out, idx == from ? to : idx)
    end
    return isempty(out) ? _empty_bound() : out
end

function _drop_ne_with(ne::Vector{NonEqualPair}, idx::Index)
    isempty(ne) && return _empty_ne()
    out = NonEqualPair[]
    for (α, β) in ne
        α == idx && continue
        β == idx && continue
        _push_ne_unique!(out, (α, β))
    end
    return isempty(out) ? _empty_ne() : out
end

function _site_order_key(idx::Index)
    return has_index(idx) ? (1, idx.name) : (0, Symbol(""))
end

function _copy_factor(f::_SiteFactor)
    return _SiteFactor(f.index, copy(f.ops))
end

function _copy_chain(chain::_SpaceChain)
    return _SpaceChain(chain.space_index, [_copy_factor(f) for f in chain.factors])
end

function _copy_chains(chains::Vector{_SpaceChain})
    return [_copy_chain(chain) for chain in chains]
end

function _flatten_chain(chain::_SpaceChain)
    out = QSym[]
    sizehint!(out, sum(length(f.ops) for f in chain.factors))
    for factor in chain.factors
        append!(out, factor.ops)
    end
    return out
end

function _flatten_chains(chains::Vector{_SpaceChain})
    isempty(chains) && return QSym[]
    out = QSym[]
    sizehint!(out, sum(sum(length(f.ops) for f in chain.factors) for chain in chains))
    for chain in chains
        append!(out, _flatten_chain(chain))
    end
    return out
end

function _canonical_ops_from_chains(chains::Vector{_SpaceChain})
    isempty(chains) && return QSym[]
    out = QSym[]
    for chain in chains
        factors = sort([_copy_factor(f) for f in chain.factors]; by = f -> _site_order_key(f.index), alg = Base.MergeSort)
        for factor in factors
            append!(out, factor.ops)
        end
    end
    return out
end

function _site_blocks(chains::Vector{_SpaceChain})
    isempty(chains) && return _SiteBlock[]
    blocks = _SiteBlock[]
    for chain in chains, factor in chain.factors
        push!(blocks, _SiteBlock(chain.space_index, factor.index, copy(factor.ops)))
    end
    return blocks
end

_site_blocks(ops::Vector{QSym}) = _site_blocks(_chains_from_ops(ops))

function _chains_from_ops(ops::Vector{QSym})
    isempty(ops) && return _SpaceChain[]
    by_space = Dict{Int, Vector{_SiteFactor}}()
    seen_spaces = Int[]
    for op in ops
        si = op.space_index
        if !haskey(by_space, si)
            by_space[si] = _SiteFactor[]
            push!(seen_spaces, si)
        end
        factors = by_space[si]
        if !isempty(factors) && factors[end].index == op.index
            push!(factors[end].ops, op)
        else
            push!(factors, _SiteFactor(op.index, QSym[op]))
        end
    end
    sort!(seen_spaces)
    chains = Vector{_SpaceChain}(undef, length(seen_spaces))
    for (i, si) in enumerate(seen_spaces)
        chains[i] = _SpaceChain(si, by_space[si])
    end
    return chains
end

function _term_scope(
        bound::Vector{Index} = _empty_bound(),
        ne::Vector{NonEqualPair} = _empty_ne()
    )
    return _TermScope(_canonical_bound(bound), _canonical_ne(ne))
end

function _monomial(
        ops::Vector{QSym},
        bound::Vector{Index} = _empty_bound(),
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    chains = _chains_from_ops(copy(ops))
    flat_ops = _canonical_ops_from_chains(chains)
    return _Monomial(_term_scope(bound, ne), chains, flat_ops)
end

function _monomial_from_chains(
        chains::Vector{_SpaceChain},
        bound::Vector{Index} = _empty_bound(),
        ne::Vector{NonEqualPair} = _empty_ne()
    )
    chains_copy = _copy_chains(chains)
    return _Monomial(_term_scope(bound, ne), chains_copy, _canonical_ops_from_chains(chains_copy))
end

function Base.getproperty(term::QTerm, s::Symbol)
    if s === :ops
        return getfield(term, :monomial).ops
    elseif s === :ne
        return getfield(term, :monomial).scope.ne
    elseif s === :phys_ops
        return getfield(term, :monomial).ops
    elseif s === :bound
        return getfield(term, :monomial).scope.bound
    elseif s === :blocks
        return _site_blocks(getfield(term, :monomial).chains)
    elseif s === :chains
        return getfield(term, :monomial).chains
    end
    return getfield(term, s)
end

function Base.propertynames(::QTerm, private::Bool = false)
    props = (:monomial, :ops, :ne, :phys_ops, :bound, :blocks, :chains)
    return private ? props : Base.tail(props)
end

function _alpha_placeholder(idx::Index, slot::Int)
    name = Symbol("_b", slot)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, idx.range, idx.space_index, Num(sym_var))
end

function _alpha_bound_pairs(bound::Vector{Index})
    isempty(bound) && return Pair{Index, Index}[]
    pairs = Pair{Index, Index}[]
    sizehint!(pairs, length(bound))
    for (slot, idx) in enumerate(bound)
        push!(pairs, idx => _alpha_placeholder(idx, slot))
    end
    return pairs
end

function _rename_ops(ops::Vector{QSym}, pairs::Vector{Pair{Index, Index}})
    isempty(pairs) && return ops
    out = copy(ops)
    for (from, to) in pairs
        out = QSym[change_index(op, from, to) for op in out]
    end
    return out
end

function _rename_bound(bound::Vector{Index}, pairs::Vector{Pair{Index, Index}})
    isempty(pairs) && return bound
    out = bound
    for (from, to) in pairs
        out = _substitute_bound(out, from, to)
    end
    return out
end

function _rename_ne(ne::Vector{NonEqualPair}, pairs::Vector{Pair{Index, Index}})
    isempty(pairs) && return ne
    out = ne
    for (from, to) in pairs
        out = _substitute_ne(out, from, to)
    end
    return out
end

function _rename_factor_signature(f::_SiteFactor, pairs::Vector{Pair{Index, Index}})
    idx = f.index
    ops = f.ops
    for (from, to) in pairs
        idx = idx == from ? to : idx
        ops = QSym[change_index(op, from, to) for op in ops]
    end
    return (idx, Tuple(ops))
end

function _rename_chain_signature(chains::Vector{_SpaceChain}, pairs::Vector{Pair{Index, Index}})
    isempty(chains) && return ()
    return Tuple(
        (chain.space_index, Tuple(_rename_factor_signature(f, pairs) for f in chain.factors))
        for chain in chains
    )
end

function _alpha_signature(term::QTerm)
    pairs = _alpha_bound_pairs(term.bound)
    bound = _rename_bound(term.bound, pairs)
    ne = _rename_ne(term.ne, pairs)
    chains = _rename_chain_signature(term.chains, pairs)
    return bound, ne, chains
end

function Base.isequal(a::QTerm, b::QTerm)
    a_bound, a_ne, a_chains = _alpha_signature(a)
    b_bound, b_ne, b_chains = _alpha_signature(b)
    return isequal(a_bound, b_bound) && isequal(a_ne, b_ne) && isequal(a_chains, b_chains)
end
Base.:(==)(a::QTerm, b::QTerm) = isequal(a, b)

function Base.hash(term::QTerm, h::UInt)
    bound, ne, chains = _alpha_signature(term)
    return hash(:QTerm, hash(bound, hash(ne, hash(chains, h))))
end

"""
    QTermDict

Alias for `Dict{QTerm, CNum}` — the canonical sum storage backing every
[`QAdd`](@ref).
"""
const QTermDict = Dict{QTerm, CNum}

function _copy_key(term::QTerm)
    return QTerm(_monomial_from_chains(term.chains, term.bound, term.ne))
end

function _term_key(
        ops::Vector{QSym},
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    return QTerm(_monomial(ops, _empty_bound(), ne, trace_ops))
end

function _term_key(
        ops::Vector{QSym},
        bound::Vector{Index},
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    return QTerm(_monomial(ops, bound, ne, trace_ops))
end

function _term_key(
        chains::Vector{_SpaceChain},
        bound::Vector{Index},
        ne::Vector{NonEqualPair}
    )
    return QTerm(_monomial_from_chains(chains, bound, ne))
end

function _push_key_unique!(out::Vector{QTerm}, term::QTerm)
    term in out || push!(out, term)
    return out
end

function _used_indices(c::CNum, ops::Vector{QSym}, candidates::Vector{Index})
    isempty(candidates) && return _empty_bound()
    used = Index[]
    for idx in candidates
        _depends_on_index_term(c, ops, idx) && _push_index_unique!(used, idx)
    end
    return isempty(used) ? _empty_bound() : used
end

function _normalize_scope(
        c::CNum, ops::Vector{QSym},
        bound::Vector{Index}, ne::Vector{NonEqualPair}
    )
    candidates = copy(bound)
    for op in ops
        has_index(op.index) && _push_index_unique!(candidates, op.index)
    end
    used = _used_indices(c, ops, candidates)
    for (α, β) in ne
        _push_index_unique!(used, α)
        _push_index_unique!(used, β)
    end
    kept_bound = _empty_bound()
    if !isempty(bound)
        tmp = Index[]
        for idx in bound
            idx in used && push!(tmp, idx)
        end
        kept_bound = isempty(tmp) ? _empty_bound() : tmp
    end
    kept_ne = _canonical_ne(ne)
    return kept_bound, kept_ne
end

function _rename_coeff_bound(c::CNum, from_bound::Vector{Index}, to_bound::Vector{Index})
    from_bound == to_bound && return c
    isempty(from_bound) && return c
    temps = Pair{Index, Index}[]
    sizehint!(temps, length(from_bound))
    renamed = c
    for (slot, idx) in enumerate(from_bound)
        tmp = _alpha_placeholder(idx, length(from_bound) + slot)
        push!(temps, idx => tmp)
        renamed = change_index(renamed, idx, tmp)
    end
    for ((_, tmp), idx) in zip(temps, to_bound)
        renamed = change_index(renamed, tmp, idx)
    end
    return renamed
end

function _addto!(
        d::QTermDict, ops::Vector{QSym}, c::CNum,
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    bound, normalized_ne = _normalize_scope(c, ops, _empty_bound(), ne)
    return _addto_key!(d, _term_key(ops, bound, normalized_ne, trace_ops), c)
end

function _addto!(
        d::QTermDict, ops::Vector{QSym}, c::CNum,
        bound::Vector{Index},
        ne::Vector{NonEqualPair} = _empty_ne(),
        trace_ops::Vector{QSym} = ops
    )
    normalized_bound, normalized_ne = _normalize_scope(c, ops, bound, ne)
    return _addto_key!(d, _term_key(ops, normalized_bound, normalized_ne, trace_ops), c)
end

function _addto_key!(d::QTermDict, term::QTerm, c::CNum)
    existing_key = nothing
    for key in keys(d)
        if isequal(key, term)
            existing_key = key
            break
        end
    end
    if existing_key === nothing
        _iszero_cnum(c) && return d
        d[term] = c
        return d
    end
    c = _rename_coeff_bound(c, term.bound, existing_key.bound)
    existing = d[existing_key]
    new_c = _add_cnum(existing, c)
    if _iszero_cnum(new_c)
        delete!(d, existing_key)
    else
        d[existing_key] = new_c
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

function _bound_indices(args::QTermDict)
    out = Index[]
    for term in keys(args), idx in term.bound
        _push_index_unique!(out, idx)
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

_uses_trace_key(::QTerm) = false
