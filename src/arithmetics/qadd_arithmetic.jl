# ============================================================================
#  Multiplication — canonical chain normalization with explicit scope branching
# ============================================================================

struct _ReductionMode end
const _REDUCE_ONLY = _ReductionMode()

struct _NormState
    prefactor::CNum
    bound::Vector{Index}
    ne::Vector{NonEqualPair}
    chains::Vector{_SpaceChain}
end

function _local_rewrite(c::CNum, ops::Vector{QSym}, ord::OrderingConvention)
    return _apply_ordering(c, ops, ord)
end
_local_rewrite(c::CNum, ops::Vector{QSym}, ::_ReductionMode) = _apply_reductions(c, ops)

function _copy_state(state::_NormState)
    return _NormState(state.prefactor, _copy_bound(state.bound), _copy_ne(state.ne), _copy_chains(state.chains))
end

function _substitute_factor(f::_SiteFactor, from::Index, to::Index)
    new_index = f.index == from ? to : f.index
    new_ops = QSym[change_index(op, from, to) for op in f.ops]
    return _SiteFactor(new_index, new_ops)
end

function _substitute_chain(chain::_SpaceChain, from::Index, to::Index)
    return _SpaceChain(chain.space_index, [_substitute_factor(f, from, to) for f in chain.factors])
end

function _substitute_state(state::_NormState, from::Index, to::Index)
    from_bound = _is_bound_index(state.bound, from)
    to_bound = _is_bound_index(state.bound, to)
    new_bound = if from_bound && !to_bound
        _drop_bound_with(state.bound, from)
    else
        _substitute_bound(state.bound, from, to)
    end
    return _NormState(
        change_index(state.prefactor, from, to),
        new_bound,
        _substitute_ne(state.ne, from, to),
        [_substitute_chain(chain, from, to) for chain in state.chains],
    )
end

function _factor_out_of_order(left::_SiteFactor, right::_SiteFactor)
    return _site_order_key(right.index) < _site_order_key(left.index)
end

function _sites_known_disequal(ne::Vector{NonEqualPair}, left::Index, right::Index)
    left == right && return false
    (!has_index(left) || !has_index(right)) && return true
    return _ne_contains(ne, left, right)
end

function _binding_slot(bound::Vector{Index}, idx::Index)
    for (i, x) in enumerate(bound)
        x == idx && return i
    end
    return 0
end

_is_bound_index(bound::Vector{Index}, idx::Index) = _binding_slot(bound, idx) != 0

function _choose_unify(left::Index, right::Index, bound::Vector{Index})
    left_slot = _binding_slot(bound, left)
    right_slot = _binding_slot(bound, right)
    if left_slot != 0 && right_slot == 0
        return left, right
    elseif right_slot != 0 && left_slot == 0
        return right, left
    elseif left_slot != 0 && right_slot != 0
        return left_slot < right_slot ? (right, left) : (left, right)
    else
        return isless(left, right) ? (right, left) : (left, right)
    end
end

function _replace_factor_pair(
        factors::Vector{_SiteFactor}, i::Int, idx::Index, ops::Vector{QSym}
    )
    out = _SiteFactor[]
    sizehint!(out, length(factors) - 1 + (isempty(ops) ? 0 : 1))
    for k in 1:(i - 1)
        push!(out, _copy_factor(factors[k]))
    end
    isempty(ops) || push!(out, _SiteFactor(idx, copy(ops)))
    for k in (i + 2):length(factors)
        push!(out, _copy_factor(factors[k]))
    end
    return out
end

function _replace_single_factor(
        factors::Vector{_SiteFactor}, i::Int, idx::Index, ops::Vector{QSym}
    )
    out = _SiteFactor[]
    sizehint!(out, length(factors) - 1 + (isempty(ops) ? 0 : 1))
    for k in 1:(i - 1)
        push!(out, _copy_factor(factors[k]))
    end
    isempty(ops) || push!(out, _SiteFactor(idx, copy(ops)))
    for k in (i + 1):length(factors)
        push!(out, _copy_factor(factors[k]))
    end
    return out
end

function _swap_factor_pair(factors::Vector{_SiteFactor}, i::Int)
    out = [_copy_factor(f) for f in factors]
    out[i], out[i + 1] = out[i + 1], out[i]
    return out
end

function _with_chain_replaced(state::_NormState, chain_idx::Int, new_factors::Vector{_SiteFactor})
    chains = _copy_chains(state.chains)
    if isempty(new_factors)
        deleteat!(chains, chain_idx)
    else
        chains[chain_idx] = _SpaceChain(chains[chain_idx].space_index, new_factors)
    end
    return _NormState(state.prefactor, _copy_bound(state.bound), _copy_ne(state.ne), chains)
end

function _normalize_state(state::_NormState, mode)
    for chain_idx in eachindex(state.chains)
        chain = state.chains[chain_idx]
        factors = chain.factors
        for factor_idx in eachindex(factors)
            factor = factors[factor_idx]
            length(factor.ops) <= 1 && continue
            local_terms = _local_rewrite(state.prefactor, factor.ops, mode)
            if length(local_terms) == 1
                only_term = only(local_terms)
                if isequal(only_term.prefactor, state.prefactor) && isequal(only_term.ops, factor.ops)
                    continue
                end
            end
            out = _NormState[]
            for t in local_terms
                next_state = _with_chain_replaced(
                    _NormState(t.prefactor, state.bound, state.ne, state.chains),
                    chain_idx,
                    _replace_single_factor(factors, factor_idx, factor.index, t.ops),
                )
                append!(out, _normalize_state(next_state, mode))
            end
            return out
        end
        length(factors) <= 1 && continue
        for i in 1:(length(factors) - 1)
            left = factors[i]
            right = factors[i + 1]

            if left.index == right.index
                merged = vcat(left.ops, right.ops)
                out = _NormState[]
                for t in _local_rewrite(state.prefactor, merged, mode)
                    next_state = _with_chain_replaced(
                        _NormState(t.prefactor, state.bound, state.ne, state.chains),
                        chain_idx,
                        _replace_factor_pair(factors, i, left.index, t.ops),
                    )
                    append!(out, _normalize_state(next_state, mode))
                end
                return out
            end

            if _sites_known_disequal(state.ne, left.index, right.index)
                if _factor_out_of_order(left, right)
                    swapped = _with_chain_replaced(state, chain_idx, _swap_factor_pair(factors, i))
                    return _normalize_state(swapped, mode)
                end
                continue
            end

            left_bound = _is_bound_index(state.bound, left.index)
            right_bound = _is_bound_index(state.bound, right.index)
            if !(left_bound || right_bound)
                continue
            end

            eq_from, eq_to = _choose_unify(left.index, right.index, state.bound)
            eq_state = _substitute_state(state, eq_from, eq_to)

            diseq_state = _NormState(
                state.prefactor,
                _copy_bound(state.bound),
                _merge_ne_pair(state.ne, left.index, right.index),
                _copy_chains(state.chains),
            )
            if _factor_out_of_order(left, right)
                diseq_state = _with_chain_replaced(diseq_state, chain_idx, _swap_factor_pair(factors, i))
            end

            out = _NormState[]
            append!(out, _normalize_state(eq_state, mode))
            append!(out, _normalize_state(diseq_state, mode))
            return out
        end
    end
    return _NormState[state]
end

function _emit_state!(d::QTermDict, state::_NormState)
    ops = _canonical_ops_from_chains(state.chains)
    bound, ne = _normalize_scope(state.prefactor, ops, state.bound, state.ne)
    _addto_key!(d, _term_key(state.chains, bound, ne), state.prefactor)
    return d
end

function _accumulate_normalized!(
        d::QTermDict, c::CNum, ops::Vector{QSym},
        bound::Vector{Index}, ne::Vector{NonEqualPair}, mode
    )
    _iszero_cnum(c) && return d
    state = _NormState(c, _copy_bound(bound), _copy_ne(ne), _chains_from_ops(ops))
    for final_state in _normalize_state(state, mode)
        _emit_state!(d, final_state)
    end
    return d
end

function Base.:*(a::QSym, b::QSym)
    d = QTermDict()
    _accumulate_normalized!(d, _CNUM_ONE, QSym[a, b], _empty_bound(), _empty_ne(), get_ordering())
    return _qadd(d)
end

Base.:*(a::QSym, b::Number) = _single_qadd(_to_cnum(b), QSym[a])
Base.:*(b::Number, a::QSym) = a * b

function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict()
    for (term, c) in a.arguments
        new_c = _mul_cnum(c, cb)
        _iszero_cnum(new_c) && continue
        d[_copy_key(term)] = new_c
    end
    return _qadd(d, copy(a.indices))
end
Base.:*(a::Number, b::QAdd) = b * a

function _term_indices(term::QTerm)
    out = copy(term.bound)
    for op in term.ops
        has_index(op.index) && _push_index_unique!(out, op.index)
    end
    for (α, β) in term.ne
        _push_index_unique!(out, α)
        _push_index_unique!(out, β)
    end
    return out
end

function _fresh_index(idx::Index, taken::Vector{Index})
    n = 1
    while true
        name = Symbol(idx.name, "_", n)
        sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
        candidate = Index(name, idx.range, idx.space_index, Num(sym_var))
        candidate ∉ taken && return candidate
        n += 1
    end
end

function _alpha_rename_collisions(term::QTerm, c::CNum, taken::Vector{Index})
    new_c = c
    new_chains = _copy_chains(term.chains)
    new_bound = copy(term.bound)
    new_ne = copy(term.ne)
    seen = copy(taken)
    changed = false
    for idx in term.bound
        idx in seen || continue
        fresh = _fresh_index(idx, seen)
        new_c = change_index(new_c, idx, fresh)
        new_chains = [_substitute_chain(chain, idx, fresh) for chain in new_chains]
        new_bound = _substitute_bound(new_bound, idx, fresh)
        new_ne = _substitute_ne(new_ne, idx, fresh)
        push!(seen, fresh)
        changed = true
    end
    changed || return term, c
    return _term_key(new_chains, new_bound, new_ne), new_c
end

function Base.:*(a::QAdd, b::QSym)
    d = QTermDict()
    for (term, c) in a.arguments
        _accumulate_normalized!(d, c, vcat(_flatten_chains(term.chains), QSym[b]), term.bound, term.ne, get_ordering())
    end
    return _qadd(d)
end

function Base.:*(a::QSym, b::QAdd)
    d = QTermDict()
    for (term, c) in b.arguments
        _accumulate_normalized!(d, c, vcat(QSym[a], _flatten_chains(term.chains)), term.bound, term.ne, get_ordering())
    end
    return _qadd(d)
end

function Base.:*(a::QAdd, b::QAdd)
    d = QTermDict()
    for (term_a, c_a) in a.arguments, (term_b, c_b) in b.arguments
        term_b_eff, c_b_eff = _alpha_rename_collisions(term_b, c_b, _term_indices(term_a))
        c_prod = _mul_cnum(c_a, c_b_eff)
        inherited = isempty(term_a.ne) ? term_b_eff.ne :
            (isempty(term_b_eff.ne) ? term_a.ne : _merge_ne(term_a.ne, term_b_eff.ne))
        bound = isempty(term_a.bound) ? term_b_eff.bound :
            (isempty(term_b_eff.bound) ? term_a.bound : _merge_bound(term_a.bound, term_b_eff.bound))
        _accumulate_normalized!(
            d,
            c_prod,
            vcat(_flatten_chains(term_a.chains), _flatten_chains(term_b_eff.chains)),
            bound,
            inherited,
            get_ordering(),
        )
    end
    return _qadd(d)
end

# ============================================================================
#  Addition — always returns QAdd
# ============================================================================

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[b], _CNUM_ONE)
    return _qadd(d)
end

function Base.:+(a::QAdd, b::QSym)
    d = _copy_args(a.arguments)
    _addto!(d, QSym[b], _CNUM_ONE)
    return _qadd(d, copy(a.indices))
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = _copy_args(a.arguments)
    for (term, c) in b.arguments
        _addto_key!(d, _copy_key(term), c)
    end
    return _qadd(d, _merge_unique(a.indices, b.indices))
end

function Base.:+(a::QSym, b::Number)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[], _to_cnum(b))
    return _qadd(d)
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QAdd, b::Number)
    d = _copy_args(a.arguments)
    _addto!(d, QSym[], _to_cnum(b))
    return _qadd(d, copy(a.indices))
end
Base.:+(a::Number, b::QAdd) = b + a

Base.zero(::Type{QAdd}) = _zero_qadd()
Base.zero(::QAdd) = _zero_qadd()

# ============================================================================
#  Subtraction, negation, division, power
# ============================================================================

Base.:-(a::QSym) = _single_qadd(_CNUM_NEG1, QSym[a])

function Base.:-(a::QAdd)
    d = QTermDict()
    for (term, c) in a.arguments
        d[_copy_key(term)] = _neg_cnum(c)
    end
    return _qadd(d, copy(a.indices))
end

Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

Base.:/(a::QSym, b::Number) = a * inv(b)
Base.:/(a::QAdd, b::Number) = a * inv(b)

Base.://(a::QSym, b::Integer) = a * (1 // b)
Base.://(a::QAdd, b::Integer) = a * (1 // b)

function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, QSym[])
    d = QTermDict()
    _accumulate_normalized!(d, _CNUM_ONE, QSym[a for _ in 1:n], _empty_bound(), _empty_ne(), get_ordering())
    return _qadd(d)
end

function Base.:^(a::QAdd, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, QSym[])
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end
