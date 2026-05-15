# --- Algebraic reductions (ordering-independent) ---

"""
    _apply_reductions(arg_c, ops) -> Vector{OrderedTerm}

Apply ordering-independent algebraic identities to a product. Only fires
Transition composition (`|i⟩⟨j| · |k⟩⟨l| → δⱼₖ |i⟩⟨l|`) and the Pauli product
rule (`σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ`); never reorders via commutation. Drives the
shared [`_process_product!`](@ref) under `LazyOrder` (skips swaps) with
`expand_gs = false` (keeps `σᵍᵍ` atomic — `simplify(expr, h)` opts back in).
"""
_apply_reductions(arg_c::CNum, ops::Vector{QSym}) =
    _drive_worklist(arg_c, ops, LazyOrder(), false)

function _extra_scope_pair(base::QTerm, constrained::QTerm)
    isequal(base.ops, constrained.ops) || return nothing
    isequal(base.bound, constrained.bound) || return nothing
    length(constrained.ne) == length(base.ne) + 1 || return nothing
    for p in constrained.ne
        p in base.ne && continue
        trial = _canonical_ne(_merge_ne_pair(base.ne, p[1], p[2]))
        isequal(trial, constrained.ne) && return p
    end
    return nothing
end

function _is_transition_projector_factor(f::_SiteFactor)
    length(f.ops) == 1 || return false
    op = only(f.ops)
    return op isa Transition && op.i == op.j
end

function _single_transition(f::_SiteFactor)
    length(f.ops) == 1 || return nothing
    op = only(f.ops)
    return op isa Transition ? op : nothing
end

function _adjacent_pair_factors(term::QTerm, p::NonEqualPair)
    for chain in term.chains
        factors = chain.factors
        length(factors) <= 1 && continue
        for i in 1:(length(factors) - 1)
            left = factors[i]
            right = factors[i + 1]
            if (left.index == p[1] && right.index == p[2]) ||
                    (left.index == p[2] && right.index == p[1])
                return left, right
            end
        end
    end
    return nothing
end

function _choose_split_unify(term::QTerm, p::NonEqualPair)
    if !isempty(term.bound)
        return _choose_unify(p[1], p[2], term.bound)
    end

    pair = _adjacent_pair_factors(term, p)
    if pair !== nothing
        left, right = pair
        left_proj = _is_transition_projector_factor(left)
        right_proj = _is_transition_projector_factor(right)
        if xor(left_proj, right_proj)
            return left_proj ? (left.index, right.index) : (right.index, left.index)
        end

        left_t = _single_transition(left)
        right_t = _single_transition(right)
        if left_t !== nothing && right_t !== nothing &&
                left_t.name == right_t.name &&
                left_t.j == right_t.i &&
                left_t.i == right_t.j
            # When the diagonal branch collapses to a projector, keep the
            # right factor's site label so later completeness terms line up.
            return left.index, right.index
        end
    end

    return _choose_unify(p[1], p[2], term.bound)
end

function _suppress_free_transition_projector_diag(term::QTerm, p::NonEqualPair)
    isempty(term.bound) || return false
    pair = _adjacent_pair_factors(term, p)
    pair === nothing && return false
    left, right = pair
    left_t = _single_transition(left)
    right_t = _single_transition(right)
    left_t === nothing && return false
    right_t === nothing && return false
    left_t.name == right_t.name || return false
    return left_t.j == right_t.i && left_t.i == right_t.j
end

function _split_term_on_pair(term::QTerm, c::CNum, p::NonEqualPair, mode)
    offdiag = QTermDict()
    diag = QTermDict()
    _accumulate_normalized!(
        offdiag,
        c,
        _flatten_chains(term.chains),
        term.bound,
        _merge_ne_pair(term.ne, p[1], p[2]),
        mode,
    )
    _suppress_free_transition_projector_diag(term, p) && return offdiag, diag

    from, to = _choose_split_unify(term, p)
    state = _NormState(c, _copy_bound(term.bound), _copy_ne(term.ne), _copy_chains(term.chains))
    diag_state = _substitute_state(state, from, to)
    _accumulate_normalized!(
        diag,
        diag_state.prefactor,
        _flatten_chains(diag_state.chains),
        diag_state.bound,
        diag_state.ne,
        mode,
    )
    return offdiag, diag
end

function _free_pair_to_branch(term::QTerm, c::CNum)
    isempty(term.bound) || return nothing
    for chain in term.chains
        factors = chain.factors
        length(factors) <= 1 && continue
        for i in 1:(length(factors) - 1)
            left = factors[i]
            right = factors[i + 1]
            left.index == right.index && continue
            _sites_known_disequal(term.ne, left.index, right.index) && continue
            _depends_on_index_term(c, term.ops, left.index) || continue
            _depends_on_index_term(c, term.ops, right.index) || continue
            return (left.index, right.index)
        end
    end
    return nothing
end

function _resolve_scope_terms!(d::QTermDict, mode)
    changed = true
    while changed
        changed = false
        terms = collect(keys(d))
        for base in terms
            haskey(d, base) || continue
            for constrained in terms
                base === constrained && continue
                haskey(d, constrained) || continue
                extra = _extra_scope_pair(base, constrained)
                extra === nothing && continue

                base_c = d[base]
                _depends_on_index_term(base_c, base.ops, extra[1]) || continue
                _depends_on_index_term(base_c, base.ops, extra[2]) || continue
                constrained_c = d[constrained]
                delete!(d, base)
                delete!(d, constrained)

                offdiag, diag = _split_term_on_pair(base, base_c, extra, mode)
                for (term, c) in offdiag
                    _addto_key!(d, term, c)
                end
                for (term, c) in diag
                    _addto_key!(d, term, c)
                end
                _addto_key!(d, constrained, constrained_c)
                changed = true
                break
            end
            changed && break
        end
    end
    return d
end

function _resolve_free_pair_terms!(d::QTermDict, mode)
    changed = true
    while changed
        changed = false
        terms = collect(keys(d))
        for term in terms
            haskey(d, term) || continue
            c = d[term]
            pair = _free_pair_to_branch(term, c)
            pair === nothing && continue

            delete!(d, term)
            offdiag, diag = _split_term_on_pair(term, c, pair, mode)
            for (next_term, next_c) in offdiag
                _addto_key!(d, next_term, next_c)
            end
            for (next_term, next_c) in diag
                _addto_key!(d, next_term, next_c)
            end
            changed = true
            break
        end
    end
    return d
end

function _simplify_num_part(x::Num)
    s0 = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(x)))
    _iszero_num(s0) && return _NUM_ZERO
    sx = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(Symbolics.expand(s0))))
    return _iszero_num(sx) ? _NUM_ZERO : sx
end

function _simplify_prefactor(x::CNum)
    re = _simplify_num_part(real(x))
    im = _simplify_num_part(imag(x))
    return Complex(re, im)
end
_simplify_prefactor(x::Number) = x

# --- Simplification (ordering-independent) ---

function _needs_scope_resolution(d::QTermDict)
    for term in keys(d)
        isempty(term.bound) || return true
        isempty(term.ne) || return true
        for op in term.ops
            has_index(op.index) && return true
        end
    end
    return false
end

function _qsimplify(op::QSym)
    return _single_qadd(_CNUM_ONE, QSym[op])
end

function _qsimplify(s::QAdd)
    reduced_d = QTermDict()
    for (term, c) in s.arguments
        _iszero_cnum(c) && continue
        _accumulate_normalized!(reduced_d, c, _flatten_chains(term.chains), term.bound, term.ne, _REDUCE_ONLY)
    end
    if _needs_scope_resolution(reduced_d)
        changed = true
        while changed
            before = _copy_args(reduced_d)
            _resolve_scope_terms!(reduced_d, get_ordering())
            _resolve_free_pair_terms!(reduced_d, get_ordering())
            changed = !isequal(before, reduced_d)
        end
    end

    d = QTermDict()
    for (term, c) in reduced_d
        new_c = _simplify_prefactor(c)
        _iszero_cnum(new_c) && continue
        _addto_key!(d, _copy_key(term), new_c)
    end
    # Drop summation indices that no surviving term depends on.
    indices = _drop_unused_indices(d, s.indices)
    return _qadd(d, indices)
end

"""
    _drop_unused_indices(d, indices) -> Vector{Index}

Filter `indices` to those some term in `d` actually depends on. Algebraic
reductions can collapse all index-bearing terms (e.g. `σ¹² · σ³¹ = 0`),
leaving the surrounding `Σ` over a now-vacuous index — this drops it so
the resulting `QAdd` doesn't claim a phantom summation.
"""
function _drop_unused_indices(d::QTermDict, indices::Vector{Index})
    isempty(indices) && return indices
    used = Index[]
    for idx in indices
        for (term, c) in d
            if idx in term.bound || _depends_on_index_term(c, term.ops, idx)
                push!(used, idx)
                break
            end
        end
    end
    length(used) == length(indices) && return indices
    return used
end

# With Hilbert space (ground state rewriting)
function _qsimplify(op::QSym, h::HilbertSpace)
    return _apply_ground_state(_qsimplify(op), h)
end
function _qsimplify(s::QAdd, h::HilbertSpace)
    return _apply_ground_state(_qsimplify(s), h)
end

# Public API

"""
    simplify(expr::QField)
    simplify(expr::QField, h::HilbertSpace)

Simplify a quantum operator expression by applying ordering-independent algebraic identities.

Applied rules:
- **Transition composition**: ``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``
- **Pauli product rule**: ``\\sigma_j \\sigma_k = \\delta_{jk} I + i\\epsilon_{jkl} \\sigma_l``
- Like-term collection and zero-term elimination

Does **not** apply commutation-based reordering (Fock `[a, a†]`, Spin `[Sⱼ, Sₖ]`,
PhaseSpace `[p, x]`). For that, use [`normal_order`](@ref).

The `h`-argument overload is the [`LazyOrder`](@ref) opt-in for ground-state
completeness: it additionally rewrites ground-state projectors
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|`` of every
[`NLevelSpace`](@ref) subspace. Without `h`, ``|g\\rangle\\langle g|`` is kept
atomic — useful for inspecting the un-expanded form.

Under [`NormalOrder`](@ref) (the default ordering convention) the eager arithmetic
already produces canonical form including completeness, so `simplify` is typically
idempotent and the `h`-overload is a no-op cleanup pass.

See also [`normal_order`](@ref), [`expand`](@ref).
"""
SymbolicUtils.simplify(expr::QField; kwargs...) = _qsimplify(expr)
SymbolicUtils.simplify(expr::QField, h::HilbertSpace; kwargs...) = _qsimplify(expr, h)

"""
    expand(expr::QField)

Expand symbolic prefactors in quantum operator expressions.

Applies `Symbolics.expand` to each prefactor in a [`QAdd`](@ref) sum,
e.g. expanding `(a+b)^2` into `a^2 + 2ab + b^2`. Zero terms are dropped after expansion.

See also [`simplify`](@ref).
"""
function Symbolics.expand(s::QAdd; kwargs...)
    d = QTermDict()
    for (term, c) in s.arguments
        new_c = _expand_prefactor(c; kwargs...)
        _iszero_cnum(new_c) && continue
        _addto_key!(d, _copy_key(term), new_c)
    end
    return _qadd(d, copy(s.indices))
end
function Symbolics.expand(op::QSym; kwargs...)
    return _single_qadd(_CNUM_ONE, QSym[op])
end

_expand_prefactor(x::CNum; kwargs...) = _iszero_cnum(x) ? x : Symbolics.expand(x; kwargs...)
_expand_prefactor(x::Number; kwargs...) = x
