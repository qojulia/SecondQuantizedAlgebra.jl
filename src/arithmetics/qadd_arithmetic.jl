# ============================================================================
#  Multiplication — always returns QAdd, eagerly ordered
# ============================================================================

function Base.:*(a::QSym, b::QSym)
    phys_ops = QSym[a, b]
    _phys_sort!(phys_ops)
    ops = copy(phys_ops)
    _site_sort!(ops)
    return _ordered_qadd(_CNUM_ONE, ops, _EMPTY_NE, phys_ops)
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
    return QAdd(d, copy(a.indices))
end
Base.:*(a::Number, b::QAdd) = b * a

function _apply_diag_terms!(
        d::QTermDict, terms::Vector{OrderedTerm}, ne::Vector{NonEqualPair},
        phys_ops::Vector{QSym}
    )
    for t in terms
        _addto!(d, t.ops, t.prefactor, ne, phys_ops)
    end
    return d
end

"""
    _substitute_unsorted(c, ops, α, β) -> Vector{OrderedTerm}

Substitute `α → β` on the *unsorted* operator sequence `ops` (the original
physical order before [`_site_sort!`](@ref)) and then site-sort + apply
ordering. Pre-sort substitution preserves the physical operator order so
that same-site composition rules (e.g. ``\\sigma_{ij}\\sigma_{jk} =
\\sigma_{ik}``) fire on the correct adjacency. See `devdocs.md` "Diagonal
splitting" for why post-sort substitution is wrong.
"""
function _substitute_unsorted(
        c::CNum, ops::Vector{QSym}, α::Index, β::Index
    )
    sub_ops = QSym[change_index(o, α, β) for o in ops]
    sub_c = change_index(c, α, β)
    _site_sort!(sub_ops)
    return _apply_ordering(sub_c, sub_ops)
end

function _distinct_op_indices(ops::Vector{QSym})
    out = Index[]
    for op in ops
        idx = op.index
        has_index(idx) || continue
        idx in out && continue
        push!(out, idx)
    end
    return out
end

"""
    _emit_diagonal!(d, diag_terms, off_diag_terms, α, β) -> Vector{QTerm}

Add the diagonal contribution `diag_terms` (from `α = β`) into `d` using the
shared `ne` of the off-diagonal entries, then re-key those off-diagonal entries
under the augmented constraint `ne ∪ {(α, β)}`. Returns the new key list. All
entries in `off_diag_terms` are assumed to share the same `ne` — an invariant
maintained by [`_accumulate_with_diag!`](@ref).
"""
function _emit_diagonal!(
        d::QTermDict, diag_terms::Vector{OrderedTerm},
        off_diag_terms::Vector{QTerm}, diag_phys_ops::Vector{QSym},
        α::Index, β::Index
    )
    isempty(off_diag_terms) && return off_diag_terms
    current_ne = first(off_diag_terms).ne
    _apply_diag_terms!(d, diag_terms, _drop_ne_with(current_ne, α), diag_phys_ops)

    new_ne = _merge_ne_pair(current_ne, α, β)
    moved = QTerm[]
    for old_term in off_diag_terms
        c = get(d, old_term, nothing)
        c === nothing && continue
        delete!(d, old_term)
        new_term = _term_key(old_term.ops, new_ne, old_term.phys_ops)
        _addto_key!(d, new_term, c)
        _push_key_unique!(moved, new_term)
    end
    return moved
end

"""
    _accumulate_with_diag!(d, c, unsorted_ops, sum_idxes, inherited_ne)

Accumulate the off-diagonal ordered contribution from `(c, unsorted_ops)` into
`d`, then emit any required diagonal substitutions using the unsorted physical
operator order.
"""
function _accumulate_with_diag!(
        d::QTermDict, c::CNum, unsorted_ops::Vector{QSym},
        sum_idxes::Vector{Index},
        inherited_ne::Vector{NonEqualPair}
    )
    phys_ops = copy(unsorted_ops)
    _phys_sort!(phys_ops)
    sorted = copy(unsorted_ops)
    _site_sort!(sorted)
    sorted_terms = _apply_ordering(c, sorted)

    distinct = _distinct_op_indices(unsorted_ops)
    if length(distinct) < 2
        for t in sorted_terms
            _addto!(d, t.ops, t.prefactor, inherited_ne, phys_ops)
        end
        return nothing
    end

    off_diag_terms = QTerm[]
    for t in sorted_terms
        term = _term_key(t.ops, inherited_ne, phys_ops)
        _addto_key!(d, term, t.prefactor)
        _push_key_unique!(off_diag_terms, term)
    end

    for sum_idx in sum_idxes
        _depends_on_index_term(c, unsorted_ops, sum_idx) || continue
        for ext_idx in distinct
            ext_idx == sum_idx && continue
            ext_idx.space_index == sum_idx.space_index || continue
            _ne_contains(first(off_diag_terms).ne, sum_idx, ext_idx) && continue
            diag_phys_ops = QSym[change_index(o, sum_idx, ext_idx) for o in phys_ops]
            _phys_sort!(diag_phys_ops)
            off_diag_terms = _emit_diagonal!(
                d, _substitute_unsorted(c, unsorted_ops, sum_idx, ext_idx),
                off_diag_terms, diag_phys_ops, sum_idx, ext_idx,
            )
        end
    end

    return nothing
end

function Base.:*(a::QAdd, b::QSym)
    d = QTermDict()
    for (term, c) in a.arguments
        n = length(term.phys_ops)
        unsorted = Vector{QSym}(undef, n + 1)
        copyto!(unsorted, 1, term.phys_ops, 1, n)
        unsorted[n + 1] = b
        _accumulate_with_diag!(d, c, unsorted, a.indices, term.ne)
    end
    return QAdd(d, copy(a.indices))
end

function Base.:*(a::QSym, b::QAdd)
    d = QTermDict()
    for (term, c) in b.arguments
        n = length(term.phys_ops)
        unsorted = Vector{QSym}(undef, n + 1)
        unsorted[1] = a
        copyto!(unsorted, 2, term.phys_ops, 1, n)
        _accumulate_with_diag!(d, c, unsorted, b.indices, term.ne)
    end
    return QAdd(d, copy(b.indices))
end

function Base.:*(a::QAdd, b::QAdd)
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            idx in b.indices && throw(
                ArgumentError(
                    "Summation index $(idx.name) appears in both factors. " *
                        "Use `change_index` to re-index one side, or use different indices."
                )
            )
        end
    end

    d = QTermDict()
    sum_idxes = isempty(b.indices) ? a.indices :
        (isempty(a.indices) ? b.indices : vcat(a.indices, b.indices))

    for (term_a, c_a) in a.arguments, (term_b, c_b) in b.arguments
        unsorted = vcat(term_a.phys_ops, term_b.phys_ops)
        c_prod = _mul_cnum(c_a, c_b)
        inherited = isempty(term_a.ne) ? term_b.ne :
            (isempty(term_b.ne) ? term_a.ne : _merge_ne(term_a.ne, term_b.ne))
        _accumulate_with_diag!(d, c_prod, unsorted, sum_idxes, inherited)
    end

    return QAdd(d, _merge_unique(a.indices, b.indices))
end

# ============================================================================
#  Addition — always returns QAdd
# ============================================================================

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, Index[])
end

function Base.:+(a::QAdd, b::QSym)
    d = _copy_args(a.arguments)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, copy(a.indices))
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = _copy_args(a.arguments)
    for (term, c) in b.arguments
        _addto_key!(d, _copy_key(term), c)
    end
    return QAdd(d, _merge_unique(a.indices, b.indices))
end

function Base.:+(a::QSym, b::Number)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[], _to_cnum(b))
    return QAdd(d, Index[])
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QAdd, b::Number)
    d = _copy_args(a.arguments)
    _addto!(d, QSym[], _to_cnum(b))
    return QAdd(d, copy(a.indices))
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
    return QAdd(d, copy(a.indices))
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
    ops = QSym[a for _ in 1:n]
    return _ordered_qadd(_CNUM_ONE, ops)
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
