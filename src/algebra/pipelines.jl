"""
    stream!(out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair})

The default canonicalization pipeline. Order: `reduce -> commute -> reduce -> sink`.

The first reduce folds Transition/Pauli same-site pairs (these never commute,
they compose). Commute then operates only on Fock/Spin/PhaseSpace ladder pairs.
The trailing reduce catches any same-site composition surfaced by the commute
residual (e.g. a Spin commutator's contracted op meeting a same-site neighbor).
"""
@inline function stream!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
    )
    reduce_ops(ops, c) do ops1, c1
        commute_ops(ops1, c1) do ops2, c2
            reduce_ops(ops2, c2) do ops3, c3
                canonicalize_to_dict!(out, ops3, c3, ne)
            end
        end
    end
    return out
end

"""
    canonicalize!(out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair})

Re-establish canonical form: sort, then run the pipeline.
"""
@inline function canonicalize!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
    )
    partial_sort!(ops, ne)
    stream!(out, ops, c, ne)
    return out
end

"""
    emit_product!(out, ta_ops, ca, ta_ne, tb_ops, cb, tb_ne, sum_indices, needs_diag)

Concatenate the operator strings of two terms, merge their constraint sets, and
route through `canonicalize!` (or `accumulate_with_diag!` when summation
indices may collapse).
"""
@inline function emit_product!(
        out::QTermDict,
        ta_ops::Vector{Op}, ca::CNum, ta_ne::Vector{NonEqualPair},
        tb_ops::Vector{Op}, cb::CNum, tb_ne::Vector{NonEqualPair},
        sum_indices::Vector{Index}, needs_diag_split::Bool,
    )
    n = length(ta_ops) + length(tb_ops)
    ops = Vector{Op}(undef, n)
    copyto!(ops, 1, ta_ops, 1, length(ta_ops))
    copyto!(ops, length(ta_ops) + 1, tb_ops, 1, length(tb_ops))
    ne = merge_ne(ta_ne, tb_ne)
    c = mul_cnum(ca, cb)
    if needs_diag_split
        accumulate_with_diag!(out, ops, c, sum_indices, ne)
    else
        canonicalize!(out, ops, c, ne)
    end
    return out
end

function distinct_op_indices(ops::Vector{Op})
    out = Index[]
    for op in ops
        idx = op.index
        has_index(idx) || continue
        idx in out && continue
        push!(out, idx)
    end
    return out
end

function depends_on_index_ops(c::CNum, ops::Vector{Op}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    return depends_on_index_term(c, ops, idx)
end

# When both indices are bound, the two directions (i→j and j→i) describe
# the same diagonal slice; keeping one (later→earlier) avoids double-counting
# constant residuals (the `+1` from `a*a' = a'a + 1`) and keeps the survivor
# invariant under factor order so `Σa * Σa'` and `Σa' * Σa` agree.
function diagonal_pairs(
        ops::Vector{Op}, c::CNum, sum_indices::Vector{Index},
        ne::Vector{NonEqualPair}, distinct::Vector{Index},
    )
    diag_pairs = Tuple{Index, Index}[]
    for sum_idx in sum_indices
        depends_on_index_ops(c, ops, sum_idx) || continue
        for ext_idx in distinct
            may_coincide(sum_idx, ext_idx, ne) || continue
            # Keep one direction of a bound-bound pair.
            if ext_idx in sum_indices
                name_rank(sum_idx.name_id) > name_rank(ext_idx.name_id) || continue
            end
            push!(diag_pairs, (sum_idx, ext_idx))
        end
    end
    return diag_pairs
end

# Two distinct indices on the same space, not already constrained apart.
function may_coincide(sum_idx::Index, ext_idx::Index, ne::Vector{NonEqualPair})
    ext_idx == sum_idx && return false
    ext_idx.space_index == sum_idx.space_index || return false
    return !ne_contains(ne, sum_idx, ext_idx)
end

function sibling_exclusions(diag_pairs::Vector{Tuple{Index, Index}}, sum_idx::Index, ext_idx::Index)
    sibling_ne = NonEqualPair[]
    for (other_sum, other_ext) in diag_pairs
        other_sum == sum_idx || continue
        other_ext == ext_idx && continue
        other_ext.space_index == ext_idx.space_index || continue
        push!(sibling_ne, (ext_idx, other_ext))
    end
    return sibling_ne
end

"""
    accumulate_with_diag!(out, ops, c, sum_indices, ne) -> out

When summing over indices that may coincide with another operator's index, emit
both the off-diagonal contribution (under `ne ∪ {(sum_idx, ext_idx)}` enforcing
the indices differ) and each diagonal contribution (substituting
`sum_idx -> ext_idx`, with each constraint on `sum_idx` rewritten onto
`ext_idx`, mirroring the diagonal split in [`Σ`](@ref)).
"""
function accumulate_with_diag!(
        out::QTermDict, ops::Vector{Op}, c::CNum,
        sum_indices::Vector{Index}, ne::Vector{NonEqualPair},
    )
    distinct = distinct_op_indices(ops)
    if length(distinct) < 2
        canonicalize!(out, copy(ops), c, ne)
        return out
    end

    diag_pairs = diagonal_pairs(ops, c, sum_indices, ne, distinct)

    if isempty(diag_pairs)
        canonicalize!(out, copy(ops), c, ne)
        return out
    end

    augmented_ne = ne
    for (a, b) in diag_pairs
        augmented_ne = merge_ne_pair(augmented_ne, a, b)
    end
    canonicalize!(out, copy(ops), c, augmented_ne)

    for (sum_idx, ext_idx) in diag_pairs
        sub_ops = Op[change_index(o, sum_idx, ext_idx) for o in ops]
        sub_c = change_index(c, sum_idx, ext_idx)
        sub_ne = substitute_ne(ne, sum_idx, ext_idx)
        sibling_ne = sibling_exclusions(diag_pairs, sum_idx, ext_idx)
        isempty(sibling_ne) || (sub_ne = merge_ne(sub_ne, sibling_ne))
        if ext_idx in sum_indices
            emit_scaled_by_scope!(out, sub_ops, sub_c, sub_ne, Index[ext_idx])
        else
            canonicalize!(out, sub_ops, sub_c, sub_ne)
        end
    end
    return out
end

"""
    emit_scaled_by_scope!(out, ops, c, ne, scope)

Canonicalize `(ops, c, ne)` into `out`, then settle each `scope` index through
`settle_scope!`. Use when residuals that lose a bound-index dependence must absorb the
surviving sum's count (e.g. the `+1` from `a*a' = a'a + 1` inside a `Σ_i` scope).
"""
function emit_scaled_by_scope!(
        out::QTermDict, ops::Vector{Op}, c::CNum,
        ne::Vector{NonEqualPair}, scope::Vector{Index},
    )
    temp = QTermDict()
    canonicalize!(temp, ops, c, ne)
    return settle_scope!(out, temp, scope)
end

"""
    settle_scope!(out, temp, scope)

Add the terms of `temp`, all descended from one term bound by `scope`, into `out`. Each
scope index that a term's operators and coefficient no longer mention is eliminated by
`eliminate_binder`: the term is multiplied by the index's constrained count
(`range - |partners|`), or keeps the binder and its constraints when that count is unknown.
"""
function settle_scope!(out::QTermDict, temp::QTermDict, scope::Vector{Index})
    for (term, cv) in temp
        coef = cv
        ne = term.ne
        for scope_idx in scope
            depends_on_index_term(coef, term.ops, scope_idx) && continue
            eliminated, factor, ne = eliminate_binder(ne, scope_idx)
            eliminated && (coef = mul_cnum(coef, factor))
        end
        key = ne === term.ne ? copy_key(term) : QTerm(copy(term.ops), ne)
        addto_key!(out, key, coef)
    end
    return out
end

function Base.:*(a::QSym, b::QSym)
    out = QTermDict()
    emit_product!(
        out,
        Op[a], CNUM_ONE, EMPTY_NE,
        Op[b], CNUM_ONE, EMPTY_NE,
        EMPTY_INDICES, false
    )
    return QAdd(out, EMPTY_INDICES)
end

function Base.:*(a::QAdd, b::QSym)
    out = QTermDict()
    needs = !isempty(a.indices)
    for (ta, ca) in a
        emit_product!(
            out, ta.ops, ca, ta.ne,
            Op[b], CNUM_ONE, EMPTY_NE,
            a.indices, needs
        )
    end
    return QAdd(out, absorb_pinned_sums(a.indices, a, b))
end

function Base.:*(a::QSym, b::QAdd)
    out = QTermDict()
    needs = !isempty(b.indices)
    for (tb, cb) in b
        emit_product!(
            out, Op[a], CNUM_ONE, EMPTY_NE,
            tb.ops, cb, tb.ne,
            b.indices, needs
        )
    end
    return QAdd(out, absorb_pinned_sums(b.indices, a, b))
end

function Base.:*(a::QAdd, b::QAdd)
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            idx in b.indices && throw(
                ArgumentError(
                    "Summation index $(index_name(idx)) appears in both factors. " *
                        "Bound variables on the two sides of a product must be " *
                        "distinct so the resulting double sum is unambiguous. " *
                        "Use `change_index` to rename one side, or construct the " *
                        "two sums with different indices."
                )
            )
        end
    end
    out = QTermDict()
    sum_indices = merge_unique(a.indices, b.indices)
    needs = !isempty(sum_indices)
    for (ta, ca) in a, (tb, cb) in b
        emit_product!(
            out, ta.ops, ca, ta.ne,
            tb.ops, cb, tb.ne,
            sum_indices, needs
        )
    end
    return QAdd(out, absorb_pinned_sums(sum_indices, a, b))
end

function absorb_pinned_sums(
        sum_indices::Vector{Index}, a::Union{QAdd, QSym}, b::Union{QAdd, QSym}
    )
    isempty(sum_indices) && return sum_indices
    a_indices = a isa QAdd ? a.indices : EMPTY_INDICES
    b_indices = b isa QAdd ? b.indices : EMPTY_INDICES
    keep = Index[]
    @inbounds for s in sum_indices
        from_a = s in a_indices
        from_b = s in b_indices
        pinned = (from_a && !from_b && every_term_has_op_index(b, s)) ||
            (from_b && !from_a && every_term_has_op_index(a, s))
        pinned || push!(keep, s)
    end
    return length(keep) == length(sum_indices) ? sum_indices : keep
end

every_term_has_op_index(q::QSym, idx::Index) = q.index == idx
function every_term_has_op_index(q::QAdd, idx::Index)
    isempty(q.arguments) && return false
    @inbounds for t in keys(q.arguments)
        found = false
        for op in t.ops
            if op.index == idx
                found = true
                break
            end
        end
        found || return false
    end
    return true
end
