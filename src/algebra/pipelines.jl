"""
    stream!(out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair})

The default canonicalization pipeline: reduce, commute, reduce, then store in `out`.

The first reduce composes same-site Transition and Pauli pairs, which compose rather than
commute. Commute then swaps the remaining same-site pairs (Fock, Spin, PhaseSpace and
CollectiveTransition) and emits their commutator residuals. The trailing reduce catches a
composition that a residual creates next to a same-site neighbour, such as the contracted
operator of a Spin commutator.
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
    c = ca * cb
    if needs_diag_split
        accumulate_with_diag!(out, ops, c, sum_indices, ne)
    else
        canonicalize!(out, ops, c, ne)
    end
    return out
end

function push_op_indices!(inds::Vector{Index}, ops::Vector{Op})
    for op in ops
        has_index(op.index) && op.index ∉ inds && push!(inds, op.index)
    end
    return inds
end

distinct_op_indices(ops::Vector{Op}) = push_op_indices!(Index[], ops)

function has_two_indices(ops::Vector{Op})
    first_idx = NO_INDEX
    for op in ops
        has_index(op.index) || continue
        first_idx === NO_INDEX && (first_idx = op.index; continue)
        op.index == first_idx || return true
    end
    return false
end

const NO_DIAGONALS = Tuple{Index, Index}[]

function binder_diagonals(
        ops::Vector{Op}, c::CNum, s::Index, bound::Vector{Index}, ne::Vector{NonEqualPair},
    )
    pairs = NO_DIAGONALS
    depends_on_index_term(c, ops, s) || return pairs
    for e in distinct_op_indices(ops)
        may_coincide(s, e, ne) || continue
        e in bound && name_rank(s.name_id) <= name_rank(e.name_id) && continue
        pairs === NO_DIAGONALS && (pairs = Tuple{Index, Index}[])
        push!(pairs, (s, e))
    end
    return pairs
end

merge_ne_pairs(ne::Vector{NonEqualPair}, pairs::Vector{NonEqualPair}) = merge_ne(ne, pairs)
merge_ne_pairs(ne::Vector{NonEqualPair}, i::Index, others::Vector{Index}) =
    isempty(others) ? ne : merge_ne(ne, NonEqualPair[(i, j) for j in others])

"""
    partition_binders!(leaf, sides, c, ne, scope, bound) -> nothing

Split the sum over `scope` by the coincidences of its binders, one binder at a time: the
off-diagonal part keeps the binder distinct from every index it may equal, each diagonal
substitutes it, and the next binder is split inside every part. Splitting all binders at
once would count a joint coincidence such as `i = k, j = l` once per binder.

`sides` are operator strings that share the substitution (one for an operator, two for a
sandwich). Each part reaches `leaf(sides, c, ne, scope, absorbed)`, with `scope` the binders
left after substitution and `absorbed` the bound indices a binder was substituted onto.
"""
function partition_binders!(
        leaf::F, sides::NTuple{M, Vector{Op}}, c::CNum, ne::Vector{NonEqualPair},
        scope::Vector{Index}, bound::Vector{Index},
        k::Int = 1, absorbed::Vector{Index} = EMPTY_INDICES,
    ) where {F, M}
    k > length(scope) && return leaf(sides, c, ne, scope, absorbed)
    s = scope[k]
    ops = M == 1 ? sides[1] : reduce(vcat, sides)
    pairs = binder_diagonals(ops, c, s, bound, ne)
    partition_binders!(leaf, sides, c, merge_ne_pairs(ne, pairs), scope, bound, k + 1, absorbed)
    isempty(pairs) && return nothing
    rest = Index[x for x in scope if x != s]
    for (_, e) in pairs
        sub = map(side -> Op[change_index(o, s, e) for o in side], sides)
        sub_absorbed = e in bound ? vcat(absorbed, Index[e]) : absorbed
        partition_binders!(
            leaf, sub, change_index(c, s, e), diagonal_ne(ne, pairs, s, e), rest, bound, k,
            sub_absorbed,
        )
    end
    return nothing
end

function may_coincide(sum_idx::Index, ext_idx::Index, ne::Vector{NonEqualPair})
    ext_idx == sum_idx && return false
    ext_idx.space_index == sum_idx.space_index || return false
    return !ne_contains(ne, sum_idx, ext_idx)
end

coinciding_pairs(i::Index, ops::Vector{Op}, ne::Vector{NonEqualPair}) =
    Tuple{Index, Index}[(i, ext) for ext in distinct_op_indices(ops) if may_coincide(i, ext, ne)]

function sibling_exclusions(diag_pairs::Vector{Tuple{Index, Index}}, sum_idx::Index, ext_idx::Index)
    sibling_ne = NonEqualPair[]
    for (other_sum, other_ext) in diag_pairs
        other_sum == sum_idx || continue
        other_ext == ext_idx && break
        other_ext.space_index == ext_idx.space_index || continue
        push!(sibling_ne, (ext_idx, other_ext))
    end
    return sibling_ne
end

function diagonal_ne(
        ne::Vector{NonEqualPair}, diag_pairs::Vector{Tuple{Index, Index}}, sum_idx::Index,
        ext_idx::Index,
    )
    sub_ne = substitute_ne(ne, sum_idx => ext_idx)
    sibling_ne = sibling_exclusions(diag_pairs, sum_idx, ext_idx)
    return isempty(sibling_ne) ? sub_ne : merge_ne(sub_ne, sibling_ne)
end

"""
    accumulate_with_diag!(out, ops, c, sum_indices, ne) -> out

When summing over indices that may coincide with another operator's index, emit the
off-diagonal contribution (under `ne ∪ {(sum_idx, ext_idx)}` enforcing the indices differ)
and each diagonal contribution (substituting `sum_idx -> ext_idx`, with each constraint on
`sum_idx` rewritten onto `ext_idx`), partitioned binder by binder with
[`partition_binders!`](@ref).
"""
function accumulate_with_diag!(
        out::QTermDict, ops::Vector{Op}, c::CNum,
        sum_indices::Vector{Index}, ne::Vector{NonEqualPair},
    )
    if !has_two_indices(ops)
        canonicalize!(out, copy(ops), c, ne)
        return out
    end
    partition_binders!((ops,), c, ne, sum_indices, sum_indices) do (part,), pc, pne, _, absorbed
        if isempty(absorbed)
            canonicalize!(out, copy(part), pc, pne)
        else
            emit_scaled_by_scope!(out, copy(part), pc, pne, absorbed)
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
`eliminate_binders`: the term is multiplied by the index's constrained count
(`range - |partners|`), or keeps the binder and its constraints when that count is unknown.
"""
function settle_scope!(out::QTermDict, temp::QTermDict, scope::Vector{Index})
    for (term, cv) in temp
        mentions(c, idx) = depends_on_index_term(c, term.ops, idx)
        coef, ne = eliminate_binders(mentions, cv, term.ne, scope)
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
