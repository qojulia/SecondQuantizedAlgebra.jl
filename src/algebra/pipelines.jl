"""
    _stream!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

The default canonicalization pipeline. Order: `reduce -> commute -> reduce -> sink`.

The first reduce folds Transition/Pauli same-site pairs (these never commute,
they compose). Commute then operates only on Fock/Spin/PhaseSpace ladder pairs.
The trailing reduce catches any same-site composition surfaced by the commute
residual (e.g. a Spin commutator's contracted op meeting a same-site neighbor).
"""
@inline function _stream!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _reduce_ops(ops, c) do ops1, c1
        _commute_ops(ops1, c1) do ops2, c2
            _reduce_ops(ops2, c2) do ops3, c3
                _canonicalize_to_dict!(out, ops3, c3, ne)
            end
        end
    end
    return nothing
end

"""
    _canonicalize!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

Re-establish canonical form: sort, then run the pipeline.
"""
@inline function _canonicalize!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _partial_sort!(ops, ne)
    _stream!(out, ops, c, ne)
    return nothing
end

"""
    _emit_product!(out, ta_ops, ca, ta_ne, tb_ops, cb, tb_ne, sum_indices, needs_diag)

Concatenate the operator strings of two terms, merge their constraint sets, and
route through `_canonicalize!` (or `_accumulate_with_diag!` when summation
indices may collapse).
"""
@inline function _emit_product!(
        out::QTermDict,
        ta_ops::Vector{QSym}, ca::CNum, ta_ne::Vector{NonEqualPair},
        tb_ops::Vector{QSym}, cb::CNum, tb_ne::Vector{NonEqualPair},
        sum_indices::Vector{Index}, needs_diag_split::Bool,
    )
    n = length(ta_ops) + length(tb_ops)
    ops = Vector{QSym}(undef, n)
    copyto!(ops, 1, ta_ops, 1, length(ta_ops))
    copyto!(ops, length(ta_ops) + 1, tb_ops, 1, length(tb_ops))
    ne = _merge_ne(ta_ne, tb_ne)
    c = _mul_cnum(ca, cb)
    if needs_diag_split
        _accumulate_with_diag!(out, ops, c, sum_indices, ne)
    else
        _canonicalize!(out, ops, c, ne)
    end
    return nothing
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

function _depends_on_index_ops(c::CNum, ops::Vector{QSym}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    return _depends_on_index_term(c, ops, idx)
end

"""
    _accumulate_with_diag!(out, ops, c, sum_indices, ne) -> nothing

When summing over indices that may coincide with another operator's index, emit
both the off-diagonal contribution (under `ne ∪ {(sum_idx, ext_idx)}` enforcing
the indices differ) and each diagonal contribution (substituting
`sum_idx -> ext_idx`, under `ne` with any constraint on `sum_idx` dropped).
"""
function _accumulate_with_diag!(
        out::QTermDict, ops::Vector{QSym}, c::CNum,
        sum_indices::Vector{Index}, ne::Vector{NonEqualPair},
    )
    distinct = _distinct_op_indices(ops)
    if length(distinct) < 2
        _canonicalize!(out, copy(ops), c, ne)
        return nothing
    end

    # Collect every (sum_idx, ext_idx) pair that needs a diagonal contribution.
    diag_pairs = Tuple{Index, Index}[]
    for sum_idx in sum_indices
        _depends_on_index_ops(c, ops, sum_idx) || continue
        for ext_idx in distinct
            ext_idx == sum_idx && continue
            ext_idx.space_index == sum_idx.space_index || continue
            _ne_contains(ne, sum_idx, ext_idx) && continue
            push!(diag_pairs, (sum_idx, ext_idx))
        end
    end

    if isempty(diag_pairs)
        _canonicalize!(out, copy(ops), c, ne)
        return nothing
    end

    # Off-diagonal: canonicalize under ne augmented with every (sum, ext) pair,
    # so that partial_sort can use those inequalities when ordering operators.
    augmented_ne = ne
    for (a, b) in diag_pairs
        augmented_ne = _merge_ne_pair(augmented_ne, a, b)
    end
    _canonicalize!(out, copy(ops), c, augmented_ne)

    # Diagonal contributions: one per (sum_idx, ext_idx).
    for (sum_idx, ext_idx) in diag_pairs
        sub_ops = QSym[change_index(o, sum_idx, ext_idx) for o in ops]
        sub_c = change_index(c, sum_idx, ext_idx)
        sub_ne = _drop_ne_with(ne, sum_idx)
        _canonicalize!(out, sub_ops, sub_c, sub_ne)
    end
    return nothing
end

function Base.:*(a::QSym, b::QSym)
    out = QTermDict()
    _emit_product!(
        out,
        QSym[a], _CNUM_ONE, _EMPTY_NE,
        QSym[b], _CNUM_ONE, _EMPTY_NE,
        _EMPTY_INDICES, false
    )
    return QAdd(out, _EMPTY_INDICES)
end

function Base.:*(a::QAdd, b::QSym)
    out = QTermDict()
    needs = !isempty(a.indices)
    for (ta, ca) in a
        _emit_product!(
            out, ta.ops, ca, ta.ne,
            QSym[b], _CNUM_ONE, _EMPTY_NE,
            a.indices, needs
        )
    end
    return QAdd(out, _absorb_pinned_sums(a.indices, a, b))
end

function Base.:*(a::QSym, b::QAdd)
    out = QTermDict()
    needs = !isempty(b.indices)
    for (tb, cb) in b
        _emit_product!(
            out, QSym[a], _CNUM_ONE, _EMPTY_NE,
            tb.ops, cb, tb.ne,
            b.indices, needs
        )
    end
    return QAdd(out, _absorb_pinned_sums(b.indices, a, b))
end

function Base.:*(a::QAdd, b::QAdd)
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            idx in b.indices && throw(
                ArgumentError(
                    "Summation index $(idx.name) appears in both factors. " *
                        "Bound variables on the two sides of a product must be " *
                        "distinct so the resulting double sum is unambiguous. " *
                        "Use `change_index` to rename one side, or construct the " *
                        "two sums with different indices."
                )
            )
        end
    end
    out = QTermDict()
    sum_indices = _merge_unique(a.indices, b.indices)
    needs = !isempty(sum_indices)
    for (ta, ca) in a, (tb, cb) in b
        _emit_product!(
            out, ta.ops, ca, ta.ne,
            tb.ops, cb, tb.ne,
            sum_indices, needs
        )
    end
    return QAdd(out, _absorb_pinned_sums(sum_indices, a, b))
end

# Drop summation indices pinned by the product: a bound `.indices` entry
# of one factor that every term of the other factor carries as a free op
# index is no longer ranging, so the sum scope must not survive. The
# uniformity check (`_every_term_has_op_index`) is what guarantees the
# result still has consistent `.indices` semantics across every term.
function _absorb_pinned_sums(
        sum_indices::Vector{Index}, a::Union{QAdd, QSym}, b::Union{QAdd, QSym}
    )
    isempty(sum_indices) && return sum_indices
    a_indices = a isa QAdd ? a.indices : _EMPTY_INDICES
    b_indices = b isa QAdd ? b.indices : _EMPTY_INDICES
    keep = Index[]
    @inbounds for s in sum_indices
        from_a = s in a_indices
        from_b = s in b_indices
        pinned = (from_a && !from_b && _every_term_has_op_index(b, s)) ||
            (from_b && !from_a && _every_term_has_op_index(a, s))
        pinned || push!(keep, s)
    end
    return length(keep) == length(sum_indices) ? sum_indices : keep
end

_every_term_has_op_index(q::QSym, idx::Index) = q.index == idx
function _every_term_has_op_index(q::QAdd, idx::Index)
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
