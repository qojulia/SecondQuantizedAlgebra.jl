"""
    partial_sort!(ops::Vector{Op}, ne::Vector{NonEqualPair})

In-place stable partial-sort using `site_compare`. Distinct-site adjacents are
placed in canonical order; `Equal` and `Undetermined` pairs are left in
physical order so the sibling passes can interpret them.
"""
function partial_sort!(ops::Vector{Op}, ne::Vector{NonEqualPair})
    n = length(ops)
    n < 2 && return ops
    for i in 2:n
        j = i
        while j > 1
            cmp = site_compare(ops[j - 1], ops[j], ne)
            cmp === Greater || break
            ops[j - 1], ops[j] = ops[j], ops[j - 1]
            j -= 1
        end
    end
    return ops
end

"""
    canonicalize_to_dict!(out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair})

Terminal sink for the pipeline. Constructs `QTerm(ops, ne)`, looks up in `out`,
sums `c`, drops zero-coefficient entries. Takes ownership of `ops`; the caller
must not mutate after this call. Pre-condition: `ops` is in canonical form.
"""
@inline function canonicalize_to_dict!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
    )
    iszero(c) && return out
    return addto_key!(out, QTerm(ops, canonical_ne(ne)), c)
end

"""
    reduce_ops(sink, ops, c)

Fold adjacent provably-same-site pairs via `reduce_pair`. Single output per
input. Mutates `ops` in place; sink receives the same Vector.

Pre: `ops` is in canonical partial-sort order. Post: no adjacent provably-
same-site pair `(a, b)` for which `reduce_pair(a, b)` returns a reduction.
"""
@inline function reduce_ops(sink::F, ops::Vector{Op}, c::CNum) where {F}
    n = length(ops)
    n < 2 && (sink(ops, c); return)
    i = 1
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        if !may_reduce(a, b)
            i += 1
            continue
        end
        kind, new_op, factor = reduce_pair(a, b)
        if kind === NoReduction
            i += 1
        else
            iszero(factor) && return
            c = c * factor
            if kind === ScalarReduction
                deleteat!(ops, i:(i + 1))
            else
                ops[i] = new_op
                deleteat!(ops, i + 1)
            end
            i = max(i - 1, 1)
        end
    end
    sink(ops, c)
    return
end

"""
    commute_ops(sink, ops, c)

Apply commutation swaps for adjacent same-site pairs whose `can_commute` is
false. May fork: emits the swapped branch through the sink, plus up to two
residual branches where the pair is replaced by its residual operator, or
removed when the residual is the identity.

Pre-condition: pairs that reduce locally (Transition, Pauli) should have been
folded by `reduce_ops` first. The first residual branch mutates `ops` in place;
the other branches receive a `copy(ops)`.
"""
@inline function commute_ops(sink::F, ops::Vector{Op}, c::CNum) where {F}
    commute_ops_at(sink, ops, c, 1)
    return
end

function commute_ops_at(sink::F, ops::Vector{Op}, c::CNum, start::Int) where {F}
    i = start
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        if site_compare(a, b, EMPTY_NE) === Equal && !can_commute(a, b)
            c1, r1, c2, r2, n = commute_pair(a, b)
            swapped = copy(ops)
            swapped[i], swapped[i + 1] = b, a
            commute_ops_at(sink, swapped, c, max(i - 1, 1))
            iszero(c2) || commute_ops_at(sink, splice_residual!(copy(ops), i, r2, n), c * c2, max(i - 1, 1))
            iszero(c1) && return
            splice_residual!(ops, i, r1, n)
            c = c * c1
            i = max(i - 1, 1)
        else
            i += 1
        end
    end
    sink(ops, c)
    return
end

@inline function splice_residual!(ops::Vector{Op}, i::Int, r::Op, n::Int)
    n == 0 && return deleteat!(ops, i:(i + 1))
    ops[i] = r
    return deleteat!(ops, i + 1)
end

"""
    expand_gs_ops(sink, ops, c)

Apply `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every `σᵍᵍ` in `ops`. May fork by `n_levels`
per `σᵍᵍ`; recurses to handle multiple `σᵍᵍ` atoms in one term. Each branch
receives its own `copy(ops)`; the original is not mutated.
"""
function find_gs_idx(ops::Vector{Op})
    for k in eachindex(ops)
        is_gs, _, _, _ = ground_state_expand(ops[k])
        is_gs && return k
    end
    return 0
end

function expand_gs_ops(sink::F, ops::Vector{Op}, c::CNum) where {F}
    idx = find_gs_idx(ops)
    idx == 0 && (sink(ops, c); return)

    op = ops[idx]
    _, g, n_levels, site = ground_state_expand(op)
    id_ops = Op[]
    sizehint!(id_ops, length(ops) - 1)
    for k in eachindex(ops)
        k == idx || push!(id_ops, ops[k])
    end
    expand_gs_ops(sink, id_ops, c)
    neg_c = -c
    for k in 1:n_levels
        k == g && continue
        new_ops = copy(ops)
        new_ops[idx] = Op(OP_TRANSITION, op.name_id, site, op.index, k, k, g, n_levels)
        expand_gs_ops(sink, new_ops, neg_c)
    end
    return
end
