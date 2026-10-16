function Σ(S::SuperOperator, i::Index, non_equal::Vector{Index} = Index[])
    require_sum_range(i)
    i in S.indices && shared_binder_error(i)
    if isempty(non_equal) && !mentions_index(S, i)
        return S * index_range(i)
    end
    off_diag = SuperTermDict()
    diag = SuperTermDict()
    scope_i = Index[i]
    for (t, c) in S.arguments
        ne_aug = merge_ne_pairs(t.ne, i, non_equal)
        if !depends_on_sandwich(c, t.left, t.right, i)
            settle_sandwich!(off_diag, t.left, t.right, c, ne_aug, scope_i)
            continue
        end
        diag_pairs = coinciding_pairs(i, vcat(t.left, t.right), ne_aug)
        settle_sandwich!(off_diag, t.left, t.right, c, merge_ne_pairs(ne_aug, diag_pairs), scope_i)
        for (_, ext) in diag_pairs
            canonical_sandwich!(
                diag, Op[change_index(o, i, ext) for o in t.left],
                Op[change_index(o, i, ext) for o in t.right], change_index(c, i, ext),
                diagonal_ne(ne_aug, diag_pairs, i, ext),
            )
        end
    end
    binds = any(((t, c),) -> depends_on_scope_sandwich(c, t.left, t.right, t.ne, i), off_diag)
    result = superop_pruned(off_diag, binds ? vcat(S.indices, scope_i) : copy(S.indices))
    isempty(diag) && return result
    return result + superop_pruned(diag, copy(S.indices))
end

Σ(S::SuperOperator, i::Index, j::Index, rest::Index...) = Σ(Σ(S, i), j, rest...)

function change_index(S::SuperOperator, m::IndexMap)
    renames_nothing(m) && return S
    out = SuperTermDict()
    new_indices = Index[reindex(m, idx) for idx in S.indices]
    for (t, c) in S.arguments
        ne_becomes_contradictory(t.ne, m) && continue
        bound = Index[reindex(m, idx) for idx in bound_indices(t, c, S.indices)]
        emit_sandwich!(
            out, Op[change_index(o, m) for o in t.left],
            Op[change_index(o, m) for o in t.right], change_index(c, m),
            substitute_ne(t.ne, m), bound,
        )
    end
    return superop_pruned(out, unique!(new_indices))
end

function expand_completeness(S::SuperOperator)
    out = SuperTermDict()
    for (t, c) in S.arguments
        lt = QTermDict()
        expand_gs_ops(copy(t.left), CNUM_ONE) do ops1, c1
            stream!(lt, ops1, c1, t.ne)
        end
        rt = QTermDict()
        expand_gs_ops(copy(t.right), CNUM_ONE) do ops1, c1
            stream!(rt, ops1, c1, t.ne)
        end
        temp = SuperTermDict()
        for (a, ca) in lt, (b, cb) in rt
            addto_key!(temp, SuperTerm(a.ops, b.ops, t.ne), c * (ca * cb))
        end
        settle_scope!(out, temp, bound_indices(t, c, S.indices))
    end
    return superop_pruned(out, copy(S.indices))
end

function assume_distinct_index(S::SuperOperator, pairs::Vector{Tuple{Index, Index}})
    out = SuperTermDict()
    for (t, c) in S.arguments
        ne = merge_ne_pairs(t.ne, pairs)
        settle_sandwich!(out, t.left, t.right, c, ne, bound_indices(t, c, S.indices))
    end
    return expand_completeness(superop_pruned(out, copy(S.indices)))
end

function get_indices(S::SuperOperator)
    inds = copy(S.indices)
    for t in keys(S.arguments)
        push_op_indices!(inds, t.left)
        push_op_indices!(inds, t.right)
    end
    return inds
end
