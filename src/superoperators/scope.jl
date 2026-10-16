function Σ(S::SuperOperator, i::Index, non_equal::Vector{Index} = Index[])
    require_sum_range(i)
    i in S.indices && shared_binder_error(i)
    if isempty(non_equal) && !mentions_index(S, i)
        return S * index_range(i)
    end
    out = SuperTermDict()
    scope = Index[i]
    for (t, c) in S.arguments
        term_scope = push!(bound_indices(t, c, S.indices), i)
        ne = merge_ne_pairs(t.ne, i, non_equal)
        partition_binders!((t.left, t.right), c, ne, scope, scope) do (l, r), pc, pne, pscope, _
            settle_sandwich!(out, l, r, pc, pne, isempty(pscope) ? pscope : term_scope)
        end
    end
    return superop_pruned(out, push!(copy(S.indices), i))
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

function completeness_side(ops::Vector{Op}, ne::Vector{NonEqualPair})
    out = QTermDict()
    expand_gs_ops(copy(ops), CNUM_ONE) do expanded, c
        stream!(out, expanded, c, ne)
    end
    return out
end

expand_completeness(S::SuperOperator) = map_sides(completeness_side, S)

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
