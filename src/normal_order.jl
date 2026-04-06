"""
    normal_order(expr)
    normal_order(expr, spaces)

Alias for `simplify(expr, NormalOrder())`. Apply algebra rules and rewrite
the expression in normal order (creation operators left of annihilation).
The `spaces` argument enables ground state rewriting for Transitions.
Always returns `QAdd`.
"""
normal_order(expr::QField) = _qsimplify(expr, NormalOrder())
normal_order(expr::QField, h::HilbertSpace) = _qsimplify(expr, NormalOrder(), h)

# Ground state rewriting: |g⟩⟨g| = 1 - Σ_{k≠g} |k⟩⟨k|

function _apply_ground_state(expr::QAdd, h::NLevelSpace)
    d = QTermDict()
    for term in terms(expr)
        expanded = _expand_ground_state(term, h)
        for (ops, c) in expanded.arguments
            _addto!(d, ops, c)
        end
    end
    return QAdd(d, expr.indices, expr.non_equal)
end

# ProductSpace: apply ground state rewriting for each NLevelSpace sub-space
function _apply_ground_state(expr::QAdd, h::ProductSpace)
    result = expr
    for space in h.spaces
        result = _apply_ground_state(result, space)
    end
    return result
end

# No-op for non-NLevel spaces
_apply_ground_state(expr::QAdd, ::FockSpace) = expr
_apply_ground_state(expr::QAdd, ::PauliSpace) = expr
_apply_ground_state(expr::QAdd, ::SpinSpace) = expr
_apply_ground_state(expr::QAdd, ::PhaseSpace) = expr
_apply_ground_state(expr::QAdd, ::ClusterSpace) = expr

# Expand a single QMul term, recursively until no ground state projections remain
function _expand_ground_state(term::QMul, h::NLevelSpace)
    g = h.ground_state
    n = h.n
    for (idx, op) in enumerate(term.args_nc)
        if op isa Transition && op.i == g && op.j == g
            expanded_terms = QMul[]
            # Identity term (remove the transition)
            id_ops = QSym[term.args_nc[1:(idx - 1)]..., term.args_nc[(idx + 1):end]...]
            push!(expanded_terms, QMul(term.arg_c, id_ops))
            # Subtraction terms
            for k in 1:n
                k == g && continue
                new_op = Transition(op.name, k, k, op.space_index, op.copy_index, op.index)
                sub_ops = QSym[term.args_nc[1:(idx - 1)]..., new_op, term.args_nc[(idx + 1):end]...]
                push!(expanded_terms, QMul(-term.arg_c, sub_ops))
            end
            # Recurse: expanded terms may still contain ground state projections
            d = QTermDict()
            for t in expanded_terms
                sub = _expand_ground_state(t, h)
                for (ops, c) in sub.arguments
                    _addto!(d, ops, c)
                end
            end
            return QAdd(d, Index[], Tuple{Index, Index}[])
        end
    end
    d = QTermDict(term.args_nc => term.arg_c)
    return QAdd(d, Index[], Tuple{Index, Index}[])
end
