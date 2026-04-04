"""
    normal_order(expr)
    normal_order(expr, spaces)

Alias for `simplify(expr, NormalOrder())`. Apply algebra rules and rewrite
the expression in normal order (creation operators left of annihilation).
The `spaces` argument enables ground state rewriting for Transitions.
Always returns `QAdd`.
"""
normal_order(expr::QField) = simplify(expr, NormalOrder())
normal_order(expr::QField, h::HilbertSpace) = simplify(expr, NormalOrder(), h)

# Ground state rewriting: |g⟩⟨g| = 1 - Σ_{k≠g} |k⟩⟨k|

function _apply_ground_state(expr::QAdd{CT}, h::NLevelSpace) where {CT}
    all_terms = QMul{CT}[]
    for term in expr.arguments
        expanded = _expand_ground_state(term, h)
        append!(all_terms, expanded.arguments)
    end
    return QAdd(all_terms)
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

# Expand a single QMul term, recursively until no ground state projections remain
function _expand_ground_state(term::QMul{CT}, h::NLevelSpace) where {CT}
    g = h.ground_state
    n = h.n
    for (idx, op) in enumerate(term.args_nc)
        if op isa Transition && op.i == g && op.j == g
            terms = QMul{CT}[]
            # Identity term (remove the transition)
            id_ops = QSym[term.args_nc[1:(idx - 1)]..., term.args_nc[(idx + 1):end]...]
            push!(terms, QMul(term.arg_c, id_ops))
            # Subtraction terms
            for k in 1:n
                k == g && continue
                new_op = Transition(op.name, k, k, op.space_index)
                sub_ops = QSym[term.args_nc[1:(idx - 1)]..., new_op, term.args_nc[(idx + 1):end]...]
                push!(terms, QMul(-term.arg_c, sub_ops))
            end
            # Recurse: expanded terms may still contain ground state projections
            result = QMul{CT}[]
            for t in terms
                expanded = _expand_ground_state(t, h)
                append!(result, expanded.arguments)
            end
            return QAdd(result)
        end
    end
    return QAdd(QMul{CT}[term])
end
