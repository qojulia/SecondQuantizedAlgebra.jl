"""
    expand_sums(expr::QAdd) -> QAdd

Explicit diagonal splitting for symbolic sums. Called by QC's `scale()`.

Core rule: `Σ_i(A_i * B_j)` where i,j same space, j not in non_equal
→ `Σ_{i≠j}(A_i * B_j) + A_j * B_j`
"""
function expand_sums(s::QAdd)
    isempty(s.indices) && return s  # Nothing to expand

    result_terms = QMul[]
    result_ne = copy(s.non_equal)

    for term in s.arguments
        term_indices = get_indices(term)
        needs_split = false

        for idx in term_indices
            for sum_idx in s.indices
                if idx != sum_idx &&
                        idx.space_index == sum_idx.space_index &&
                        !((sum_idx, idx) in s.non_equal) &&
                        !((idx, sum_idx) in s.non_equal)
                    # Diagonal split: generate term with sum_idx → idx
                    diag_term = change_index(term, sum_idx, idx)
                    push!(result_terms, diag_term)
                    push!(result_ne, (sum_idx, idx))
                    needs_split = true
                end
            end
        end

        push!(result_terms, term)
    end

    return QAdd(result_terms, copy(s.indices), result_ne)
end

# Passthrough for non-sum types
expand_sums(m::QMul) = m
expand_sums(op::QSym) = op
