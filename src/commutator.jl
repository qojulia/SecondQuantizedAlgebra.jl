"""
    commutator(a, b)

Compute the commutator `[a, b] = a*b - b*a` and simplify the result.
Always returns `QAdd`. Short-circuits to zero when operands commute
(different sites, identical operators, or scalar arguments).
"""
function commutator end

const _ZERO_QADD = QAdd(QMul[QMul(_to_cnum(0), QSym[])])

# Scalars commute with everything
commutator(::Number, ::Number) = _ZERO_QADD
commutator(::Number, ::QField) = _ZERO_QADD
commutator(::QField, ::Number) = _ZERO_QADD

# QSym, QSym: short-circuit when on different sites or identical
function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _ZERO_QADD
    isequal(a, b) && return _ZERO_QADD
    return _qsimplify(a * b - b * a, NormalOrder())
end

# QMul, QSym: short-circuit when no shared site
function commutator(a::QMul, b::QSym)
    _has_site(a.args_nc, b) || return _ZERO_QADD
    return _qsimplify(a * b - b * a, NormalOrder())
end
function commutator(a::QSym, b::QMul)
    _has_site(b.args_nc, a) || return _ZERO_QADD
    return _qsimplify(a * b - b * a, NormalOrder())
end

# QMul, QMul: short-circuit when no shared sites
function commutator(a::QMul, b::QMul)
    _has_shared_site(a.args_nc, b.args_nc) || return _ZERO_QADD
    return _qsimplify(a * b - b * a, NormalOrder())
end

# Zero-allocation site checks
function _has_site(ops::Vector{QSym}, target::QSym)
    for op in ops
        _same_site(op, target) && return true
    end
    return false
end

function _has_shared_site(a::Vector{QSym}, b::Vector{QSym})
    for x in a, y in b
        _same_site(x, y) && return true
    end
    return false
end

# QAdd: distribute (bilinearity)
# For QAdd with indices (sums), handle diagonal collapse
function commutator(a::QAdd, b::QSym)
    # Diagonal collapse: if a is a sum and b is indexed on same space
    if !isempty(a.indices) && has_index(b.index)
        for (k, idx) in enumerate(a.indices)
            if idx.space_index == b.index.space_index
                collapsed = change_index(a, idx, b.index)
                new_indices = [a.indices[j] for j in eachindex(a.indices) if j != k]
                new_ne = filter(p -> p[1] != idx && p[2] != idx, a.non_equal)
                collapsed_qadd = QAdd(collapsed.arguments, new_indices, new_ne)
                return commutator(collapsed_qadd, b)
            end
        end
    end
    # Regular distribution
    all_terms = QMul[]
    for a_ in a.arguments
        _append_terms!(all_terms, commutator(a_, b))
    end
    return QAdd(all_terms, a.indices, a.non_equal)
end

function commutator(a::QSym, b::QAdd)
    if !isempty(b.indices) && has_index(a.index)
        for (k, idx) in enumerate(b.indices)
            if idx.space_index == a.index.space_index
                collapsed = change_index(b, idx, a.index)
                new_indices = [b.indices[j] for j in eachindex(b.indices) if j != k]
                new_ne = filter(p -> p[1] != idx && p[2] != idx, b.non_equal)
                collapsed_qadd = QAdd(collapsed.arguments, new_indices, new_ne)
                return commutator(a, collapsed_qadd)
            end
        end
    end
    all_terms = QMul[]
    for b_ in b.arguments
        _append_terms!(all_terms, commutator(a, b_))
    end
    return QAdd(all_terms, b.indices, b.non_equal)
end

function commutator(a::QAdd, b::QAdd)
    # Handle sum-sum: collapse inner first
    if !isempty(a.indices) && !isempty(b.indices)
        # Distribute outer over inner
        all_terms = QMul[]
        for b_ in b.arguments
            _append_terms!(all_terms, commutator(a, b_))
        end
        return QAdd(all_terms, b.indices, b.non_equal)
    end
    # Regular distribution
    all_terms = QMul[]
    for a_ in a.arguments, b_ in b.arguments
        _append_terms!(all_terms, commutator(a_, b_))
    end
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(all_terms, indices, non_equal)
end

function commutator(a::QAdd, b::QMul)
    # Diagonal collapse: if a is a sum and b contains an operator indexed on same space
    if !isempty(a.indices)
        for op in b.args_nc
            if has_index(op.index)
                for (k, idx) in enumerate(a.indices)
                    if idx.space_index == op.index.space_index
                        collapsed = change_index(a, idx, op.index)
                        new_indices = [a.indices[j] for j in eachindex(a.indices) if j != k]
                        new_ne = filter(p -> p[1] != idx && p[2] != idx, a.non_equal)
                        collapsed_qadd = QAdd(collapsed.arguments, new_indices, new_ne)
                        return commutator(collapsed_qadd, b)
                    end
                end
            end
        end
    end
    all_terms = QMul[]
    for a_ in a.arguments
        _append_terms!(all_terms, commutator(a_, b))
    end
    return QAdd(all_terms, a.indices, a.non_equal)
end
function commutator(a::QMul, b::QAdd)
    # Diagonal collapse: if b is a sum and a contains an operator indexed on same space
    if !isempty(b.indices)
        for op in a.args_nc
            if has_index(op.index)
                for (k, idx) in enumerate(b.indices)
                    if idx.space_index == op.index.space_index
                        collapsed = change_index(b, idx, op.index)
                        new_indices = [b.indices[j] for j in eachindex(b.indices) if j != k]
                        new_ne = filter(p -> p[1] != idx && p[2] != idx, b.non_equal)
                        collapsed_qadd = QAdd(collapsed.arguments, new_indices, new_ne)
                        return commutator(a, collapsed_qadd)
                    end
                end
            end
        end
    end
    all_terms = QMul[]
    for b_ in b.arguments
        _append_terms!(all_terms, commutator(a, b_))
    end
    return QAdd(all_terms, b.indices, b.non_equal)
end

function _append_terms!(all_terms::Vector{QMul}, result::QAdd)
    for t in result.arguments
        push!(all_terms, t)
    end
    return all_terms
end
