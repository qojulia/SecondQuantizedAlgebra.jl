"""
    commutator(a, b)

Compute the commutator `[a, b] = a*b - b*a`.
With eager ordering, the result is already in canonical form.
Always returns `QAdd`. Short-circuits to zero when operands commute
(different sites, identical operators, or scalar arguments).
"""
function commutator end

const _ZERO_QADD = QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
_zero_qadd() = _ZERO_QADD

# Scalars commute with everything
commutator(::Number, ::Number) = _zero_qadd()
commutator(::Number, ::QField) = _zero_qadd()
commutator(::QField, ::Number) = _zero_qadd()

# QSym, QSym: short-circuit when on different sites or identical
function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _zero_qadd()
    isequal(a, b) && return _zero_qadd()
    return a * b - b * a
end

# QAdd, QSym
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
    return a * b - b * a
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
    return a * b - b * a
end

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return _zero_qadd()
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            if idx in b.indices
                d = QTermDict()
                for (ops, c) in b.arguments
                    single = QAdd(QTermDict(ops => c), Index[], Tuple{Index, Index}[])
                    _merge_terms!(d, commutator(a, single))
                end
                indices = _merge_unique(a.indices, b.indices)
                non_equal = _merge_unique(a.non_equal, b.non_equal)
                return QAdd(d, indices, non_equal)
            end
        end
    end
    return a * b - b * a
end

function _merge_terms!(d::QTermDict, result::QAdd)
    for (ops, c) in result.arguments
        _addto!(d, ops, c)
    end
    return d
end
