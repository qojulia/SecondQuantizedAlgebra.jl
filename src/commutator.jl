"""
    commutator(a, b)

Compute the commutator `[a, b] = a*b - b*a` and simplify the result.
Always returns `QAdd{Complex{Int}}`. Short-circuits to zero when operands
commute (different Hilbert spaces, identical operators, or scalar arguments).
"""
function commutator end

const _ZERO_QADD = QAdd(QMul{Complex{Int}}[QMul(zero(Complex{Int}), QSym[])])

# Scalars commute with everything
commutator(::Number, ::Number) = _ZERO_QADD
commutator(::Number, ::QField) = _ZERO_QADD
commutator(::QField, ::Number) = _ZERO_QADD

# QSym, QSym: short-circuit when on different spaces or identical
function commutator(a::QSym, b::QSym)
    a.space_index == b.space_index || return _ZERO_QADD
    isequal(a, b) && return _ZERO_QADD
    return simplify(a * b - b * a)
end

# QMul, QSym: short-circuit when no shared space
function commutator(a::QMul, b::QSym)
    _has_space(a.args_nc, b.space_index) || return _ZERO_QADD
    return simplify(a * b - b * a)
end
function commutator(a::QSym, b::QMul)
    _has_space(b.args_nc, a.space_index) || return _ZERO_QADD
    return simplify(a * b - b * a)
end

# QMul, QMul: short-circuit when no shared spaces
function commutator(a::QMul, b::QMul)
    _has_shared_space(a.args_nc, b.args_nc) || return _ZERO_QADD
    return simplify(a * b - b * a)
end

# Zero-allocation space_index checks
function _has_space(ops::Vector{QSym}, si::Int)
    for op in ops
        op.space_index == si && return true
    end
    return false
end

function _has_shared_space(a::Vector{QSym}, b::Vector{QSym})
    for x in a, y in b
        x.space_index == y.space_index && return true
    end
    return false
end

# QAdd: distribute (bilinearity)
# Collect all QMul terms into a flat vector — avoids O(n²) intermediate QAdd copies.
function commutator(a::QAdd, b::QAdd)
    all_terms = QMul{Complex{Int}}[]
    for a_ in a.arguments, b_ in b.arguments
        _append_terms!(all_terms, commutator(a_, b_))
    end
    return QAdd(all_terms)
end
function commutator(a::QAdd, b::Union{QSym, QMul})
    all_terms = QMul{Complex{Int}}[]
    for a_ in a.arguments
        _append_terms!(all_terms, commutator(a_, b))
    end
    return QAdd(all_terms)
end
function commutator(a::Union{QSym, QMul}, b::QAdd)
    all_terms = QMul{Complex{Int}}[]
    for b_ in b.arguments
        _append_terms!(all_terms, commutator(a, b_))
    end
    return QAdd(all_terms)
end

function _append_terms!(all_terms::Vector{QMul{Complex{Int}}}, result::QAdd{Complex{Int}})
    for t in result.arguments
        push!(all_terms, t)
    end
    return all_terms
end
