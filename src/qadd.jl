"""
    QAdd <: QTerm

Sum of [`QMul`](@ref) terms, optionally with summation indices for symbolic sums.

Internally stores a `Dict{Vector{QSym}, CNum}` mapping operator sequences to prefactors.
Like terms are auto-collected on construction. Zero-prefactor terms are dropped.

Fields:
- `dict::QTermDict` — operator sequence → prefactor
- `indices::Vector{Index}` — summation indices (empty = regular sum)
- `non_equal::Vector{Tuple{Index,Index}}` — pairwise inequality constraints
"""
const QTermDict = Dict{Vector{QSym}, CNum}

# Add a prefactor to an operator key in a QTermDict.
# Drops the entry immediately if the sum is zero — avoids needing a full-dict scan.
function _addto!(d::QTermDict, key::Vector{QSym}, val::CNum)
    existing = get(d, key, nothing)
    if existing === nothing
        # New key — just insert, skip the expensive CNum zero-addition
        _iszero_cnum(val) || (d[key] = val)
    else
        new_val = existing + val
        if _iszero_cnum(new_val)
            delete!(d, key)
        else
            d[key] = new_val
        end
    end
    return d
end

# Drop zero-prefactor entries from a QTermDict.
# Only needed after bulk operations (e.g. scalar multiply) that may produce zeros.
# IMPORTANT: Do NOT use filter!(predicate, ::Dict) — it allocates heavily in Julia.
function _drop_zeros!(d::QTermDict)
    for (ops, c) in d
        _iszero_cnum(c) && delete!(d, ops)
    end
    return d
end

# Merge two index/non_equal vectors, avoiding allocation when one is empty
function _merge_unique(a::Vector{T}, b::Vector{T}) where {T}
    isempty(a) && return b
    isempty(b) && return a
    return unique!(vcat(a, b))
end

struct QAdd <: QTerm
    arguments::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
end

# Primary constructor from QMul vector — auto-collects like terms, drops zeros
function QAdd(
        terms::Vector{QMul}, indices::Vector{Index},
        non_equal::Vector{Tuple{Index, Index}}
    )
    d = QTermDict()
    for t in terms
        _addto!(d, t.args_nc, t.arg_c)
    end
    return QAdd(d, indices, non_equal)
end
QAdd(terms::Vector{QMul}) = QAdd(terms, Index[], Tuple{Index, Index}[])

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

# Equality and hashing — order-independent via Dict comparison
function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    a.indices == b.indices || return false
    a.non_equal == b.non_equal || return false
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
function Base.hash(q::QAdd, h::UInt)
    return hash(:QAdd, hash(q.arguments, hash(q.indices, hash(q.non_equal, h))))
end

# Adjoint
function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (ops, c) in q.arguments
        adj_ops = QSym[adjoint(op) for op in reverse(ops)]
        _site_sort!(adj_ops)
        _addto!(d, adj_ops, conj(c))
    end
    return QAdd(d, q.indices, q.non_equal)
end

# --- Iteration helpers ---

"""
    terms(q::QAdd)

Iterate over `QMul` views of each term (unordered). For internal computation only.
The returned QMul shares its `args_nc` vector with the Dict key — callers must not mutate it.
"""
terms(q::QAdd) = (QMul(c, ops) for (ops, c) in q.arguments)

"""
    sorted_terms(q::QAdd) -> Vector{QMul}

Return terms in deterministic order. For printing, TermInterface, and iteration.
Uses a richer sort key than `_sort_key` to distinguish operator types for display.
"""
function sorted_terms(q::QAdd)
    isempty(q.arguments) && return QMul[]
    pairs = sort!(collect(q.arguments); by = p -> _term_sort_key(p.first))
    return QMul[QMul(c, ops) for (ops, c) in pairs]
end

# Sort key for QAdd term ordering — richer than _sort_key to distinguish types.
# This does NOT affect _site_sort! (which preserves same-site operator order).
_term_sort_key(ops::Vector{QSym}) = (length(ops), map(_full_op_key, ops)...)
_full_op_key(op::QSym) = (_sort_key(op)..., _type_order(op), op.name)
_type_order(::Destroy) = 0
_type_order(::Create) = 1
_type_order(::Transition) = 2
_type_order(::Pauli) = 3
_type_order(::Spin) = 4
_type_order(::Position) = 5
_type_order(::Momentum) = 6

"""
    Base.getindex(q::QAdd, key::Vector{QSym}) -> CNum

Look up the prefactor for a given operator sequence. Returns zero if absent.
"""
Base.getindex(q::QAdd, key::Vector{QSym}) = get(q.arguments, key, _CNUM_ZERO)

# Iterator interface — yields QMul terms in Dict order (fast, no allocation per step).
# For deterministic order use sorted_terms() or SymbolicUtils.arguments().
function Base.iterate(q::QAdd)
    it = iterate(q.arguments)
    it === nothing && return nothing
    (ops, c), state = it
    return (QMul(c, ops), state)
end
function Base.iterate(q::QAdd, state)
    it = iterate(q.arguments, state)
    it === nothing && return nothing
    (ops, c), state = it
    return (QMul(c, ops), state)
end
Base.eltype(::Type{QAdd}) = QMul

# --- Helpers: wrap QSym/scalar as QMul ---

_to_qmul(a::QSym) = QMul(_to_cnum(1), QSym[a])
_to_qmul(a::QMul) = a
_scalar_qmul(x::Number) = QMul(_to_cnum(x), QSym[])

## Addition — always returns QAdd

# QSym + QSym
Base.:+(a::QSym, b::QSym) = QAdd(QMul[_to_qmul(a), _to_qmul(b)])

# QMul + QMul
Base.:+(a::QMul, b::QMul) = QAdd(QMul[a, b])

# QMul + QSym
Base.:+(a::QMul, b::QSym) = QAdd(QMul[a, _to_qmul(b)])
Base.:+(a::QSym, b::QMul) = b + a

# QAdd + QMul — merge into existing dict
function Base.:+(a::QAdd, b::QMul)
    d = copy(a.arguments)
    _addto!(d, b.args_nc, b.arg_c)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:+(a::QMul, b::QAdd) = b + a

# QAdd + QSym
Base.:+(a::QAdd, b::QSym) = a + _to_qmul(b)
Base.:+(a::QSym, b::QAdd) = b + a

# QAdd + QAdd — merge dicts
function Base.:+(a::QAdd, b::QAdd)
    d = copy(a.arguments)
    for (ops, c) in b.arguments
        _addto!(d, ops, c)
    end
    indices = _merge_unique(a.indices, b.indices)
    non_equal = _merge_unique(a.non_equal, b.non_equal)
    return QAdd(d, indices, non_equal)
end

# QField + Number
Base.:+(a::QSym, b::Number) = QAdd(QMul[_to_qmul(a), _scalar_qmul(b)])
Base.:+(a::Number, b::QSym) = b + a
Base.:+(a::QMul, b::Number) = QAdd(QMul[a, _scalar_qmul(b)])
Base.:+(a::Number, b::QMul) = b + a
Base.:+(a::QAdd, b::Number) = a + _scalar_qmul(b)
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QTermDict(ops => -c for (ops, c) in a.arguments), a.indices, a.non_equal)
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

## QAdd * ... (distributive)

# QAdd * Number
function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict(ops => c * cb for (ops, c) in a.arguments)
    _drop_zeros!(d)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym
function Base.:*(a::QAdd, b::QSym)
    result = QMul[QMul(c, ops) * b for (ops, c) in a.arguments]
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QSym, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.arguments]
    return QAdd(result, b.indices, b.non_equal)
end

# QAdd * QMul
function Base.:*(a::QAdd, b::QMul)
    result = QMul[QMul(c, ops) * b for (ops, c) in a.arguments]
    return QAdd(result, a.indices, a.non_equal)
end
function Base.:*(a::QMul, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.arguments]
    return QAdd(result, b.indices, b.non_equal)
end

# QAdd * QAdd
function Base.:*(a::QAdd, b::QAdd)
    result = QMul[]
    for (ops_a, c_a) in a.arguments, (ops_b, c_b) in b.arguments
        push!(result, QMul(c_a, ops_a) * QMul(c_b, ops_b))
    end
    indices = _merge_unique(a.indices, b.indices)
    non_equal = _merge_unique(a.non_equal, b.non_equal)
    return QAdd(result, indices, non_equal)
end

# QAdd / Number
Base.:/(a::QAdd, b::Number) = a * inv(b)


"""
    _depends_on_index(m::QMul, idx::Index) -> Bool

Check whether a `QMul` term depends on the given `Index`.
Checks both operators (via `.index`) and the c-number prefactor
(via `Symbolics.get_variables` for the index symbol).
"""
function _depends_on_index(m::QMul, idx::Index)
    for op in m.args_nc
        op.index == idx && return true
    end
    isym = SymbolicUtils.unwrap(idx.sym)
    c = m.arg_c
    for part in (real(c), imag(c))
        vars = Symbolics.get_variables(part)
        any(v -> isequal(v, isym), vars) && return true
    end
    return false
end

"""
    Σ(expr, i::Index, non_equal::Vector{Index}=Index[])

Create a symbolic sum over index `i`. Returns a `QAdd` with summation indices.
If the expression does not depend on `i`, returns `range * expr` instead.

    Σ(expr, i::Index, j::Index, ...)

Create a multi-index sum.
"""
function Σ(expr::QMul, i::Index, non_equal::Vector{Index} = Index[])
    if !_depends_on_index(expr, i)
        return QAdd(QMul[QMul(i.range * expr.arg_c, expr.args_nc)])
    end
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    return QAdd(QMul[expr], [i], ne_pairs)
end
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    all_indices = vcat(expr.indices, [i])
    all_ne = vcat(expr.non_equal, ne_pairs)
    return QAdd(expr.arguments, all_indices, all_ne)
end
function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_to_qmul(expr), i, non_equal)
end
function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_scalar_qmul(expr), i, non_equal)
end

# Multi-index: Σ(expr, i, j) = double sum
function Σ(expr, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ

"""
    expand_sums(expr::QAdd) -> QAdd

Explicit diagonal splitting for symbolic sums. Called by QC's `scale()`.

Core rule: `Σ_i(A_i * B_j)` where i,j same space, j not in non_equal
→ `Σ_{i≠j}(A_i * B_j) + A_j * B_j`
"""
function expand_sums(s::QAdd)
    isempty(s.indices) && return s  # Nothing to expand

    d = QTermDict()
    result_ne = copy(s.non_equal)

    for (ops, c) in s.arguments
        term = QMul(c, ops)
        term_indices = get_indices(term)

        for idx in term_indices
            for sum_idx in s.indices
                if idx != sum_idx &&
                        idx.space_index == sum_idx.space_index &&
                        !((sum_idx, idx) in s.non_equal) &&
                        !((idx, sum_idx) in s.non_equal)
                    diag_term = change_index(term, sum_idx, idx)
                    d[diag_term.args_nc] = diag_term.arg_c
                    push!(result_ne, (sum_idx, idx))
                end
            end
        end

        d[ops] = c
    end

    return QAdd(d, copy(s.indices), result_ne)
end

# Passthrough for non-sum types
expand_sums(m::QMul) = m
expand_sums(op::QSym) = op
