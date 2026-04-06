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

## Index helpers (must precede multiplication methods that use them)

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

function _any_depends_on_index(s::QAdd, idx::Index)
    for (ops, c) in s.arguments
        _depends_on_index(QMul(c, ops), idx) && return true
    end
    return false
end

"""
    _diagonal_split(terms, sum_idx::Index, non_equal::Vector{Tuple{Index,Index}})

Extract diagonal terms for summation index `sum_idx`. For each term containing a free
index `j` on the same space as `sum_idx` (and not already constrained `sum_idx != j`),
create the diagonal term via `change_index(term, sum_idx, j)` and record the constraint.

Returns `(diag_terms::Vector{QMul}, new_non_equal::Vector{Tuple{Index,Index}})`.
"""
function _diagonal_split(
        terms, sum_idx::Index,
        non_equal::Vector{Tuple{Index, Index}}
    )
    diag_terms = QMul[]
    new_ne = copy(non_equal)
    for term in terms
        term_indices = get_indices(term)
        for idx in term_indices
            if idx != sum_idx &&
                    idx.space_index == sum_idx.space_index &&
                    !((sum_idx, idx) in new_ne) &&
                    !((idx, sum_idx) in new_ne)
                diag_term = change_index(term, sum_idx, idx)
                push!(diag_terms, diag_term)
                push!(new_ne, (sum_idx, idx))
            end
        end
    end
    return diag_terms, new_ne
end

"""
    _apply_diagonal_split(s::QAdd, ext_idx::Index) -> QAdd

For a QAdd with summation indices, split out diagonal terms for any summation index
on the same space as `ext_idx`. Called after multiplication introduces `ext_idx`
into the terms.
"""
function _apply_diagonal_split(s::QAdd, ext_idx::Index)
    result = s
    for sum_idx in s.indices
        if sum_idx != ext_idx &&
                sum_idx.space_index == ext_idx.space_index &&
                !((sum_idx, ext_idx) in result.non_equal) &&
                !((ext_idx, sum_idx) in result.non_equal)
            diag_terms, new_ne = _diagonal_split(terms(result), sum_idx, result.non_equal)
            result = QAdd(result.arguments, result.indices, new_ne)
            if !isempty(diag_terms)
                result = result + QAdd(diag_terms)
            end
        end
    end
    return result
end

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
    product = QAdd(result, copy(a.indices), copy(a.non_equal))
    if has_index(b.index) && !isempty(a.indices)
        return _apply_diagonal_split(product, b.index)
    end
    return product
end
function Base.:*(a::QSym, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.arguments]
    product = QAdd(result, copy(b.indices), copy(b.non_equal))
    if has_index(a.index) && !isempty(b.indices)
        return _apply_diagonal_split(product, a.index)
    end
    return product
end

# QAdd * QMul
function Base.:*(a::QAdd, b::QMul)
    result = QMul[QMul(c, ops) * b for (ops, c) in a.arguments]
    product = QAdd(result, copy(a.indices), copy(a.non_equal))
    if !isempty(a.indices)
        for op in b.args_nc
            if has_index(op.index)
                product = _apply_diagonal_split(product, op.index)
            end
        end
    end
    return product
end
function Base.:*(a::QMul, b::QAdd)
    result = QMul[a * QMul(c, ops) for (ops, c) in b.arguments]
    product = QAdd(result, copy(b.indices), copy(b.non_equal))
    if !isempty(b.indices)
        for op in a.args_nc
            if has_index(op.index)
                product = _apply_diagonal_split(product, op.index)
            end
        end
    end
    return product
end

# QAdd * QAdd
function Base.:*(a::QAdd, b::QAdd)
    # Check for clashing summation indices
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            if idx in b.indices
                throw(ArgumentError(
                    "Summation index $(idx.name) appears in both factors. " *
                    "Use `change_index` to re-index one side, or use different indices."
                ))
            end
        end
    end
    result = QMul[]
    for (ops_a, c_a) in a.arguments, (ops_b, c_b) in b.arguments
        push!(result, QMul(c_a, ops_a) * QMul(c_b, ops_b))
    end
    indices = _merge_unique(a.indices, b.indices)
    non_equal = _merge_unique(a.non_equal, b.non_equal)
    product = QAdd(result, indices, non_equal)
    # Apply diagonal splitting for cross-space index pairs
    if !isempty(a.indices) && !isempty(b.indices)
        for idx_a in a.indices, idx_b in b.indices
            if idx_a.space_index == idx_b.space_index
                product = _apply_diagonal_split(product, idx_a)
                product = _apply_diagonal_split(product, idx_b)
            end
        end
    end
    return product
end

# QAdd / Number
Base.:/(a::QAdd, b::Number) = a * inv(b)


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
    diag_terms, ne_pairs = _diagonal_split([expr], i, ne_pairs)
    off_diag = QAdd(QMul[expr], [i], ne_pairs)
    isempty(diag_terms) && return off_diag
    return off_diag + QAdd(diag_terms)
end
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    if !_any_depends_on_index(expr, i)
        return expr * i.range
    end
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    all_ne = vcat(expr.non_equal, ne_pairs)
    diag_terms, all_ne = _diagonal_split(terms(expr), i, all_ne)
    all_indices = vcat(expr.indices, [i])
    off_diag = QAdd(expr.arguments, all_indices, all_ne)
    isempty(diag_terms) && return off_diag
    return off_diag + QAdd(diag_terms)
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
    expand_sums(expr)

No-op passthrough. Diagonal splitting is now performed eagerly at construction time
by `Σ` and `QAdd` multiplication. Kept for API compatibility.
"""
expand_sums(s::QAdd) = s
expand_sums(m::QMul) = m
expand_sums(op::QSym) = op
