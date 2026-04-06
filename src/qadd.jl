"""
    QAdd <: QTerm

Sum of [`QMul`](@ref) terms, optionally with summation indices for symbolic sums.

Fields:
- `arguments::Vector{QMul}` — terms of the sum
- `indices::Vector{Index}` — summation indices (empty = regular sum)
- `non_equal::Vector{Tuple{Index,Index}}` — pairwise inequality constraints
"""
struct QAdd <: QTerm
    arguments::Vector{QMul}
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
    function QAdd(
            args::Vector{QMul}, indices::Vector{Index},
            non_equal::Vector{Tuple{Index, Index}}
        )
        return new(args, indices, non_equal)
    end
end
QAdd(args::Vector{QMul}) = QAdd(args, Index[], Tuple{Index, Index}[])

Base.length(a::QAdd) = length(a.arguments)

# Equality and hashing
function Base.isequal(a::QAdd, b::QAdd)
    length(a.arguments) == length(b.arguments) || return false
    a.indices == b.indices || return false
    a.non_equal == b.non_equal || return false
    for (x, y) in zip(a.arguments, b.arguments)
        isequal(x, y) || return false
    end
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
Base.hash(q::QAdd, h::UInt) = hash(:QAdd, hash(q.arguments, hash(q.indices, hash(q.non_equal, h))))

# Adjoint
function Base.adjoint(q::QAdd)
    return QAdd(QMul[adjoint(t) for t in q.arguments], q.indices, q.non_equal)
end

# Helpers: wrap QSym/scalar as QMul
_to_qmul(a::QSym) = QMul(_to_cnum(1), QSym[a])
_to_qmul(a::QMul) = a
_scalar_qmul(x::Number) = QMul(_to_cnum(x), QSym[])

## Addition — always returns QAdd

# QSym + QSym
function Base.:+(a::QSym, b::QSym)
    return QAdd(QMul[_to_qmul(a), _to_qmul(b)])
end

# QMul + QMul
function Base.:+(a::QMul, b::QMul)
    return QAdd(QMul[a, b])
end

# QMul + QSym
function Base.:+(a::QMul, b::QSym)
    return QAdd(QMul[a, _to_qmul(b)])
end
Base.:+(a::QSym, b::QMul) = b + a

# QAdd + QMul
function Base.:+(a::QAdd, b::QMul)
    args = QMul[x for x in a.arguments]
    push!(args, b)
    return QAdd(args, a.indices, a.non_equal)
end
Base.:+(a::QMul, b::QAdd) = b + a

# QAdd + QSym
function Base.:+(a::QAdd, b::QSym)
    return a + _to_qmul(b)
end
Base.:+(a::QSym, b::QAdd) = b + a

# QAdd + QAdd
function Base.:+(a::QAdd, b::QAdd)
    args = QMul[x for x in a.arguments]
    for x in b.arguments
        push!(args, x)
    end
    # Merge indices if both have the same (or combine)
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(args, indices, non_equal)
end

# QField + Number
function Base.:+(a::QSym, b::Number)
    return QAdd(QMul[_to_qmul(a), _scalar_qmul(b)])
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QMul, b::Number)
    return QAdd(QMul[a, _scalar_qmul(b)])
end
Base.:+(a::Number, b::QMul) = b + a

function Base.:+(a::QAdd, b::Number)
    args = QMul[x for x in a.arguments]
    push!(args, _scalar_qmul(b))
    return QAdd(args, a.indices, a.non_equal)
end
Base.:+(a::Number, b::QAdd) = b + a

# Subtraction
Base.:-(a::QAdd) = QAdd(QMul[-t for t in a.arguments], a.indices, a.non_equal)
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

## QAdd * ... (distributive)

# QAdd * Number
function Base.:*(a::QAdd, b::Number)
    return QAdd(QMul[t * b for t in a.arguments], a.indices, a.non_equal)
end
Base.:*(a::Number, b::QAdd) = b * a

# QAdd * QSym
function Base.:*(a::QAdd, b::QSym)
    return QAdd(QMul[t * b for t in a.arguments], a.indices, a.non_equal)
end
function Base.:*(a::QSym, b::QAdd)
    return QAdd(QMul[a * t for t in b.arguments], b.indices, b.non_equal)
end

# QAdd * QMul
function Base.:*(a::QAdd, b::QMul)
    return QAdd(QMul[t * b for t in a.arguments], a.indices, a.non_equal)
end
function Base.:*(a::QMul, b::QAdd)
    return QAdd(QMul[a * t for t in b.arguments], b.indices, b.non_equal)
end

# QAdd * QAdd
function Base.:*(a::QAdd, b::QAdd)
    args = QMul[]
    for ai in a.arguments, bi in b.arguments
        push!(args, ai * bi)
    end
    indices = vcat(a.indices, b.indices) |> unique
    non_equal = vcat(a.non_equal, b.non_equal) |> unique
    return QAdd(args, indices, non_equal)
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
