"""
    Σ(expr, i::Index, non_equal::Vector{Index}=Index[])

Create a symbolic sum over index `i`. Returns a `QAdd` with summation indices.

    Σ(expr, i::Index, j::Index, ...)

Create a multi-index sum.
"""
function Σ(expr::QMul, i::Index, non_equal::Vector{Index} = Index[])
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
