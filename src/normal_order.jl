"""
    normal_order(expr)

Apply bosonic commutation relation [a, a†] = 1 to rewrite the expression
in normal order (all creation operators to the left of annihilation operators).
Always returns `QAdd`.

Does NOT collect like terms — use `simplify()` for that.
"""
function normal_order(op::QSym)
    return QAdd(QMul{Int}[QMul(1, QSym[op])])
end

function normal_order(m::QMul{T}) where {T}
    return _normal_order_qmul(m.arg_c, m.args_nc)
end

function normal_order(s::QAdd{T}) where {T}
    result = QMul{T}[]
    for term in s.arguments
        ordered = normal_order(term)
        for t in ordered.arguments
            push!(result, convert(QMul{T}, t))
        end
    end
    return QAdd(result)
end

"""
    _normal_order_qmul(prefactor, ops)

Recursively normal-order a product of operators. Scans for adjacent
Destroy-Create pairs on the same space, swaps and adds identity term.
"""
function _normal_order_qmul(arg_c::T, ops::Vector{QSym}) where {T}
    isempty(ops) && return QAdd(QMul{T}[QMul(arg_c, QSym[])])
    length(ops) == 1 && return QAdd(QMul{T}[QMul(arg_c, copy(ops))])

    # Find first out-of-normal-order pair on the same space
    for i in 1:(length(ops) - 1)
        a, b = ops[i], ops[i + 1]
        if a isa Destroy && b isa Create && a.space_index == b.space_index && a.name == b.name
            # Apply [a, a†] = 1: swap and add identity

            # Term 1: swapped pair
            swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
            sort!(swapped; lt = canonical_lt)
            term1 = _normal_order_qmul(arg_c, swapped)

            # Term 2: contracted (remove both operators)
            contracted = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
            sort!(contracted; lt = canonical_lt)
            term2 = _normal_order_qmul(arg_c, contracted)

            return QAdd(QMul{T}[term1.arguments..., term2.arguments...])
        end
    end

    # Already in normal order
    return QAdd(QMul{T}[QMul(arg_c, ops)])
end
