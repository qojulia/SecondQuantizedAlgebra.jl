"""
    simplify(s::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
Remove zero-prefactor terms. Returns `QAdd`.
"""
function simplify(s::QAdd{T}) where {T}
    # Collect unique operator patterns and their accumulated prefactors
    unique_ops = Vector{QSym}[]
    prefactors = T[]

    for term in s.arguments
        idx = _find_matching_ops(unique_ops, term.args_nc)
        if idx === nothing
            push!(unique_ops, term.args_nc)
            push!(prefactors, term.arg_c)
        else
            prefactors[idx] = prefactors[idx] + term.arg_c
        end
    end

    # Build result, skipping zeros
    # Determine output type (prefactor addition may widen type, e.g. Int -> Num)
    result_ops = Vector{QSym}[]
    result_cs = []
    for (ops, c) in zip(unique_ops, prefactors)
        iszero(c) && continue
        push!(result_ops, ops)
        push!(result_cs, c)
    end

    if isempty(result_cs)
        return QAdd(QMul{T}[QMul(zero(T), QSym[])])
    end

    TT = promote_type((typeof(c) for c in result_cs)...)
    result = QMul{TT}[QMul(convert(TT, c), ops) for (ops, c) in zip(result_ops, result_cs)]
    return QAdd(result)
end

function _find_matching_ops(unique_ops::Vector{Vector{QSym}}, target::Vector{QSym})
    for (i, ops) in enumerate(unique_ops)
        isequal(ops, target) && return i
    end
    return nothing
end

simplify(m::QMul) = QAdd([m])
simplify(op::QSym) = QAdd(QMul{Int}[QMul(1, QSym[op])])

# SymbolicUtils.simplify — also simplify each prefactor
function SymbolicUtils.simplify(s::QAdd; kwargs...)
    simplified = simplify(s)
    args = simplified.arguments
    TT = promote_type((typeof(_simplify_prefactor(t.arg_c; kwargs...)) for t in args)...)
    result = QMul{TT}[
        QMul(convert(TT, _simplify_prefactor(t.arg_c; kwargs...)), t.args_nc)
        for t in args
    ]
    return QAdd(result)
end
function SymbolicUtils.simplify(m::QMul; kwargs...)
    return SymbolicUtils.simplify(QAdd([m]); kwargs...)
end
function SymbolicUtils.simplify(op::QSym; kwargs...)
    return SymbolicUtils.simplify(simplify(op); kwargs...)
end

_simplify_prefactor(x::Number; kwargs...) = x
_simplify_prefactor(x; kwargs...) = SymbolicUtils.simplify(x; kwargs...)

# Symbolics.expand — distribute products, then expand prefactors
function Symbolics.expand(s::QAdd; kwargs...)
    args = s.arguments
    TT = promote_type((typeof(_expand_prefactor(t.arg_c; kwargs...)) for t in args)...)
    result = QMul{TT}[
        QMul(convert(TT, _expand_prefactor(t.arg_c; kwargs...)), t.args_nc)
        for t in args
    ]
    return QAdd(result)
end
function Symbolics.expand(m::QMul; kwargs...)
    return Symbolics.expand(QAdd([m]); kwargs...)
end
function Symbolics.expand(op::QSym; kwargs...)
    return QAdd(QMul{Int}[QMul(1, QSym[op])])
end

_expand_prefactor(x::Number; kwargs...) = x
_expand_prefactor(x; kwargs...) = Symbolics.expand(x; kwargs...)
