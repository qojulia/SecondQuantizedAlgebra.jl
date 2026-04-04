"""
    simplify(s::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
Remove zero-prefactor terms. Returns `QAdd`.
"""
function simplify(s::QAdd{T}) where {T}
    # Group by args_nc identity
    groups = Dict{UInt,Tuple{Vector{QSym},Any}}()
    order = UInt[]

    for term in s.arguments
        key = hash(term.args_nc)
        if haskey(groups, key)
            existing_ops, existing_c = groups[key]
            if isequal(existing_ops, term.args_nc)
                groups[key] = (existing_ops, existing_c + term.arg_c)
            else
                # Hash collision fallback
                key2 = hash(term.args_nc, key)
                if haskey(groups, key2)
                    _, ec = groups[key2]
                    groups[key2] = (term.args_nc, ec + term.arg_c)
                else
                    groups[key2] = (term.args_nc, term.arg_c)
                    push!(order, key2)
                end
            end
        else
            groups[key] = (term.args_nc, term.arg_c)
            push!(order, key)
        end
    end

    # Collect results, determining the output prefactor type
    terms = Tuple{Vector{QSym},Any}[]
    for key in order
        ops, c = groups[key]
        iszero(c) && continue
        push!(terms, (ops, c))
    end

    if isempty(terms)
        return QAdd(QMul{T}[QMul(zero(T), QSym[])])
    end

    # Determine common type for prefactors
    TT = typeof(terms[1][2])
    for i in 2:length(terms)
        TT = promote_type(TT, typeof(terms[i][2]))
    end
    result = QMul{TT}[QMul(convert(TT, c), ops) for (ops, c) in terms]
    return QAdd(result)
end

simplify(m::QMul) = QAdd([m])
simplify(op::QSym) = QAdd(QMul{Int}[QMul(1, QSym[op])])

# SymbolicUtils.simplify — also simplify each prefactor
function SymbolicUtils.simplify(s::QAdd{T}; kwargs...) where {T}
    simplified = simplify(s)
    result = QMul{T}[
        QMul(convert(T, _simplify_prefactor(t.arg_c; kwargs...)), t.args_nc)
        for t in simplified.arguments
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
function Symbolics.expand(s::QAdd{T}; kwargs...) where {T}
    result = QMul{T}[
        QMul(convert(T, _expand_prefactor(t.arg_c; kwargs...)), t.args_nc)
        for t in s.arguments
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
