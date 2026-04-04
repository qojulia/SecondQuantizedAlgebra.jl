"""
    simplify(expr, [ordering])

Simplify a quantum expression by:
1. Applying ordering-independent algebraic reductions (Transition composition/orthogonality, Pauli product rule)
2. Applying ordering-dependent commutation swaps (default: `NormalOrder()`)
3. Collecting like terms and summing prefactors

Returns `QAdd`. The `ordering` argument controls which commutation rules are applied.
Default is `NormalOrder()`.
"""
simplify(expr::QField) = simplify(expr, NormalOrder())
simplify(expr::QField, h::HilbertSpace) = simplify(expr, NormalOrder(), h)

function simplify(op::QSym, ::OrderingConvention)
    return QAdd(QMul{Int}[QMul(1, QSym[op])])
end

function simplify(m::QMul, ord::OrderingConvention)
    expanded = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return _collect_like_terms(expanded)
end

function simplify(s::QAdd, ord::OrderingConvention)
    all_terms = QMul[]
    for term in s.arguments
        expanded = _simplify_qmul(term.arg_c, term.args_nc, ord)
        append!(all_terms, expanded.arguments)
    end
    return _collect_like_terms(_collect_qadd(all_terms))
end

# With Hilbert space (ground state rewriting)
function simplify(op::QSym, ord::OrderingConvention, h::HilbertSpace)
    return _collect_like_terms(_apply_ground_state(simplify(op, ord), h))
end
function simplify(m::QMul, ord::OrderingConvention, h::HilbertSpace)
    expanded = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return _collect_like_terms(_apply_ground_state(_collect_qadd(expanded.arguments), h))
end
function simplify(s::QAdd, ord::OrderingConvention, h::HilbertSpace)
    all_terms = QMul[]
    for term in s.arguments
        expanded = _simplify_qmul(term.arg_c, term.args_nc, ord)
        append!(all_terms, expanded.arguments)
    end
    return _collect_like_terms(_apply_ground_state(_collect_qadd(all_terms), h))
end

# Helper: collect QMul[] into a properly typed QAdd
function _collect_qadd(terms::Vector{<:QMul})
    isempty(terms) && return QAdd(QMul{Int}[QMul(0, QSym[])])
    TT = promote_type((typeof(t.arg_c) for t in terms)...)
    return QAdd(QMul{TT}[convert(QMul{TT}, t) for t in terms])
end

"""
    _simplify_qmul(prefactor, ops, ordering)

Recursively simplify a product of operators. Applies:
1. Ordering-independent reductions (always)
2. Ordering-dependent commutation swaps (depends on `ordering`)
"""
# Note: sort!(... lt=canonical_lt) is a stable sort (Julia guarantee).
# Since canonical_lt only compares space_index, operators on the same space
# preserve their relative order after sorting — this is essential for correctness.
function _simplify_qmul(arg_c, ops::Vector{QSym}, ord::OrderingConvention)
    isempty(ops) && return QAdd([QMul(arg_c, QSym[])])
    length(ops) == 1 && return QAdd([QMul(arg_c, copy(ops))])

    for i in 1:(length(ops) - 1)
        a, b = ops[i], ops[i + 1]
        same_space = a.space_index == b.space_index && a.name == b.name

        ## Ordering-independent reductions (always applied)

        # Transition: |i⟩⟨j| · |k⟩⟨l|
        if a isa Transition && b isa Transition && same_space
            if a.j == b.i
                composed = Transition(a.name, a.i, b.j, a.space_index)
                new_ops = QSym[ops[1:(i - 1)]..., composed, ops[(i + 2):end]...]
                return _simplify_qmul(arg_c, new_ops, ord)
            else
                return QAdd([QMul(zero(arg_c), QSym[])])
            end
        end

        # Pauli: σⱼ·σₖ = δⱼₖ + iϵⱼₖₗσₗ
        if a isa Pauli && b isa Pauli && same_space
            if a.axis == b.axis
                new_ops = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
                return _simplify_qmul(arg_c, new_ops, ord)
            else
                l = 6 - a.axis - b.axis
                eps = levicivita([a.axis, b.axis, l])
                new_op = Pauli(a.name, l, a.space_index)
                new_ops = QSym[ops[1:(i - 1)]..., new_op, ops[(i + 2):end]...]
                return _simplify_qmul(arg_c * im * eps, new_ops, ord)
            end
        end

        ## Ordering-dependent swaps (dispatch on ordering convention)
        result = _apply_ordering_rule(a, b, same_space, i, arg_c, ops, ord)
        result !== nothing && return result
    end

    # Already fully simplified
    return QAdd([QMul(arg_c, ops)])
end

"""
    _apply_ordering_rule(a, b, same_space, i, arg_c, ops, ordering)

Apply ordering-specific commutation rules. Returns `nothing` if no rule applies,
or a `QAdd` with the expanded terms.
"""
function _apply_ordering_rule(a, b, same_space, i, arg_c, ops, ::NormalOrder)
    # Fock: [a, a†] = 1
    if a isa Destroy && b isa Create && same_space
        swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
        sort!(swapped; lt=canonical_lt)
        term1 = _simplify_qmul(arg_c, swapped, NormalOrder())
        contracted = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
        sort!(contracted; lt=canonical_lt)
        term2 = _simplify_qmul(arg_c, contracted, NormalOrder())
        all_terms = QMul[term1.arguments..., term2.arguments...]
        return _collect_qadd(all_terms)
    end

    # Spin: [Sⱼ, Sₖ] = iϵⱼₖₗSₗ (swap out-of-order axes)
    if a isa Spin && b isa Spin && same_space && a.axis > b.axis
        swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
        sort!(swapped; lt=canonical_lt)
        term1 = _simplify_qmul(arg_c, swapped, NormalOrder())
        l = 6 - a.axis - b.axis
        eps = levicivita([a.axis, b.axis, l])
        comm_op = Spin(a.name, l, a.space_index)
        comm_ops = QSym[ops[1:(i - 1)]..., comm_op, ops[(i + 2):end]...]
        sort!(comm_ops; lt=canonical_lt)
        term2 = _simplify_qmul(arg_c * im * eps, comm_ops, NormalOrder())
        all_terms = QMul[term1.arguments..., term2.arguments...]
        return _collect_qadd(all_terms)
    end

    # PhaseSpace: [X, P] = i (swap P·X → X·P - i)
    # Note: Position and Momentum have different names, so we check space_index only
    if a isa Momentum && b isa Position && a.space_index == b.space_index
        swapped = QSym[ops[1:(i - 1)]..., b, a, ops[(i + 2):end]...]
        sort!(swapped; lt=canonical_lt)
        term1 = _simplify_qmul(arg_c, swapped, NormalOrder())
        contracted = QSym[ops[1:(i - 1)]..., ops[(i + 2):end]...]
        sort!(contracted; lt=canonical_lt)
        term2 = _simplify_qmul(-im * arg_c, contracted, NormalOrder())
        all_terms = QMul[term1.arguments..., term2.arguments...]
        return _collect_qadd(all_terms)
    end

    return nothing
end

# TODO: SymmetricOrder for TWA (Truncated Wigner Approximation)
# Symmetric ordering symmetrizes operator products: (ab + ba)/2.
# For Fock: instead of swapping a·a† → a†·a + 1, symmetric ordering
# would use (a·a† + a†·a)/2 = a†·a + 1/2.
# For Spin: symmetric ordering of Sⱼ·Sₖ = (SⱼSₖ + SₖSⱼ)/2.
# Implement by adding:
#   function _apply_ordering_rule(a, b, same_space, i, arg_c, ops, ::SymmetricOrder)
#       ...
#   end

"""
    _collect_like_terms(expr::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
Remove zero-prefactor terms.
"""
function _collect_like_terms(s::QAdd{T}) where {T}
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
