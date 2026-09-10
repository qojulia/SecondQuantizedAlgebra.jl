"""
    QExpr

Cold-path representation for formal non-polynomial operator expressions.

`QAdd` remains the canonical polynomial representation. A `QExpr` is entered only
when an operation such as [`sin`](@ref), [`cos`](@ref), or [`expim`](@ref) cannot
be represented by `QAdd` without an infinite expansion. Products preserve factor
order and are never distributed over formal sums implicitly.
"""
@enum QExprKind::UInt8 begin
    QEXPR_ADD
    QEXPR_MUL
    QEXPR_SIN
    QEXPR_COS
    QEXPR_EXPIM
end

struct QExpr <: QField
    kind::QExprKind
    coeff::CNum
    args::Vector{Union{QAdd, QExpr}}
    function QExpr(kind::QExprKind, coeff::CNum, args::Vector{Union{QAdd, QExpr}})
        iszero_cnum(coeff) && throw(
            ArgumentError("QExpr outer coefficient must be nonzero; use the canonical zero"),
        )
        if kind == QEXPR_ADD
            isequal(coeff, CNUM_ONE) || throw(
                ArgumentError("QExpr additive nodes require a unit outer coefficient"),
            )
        elseif kind == QEXPR_SIN || kind == QEXPR_COS || kind == QEXPR_EXPIM
            length(args) == 1 || throw(
                ArgumentError("formal function node must have exactly one argument"),
            )
        end
        return new(kind, coeff, args)
    end
end

const QExprArg = Union{QAdd, QExpr}

qexpr_zero() = QExpr(QEXPR_ADD, CNUM_ONE, QExprArg[])
qexpr_one() = QExpr(QEXPR_MUL, CNUM_ONE, QExprArg[])
qexpr_is_scalar(q::QExpr) = q.kind == QEXPR_MUL && isempty(q.args)

Base.zero(::Type{QExpr}) = qexpr_zero()
Base.zero(::QExpr) = qexpr_zero()
Base.one(::Type{QExpr}) = qexpr_one()
Base.one(::QExpr) = qexpr_one()
Base.iszero(q::QExpr) = q.kind == QEXPR_ADD && isempty(q.args)
Base.isone(q::QExpr) = qexpr_is_scalar(q) && isequal(q.coeff, CNUM_ONE)

qexpr_arg(q::QSym) = +q
qexpr_arg(q::QAdd) = q
qexpr_arg(q::QExpr) = q

function qexpr_scalar_coeff(q::QAdd)::Union{Nothing, CNum}
    isempty(q.indices) || return nothing
    length(q.arguments) == 1 || return nothing
    term, coeff = first(q.arguments)
    isempty(term.ops) || return nothing
    isempty(term.ne) || return nothing
    return coeff
end

function qexpr_from_qadd(q::QAdd)
    iszero(q) && return qexpr_zero()
    scalar = qexpr_scalar_coeff(q)
    scalar === nothing || return QExpr(QEXPR_MUL, scalar, QExprArg[])
    return QExpr(QEXPR_MUL, CNUM_ONE, QExprArg[q])
end

qexpr_same_body(a::QExpr, b::QExpr) =
    a.kind == b.kind && isequal(a.args, b.args)

function qexpr_scale(q::QExpr, c::CNum)
    iszero_cnum(c) && return qexpr_zero()
    isequal(c, CNUM_ONE) && return q
    if q.kind == QEXPR_ADD
        args = QExprArg[]
        sizehint!(args, length(q.args))
        for arg in q.args
            push!(args, arg isa QAdd ? arg * c : qexpr_scale(arg, c))
        end
        return qexpr_sum(args)
    end
    new_c = mul_cnum(q.coeff, c)
    iszero_cnum(new_c) && return qexpr_zero()
    return QExpr(q.kind, new_c, copy(q.args))
end

function qexpr_sum(raw::Vector{QExprArg})
    flat = QExprArg[]
    for arg in raw
        if arg isa QAdd
            iszero(arg) || push!(flat, arg)
        elseif iszero(arg)
            continue
        elseif arg.kind == QEXPR_ADD
            append!(flat, arg.args)
        else
            push!(flat, arg)
        end
    end

    polynomial = zero(QAdd)
    have_polynomial = false
    formal = QExpr[]
    for arg in flat
        if arg isa QAdd
            polynomial = have_polynomial ? polynomial + arg : arg
            have_polynomial = true
            continue
        elseif qexpr_is_scalar(arg)
            scalar = single_qadd(arg.coeff, Op[])
            polynomial = have_polynomial ? polynomial + scalar : scalar
            have_polynomial = true
            continue
        end
        found = findfirst(term -> qexpr_same_body(term, arg), formal)
        if found === nothing
            push!(formal, arg)
            continue
        end
        old = formal[found]
        new_c = add_cnum(old.coeff, arg.coeff)
        if iszero_cnum(new_c)
            deleteat!(formal, found)
        else
            formal[found] = QExpr(old.kind, new_c, copy(old.args))
        end
    end

    sort!(formal)
    args = QExprArg[]
    have_polynomial && !iszero(polynomial) && push!(args, polynomial)
    append!(args, formal)

    isempty(args) && return qexpr_zero()
    if length(args) == 1
        arg = only(args)
        return arg isa QExpr ? arg : qexpr_from_qadd(arg)
    end
    return QExpr(QEXPR_ADD, CNUM_ONE, args)
end

function qexpr_push_product!(factors::Vector{QExprArg}, coeff::CNum, arg::QAdd)
    iszero(arg) && return (CNUM_ZERO, true)
    scalar = qexpr_scalar_coeff(arg)
    if scalar !== nothing
        coeff = mul_cnum(coeff, scalar)
        return (coeff, iszero_cnum(coeff))
    end
    if !isempty(factors) && last(factors) isa QAdd
        left = pop!(factors)::QAdd
        merged = left * arg
        iszero(merged) && return (CNUM_ZERO, true)
        merged_scalar = qexpr_scalar_coeff(merged)
        if merged_scalar === nothing
            push!(factors, merged)
        else
            coeff = mul_cnum(coeff, merged_scalar)
        end
    else
        push!(factors, arg)
    end
    return (coeff, iszero_cnum(coeff))
end

function qexpr_push_product!(factors::Vector{QExprArg}, coeff::CNum, arg::QExpr)
    iszero(arg) && return (CNUM_ZERO, true)
    if arg.kind == QEXPR_MUL
        coeff = mul_cnum(coeff, arg.coeff)
        iszero_cnum(coeff) && return (CNUM_ZERO, true)
        for factor in arg.args
            coeff, stopped = qexpr_push_product!(factors, coeff, factor)
            stopped && return (CNUM_ZERO, true)
        end
        return (coeff, false)
    end
    coeff = mul_cnum(coeff, arg.coeff)
    iszero_cnum(coeff) && return (CNUM_ZERO, true)
    bare = isequal(arg.coeff, CNUM_ONE) ? arg : QExpr(arg.kind, CNUM_ONE, copy(arg.args))
    push!(factors, bare)
    return (coeff, false)
end

function qexpr_product(raw::Vector{QExprArg}, coefficient::CNum = CNUM_ONE)
    factors = QExprArg[]
    coeff = coefficient
    for arg in raw
        coeff, stopped = qexpr_push_product!(factors, coeff, arg)
        stopped && return qexpr_zero()
    end
    iszero_cnum(coeff) && return qexpr_zero()

    isempty(factors) && return QExpr(QEXPR_MUL, coeff, QExprArg[])
    if length(factors) == 1
        factor = only(factors)
        return factor isa QExpr ? qexpr_scale(factor, coeff) : qexpr_from_qadd(factor * coeff)
    end
    return QExpr(QEXPR_MUL, coeff, factors)
end

function qexpr_call(kind::QExprKind, arg::QExprArg)
    if iszero(arg)
        kind == QEXPR_SIN && return qexpr_zero()
        (kind == QEXPR_COS || kind == QEXPR_EXPIM) && return qexpr_one()
    end
    return QExpr(kind, CNUM_ONE, QExprArg[arg])
end

Base.sin(q::QField) = qexpr_call(QEXPR_SIN, qexpr_arg(q))
Base.cos(q::QField) = qexpr_call(QEXPR_COS, qexpr_arg(q))

function qexpr_provably_hermitian(q::QSym)
    return iszero(simplify((+q) - (+adjoint(q))))
end
function qexpr_provably_hermitian(q::QAdd)
    return iszero(simplify(q - adjoint(q)))
end
qexpr_provably_hermitian(q::QExpr) = isequal(q, adjoint(q))

"""
    expim(A::QField)

Represent the formal unitary phase `exp(im*A)` for a provably Hermitian operator
expression `A`. Unlike [`UnitaryTransform`](@ref), this is the operator itself,
not its compiled adjoint action.
"""
function expim(q::QField)
    qexpr_provably_hermitian(q) || throw(
        ArgumentError("`expim(::QField)` requires a provably Hermitian argument"),
    )
    return qexpr_call(QEXPR_EXPIM, qexpr_arg(q))
end

Base.:+(a::QExpr) = a
Base.:+(a::QExpr, b::QExpr) = qexpr_sum(QExprArg[a, b])
Base.:+(a::QExpr, b::QAdd) = qexpr_sum(QExprArg[a, b])
Base.:+(a::QAdd, b::QExpr) = qexpr_sum(QExprArg[a, b])
Base.:+(a::QExpr, b::QSym) = qexpr_sum(QExprArg[a, +b])
Base.:+(a::QSym, b::QExpr) = qexpr_sum(QExprArg[+a, b])
Base.:+(a::QExpr, b::Coefficient) = qexpr_sum(QExprArg[a, single_qadd(to_cnum(b), Op[])])
Base.:+(a::Coefficient, b::QExpr) = b + a

Base.:-(a::QExpr) = qexpr_scale(a, CNUM_NEG1)

Base.:*(a::QExpr, b::QExpr) = qexpr_product(QExprArg[a, b])
Base.:*(a::QExpr, b::QAdd) = qexpr_product(QExprArg[a, b])
Base.:*(a::QAdd, b::QExpr) = qexpr_product(QExprArg[a, b])
Base.:*(a::QExpr, b::QSym) = qexpr_product(QExprArg[a, +b])
Base.:*(a::QSym, b::QExpr) = qexpr_product(QExprArg[+a, b])
Base.:*(a::QExpr, b::Coefficient) = qexpr_scale(a, to_cnum(b))
Base.:*(a::Coefficient, b::QExpr) = qexpr_scale(b, to_cnum(a))

Base.:/(a::QExpr, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? a * (1 // b) : a * inv(b)
Base.:/(a::QExpr, b::Coefficient) = qexpr_scale(a, inv(to_cnum(b)))
Base.://(a::QExpr, b::Integer) = a * (1 // b)
Base.://(a::QExpr, b::Coefficient) = a / b

function Base.:^(a::QExpr, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return qexpr_one()
    n == 1 && return a
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

function Base.isequal(a::QExpr, b::QExpr)
    return a.kind == b.kind && isequal(a.coeff, b.coeff) && isequal(a.args, b.args)
end
Base.:(==)(a::QExpr, b::QExpr) = isequal(a, b)
Base.hash(q::QExpr, h::UInt) = hash(q.args, hash(q.coeff, hash(q.kind, hash(:QExpr, h))))

function qexpr_arg_less(a::QExprArg, b::QExprArg)
    isequal(a, b) && return false
    a isa QAdd && b isa QExpr && return true
    a isa QExpr && b isa QAdd && return false
    if a isa QAdd
        return isless(a, b::QAdd)
    end
    return isless(a::QExpr, b::QExpr)
end

function qexpr_args_less(a::Vector{QExprArg}, b::Vector{QExprArg})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        isequal(a[i], b[i]) && continue
        return qexpr_arg_less(a[i], b[i])
    end
    return length(a) < length(b)
end

function Base.isless(a::QExpr, b::QExpr)
    a.kind != b.kind && return Int(a.kind) < Int(b.kind)
    isequal(a.args, b.args) || return qexpr_args_less(a.args, b.args)
    return isless(coeff_key(a.coeff), coeff_key(b.coeff))
end

qexpr_adjoint_arg(arg::QAdd) = adjoint(arg)
qexpr_adjoint_arg(arg::QExpr) = adjoint(arg)

function Base.adjoint(q::QExpr)
    c = conj_cnum(q.coeff)
    if q.kind == QEXPR_ADD
        args = QExprArg[qexpr_adjoint_arg(arg) for arg in q.args]
        return qexpr_scale(qexpr_sum(args), c)
    elseif q.kind == QEXPR_MUL
        args = QExprArg[]
        sizehint!(args, length(q.args))
        for arg in Iterators.reverse(q.args)
            push!(args, qexpr_adjoint_arg(arg))
        end
        return qexpr_product(args, c)
    elseif q.kind == QEXPR_SIN || q.kind == QEXPR_COS
        arg = qexpr_adjoint_arg(only(q.args))
        return qexpr_scale(qexpr_call(q.kind, arg), c)
    elseif q.kind == QEXPR_EXPIM
        arg = qexpr_adjoint_arg(only(q.args))
        return qexpr_scale(qexpr_call(QEXPR_EXPIM, -arg), c)
    end
    error("unknown QExpr kind $(q.kind)")
end
