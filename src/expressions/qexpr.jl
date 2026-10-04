@enum QExprKind::UInt8 begin
    QEXPR_ADD
    QEXPR_MUL
    QEXPR_SIN
    QEXPR_COS
    QEXPR_EXPIM
end

"""
    QExpr <: QField

Cold-path representation for formal non-polynomial operator expressions.

`QAdd` remains the canonical polynomial representation. A `QExpr` is entered only
when an operation such as `sin`, `cos`, or [`expim`](@ref) cannot
be represented by `QAdd` without an infinite expansion. Products keep the order of
factors that share a Hilbert space, while factors on disjoint spaces commute into a
canonical order. Products are never distributed over formal sums implicitly.

Unary formal nodes store their single `QAdd` or `QExpr` child directly; sums and
products retain dynamic vector storage. This keeps `QExpr` immutable and
non-parametric while avoiding a one-element heap container for unary nodes.
"""
struct QExpr <: QField
    kind::QExprKind
    coeff::CNum
    storage::Union{QAdd, QExpr, Vector{Union{QAdd, QExpr}}}
    function QExpr(
            kind::QExprKind,
            coeff::CNum,
            storage::Union{QAdd, QExpr, Vector{Union{QAdd, QExpr}}},
        )
        iszero_cnum(coeff) && throw(
            ArgumentError("QExpr outer coefficient must be nonzero; use the canonical zero"),
        )
        if kind == QEXPR_ADD
            storage isa Vector || throw(
                ArgumentError("QExpr additive nodes require dynamic argument storage"),
            )
            isequal(coeff, CNUM_ONE) || throw(
                ArgumentError("QExpr additive nodes require a unit outer coefficient"),
            )
        elseif kind == QEXPR_MUL
            storage isa Vector || throw(
                ArgumentError("QExpr product nodes require dynamic argument storage"),
            )
        else
            storage isa Vector && throw(
                ArgumentError("formal function nodes require direct unary argument storage"),
            )
        end
        return new(kind, coeff, storage)
    end
end

const QExprArg = Union{QAdd, QExpr}
const QExprArgs = Union{Tuple{QExprArg}, Vector{QExprArg}}

QExpr(kind::QExprKind, coeff::CNum, args::Tuple{QExprArg}) =
    QExpr(kind, coeff, only(args))

@inline function Base.getproperty(q::QExpr, name::Symbol)
    name === :args || return getfield(q, name)
    storage = getfield(q, :storage)
    return storage isa Vector ? storage : (storage,)
end

Base.propertynames(::QExpr, private::Bool = false) =
    private ? (:kind, :coeff, :args, :storage) : (:kind, :coeff, :args)

@inline function qexpr_copy_storage(q::QExpr)
    storage = getfield(q, :storage)
    return storage isa Vector ? copy(storage) : storage
end

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

# A product node holding a single polynomial factor is a polynomial in formal clothing. Its
# canonical form keeps the scale inside the `QAdd`, and sums fold it into the polynomial part.
function qexpr_polynomial_body(q::QExpr)::Union{Nothing, QAdd}
    q.kind == QEXPR_MUL || return nothing
    storage = getfield(q, :storage)::Vector{QExprArg}
    length(storage) == 1 || return nothing
    body = only(storage)
    return body isa QAdd ? body : nothing
end

qexpr_same_body(a::QExpr, b::QExpr) =
    a.kind == b.kind && isequal(getfield(a, :storage), getfield(b, :storage))

function qexpr_scale(q::QExpr, c::CNum)
    iszero_cnum(c) && return qexpr_zero()
    isequal(c, CNUM_ONE) && return q
    body = qexpr_polynomial_body(q)
    body === nothing || return qexpr_from_qadd(body * mul_cnum(q.coeff, c))
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
    return QExpr(q.kind, new_c, qexpr_copy_storage(q))
end

function qexpr_collect_sum!(
        formal::Vector{QExpr}, polynomial::Union{Nothing, QAddBuilder}, arg::QAdd,
    )
    iszero(arg) && return polynomial
    polynomial === nothing && (polynomial = QAddBuilder())
    accumulate!(polynomial, arg)
    return polynomial
end

function qexpr_collect_sum!(
        formal::Vector{QExpr}, polynomial::Union{Nothing, QAddBuilder}, arg::QExpr,
    )
    iszero(arg) && return polynomial
    if arg.kind == QEXPR_ADD
        for child in arg.args
            polynomial = qexpr_collect_sum!(formal, polynomial, child)
        end
    elseif qexpr_is_scalar(arg)
        polynomial === nothing && (polynomial = QAddBuilder())
        accumulate!(polynomial, single_qadd(arg.coeff, EMPTY_OPS))
    elseif (body = qexpr_polynomial_body(arg)) !== nothing
        polynomial === nothing && (polynomial = QAddBuilder())
        accumulate!(polynomial, body, arg.coeff)
    else
        push!(formal, arg)
    end
    return polynomial
end

function qexpr_coalesce_formal(formal::Vector{QExpr})
    length(formal) <= 1 && return formal
    if length(formal) == 2
        a, b = formal
        if qexpr_same_body(a, b)
            new_c = add_cnum(a.coeff, b.coeff)
            return iszero_cnum(new_c) ? QExpr[] :
                QExpr[QExpr(a.kind, new_c, qexpr_copy_storage(a))]
        end
        isless(b, a) && ((formal[1], formal[2]) = (b, a))
        return formal
    end

    sort!(formal)
    result = QExpr[]
    sizehint!(result, length(formal))
    for arg in formal
        if isempty(result) || !qexpr_same_body(last(result), arg)
            push!(result, arg)
            continue
        end
        old = pop!(result)
        new_c = add_cnum(old.coeff, arg.coeff)
        iszero_cnum(new_c) ||
            push!(result, QExpr(old.kind, new_c, qexpr_copy_storage(old)))
    end
    return result
end

function qexpr_sum(raw::Vector{QExprArg})
    polynomial = nothing
    formal = QExpr[]
    sizehint!(formal, length(raw))
    for arg in raw
        polynomial = qexpr_collect_sum!(formal, polynomial, arg)
    end

    polynomial_value = polynomial === nothing ? zero_qadd() : build(polynomial)
    formal = qexpr_coalesce_formal(formal)
    args = QExprArg[]
    iszero(polynomial_value) || push!(args, polynomial_value)
    append!(args, formal)

    isempty(args) && return qexpr_zero()
    if length(args) == 1
        arg = only(args)
        return arg isa QExpr ? arg : qexpr_from_qadd(arg)
    end
    return QExpr(QEXPR_ADD, CNUM_ONE, args)
end

# A scalar can sit in the outer coefficient or inside a polynomial factor, so `(θa) * f` and
# `θ * (a * f)` would differ. Polynomial factors therefore hand their content to the outer
# coefficient: a coefficient shared by every term, else a native leading coefficient (exact to
# divide by). A factor with non-uniform symbolic coefficients and a symbolic leading term is
# left as is, since symbolic division is not exact.
function qexpr_split_content(q::QAdd)::Tuple{CNum, QAdd}
    content = CNUM_ONE
    uniform = true
    for (i, c) in enumerate(values(q.arguments))
        if i == 1
            content = c
        elseif !isequal(c, content)
            uniform = false
            break
        end
    end
    if !uniform
        lead = argmin(term_order_key, keys(q.arguments))
        content = q.arguments[lead]
        is_native(content) || return (CNUM_ONE, q)
    end
    isequal(content, CNUM_ONE) && return (CNUM_ONE, q)
    scale = uniform ? CNUM_ONE : inv(content)
    d = QTermDict()
    sizehint!(d, length(q.arguments))
    for (term, c) in q.arguments
        d[copy_key(term)] = uniform ? CNUM_ONE : mul_cnum(c, scale)
    end
    return (content, QAdd(d, q.indices))
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
            content, merged = qexpr_split_content(merged)
            coeff = mul_cnum(coeff, content)
            push!(factors, merged)
        else
            coeff = mul_cnum(coeff, merged_scalar)
        end
    else
        content, arg = qexpr_split_content(arg)
        coeff = mul_cnum(coeff, content)
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
    bare = isequal(arg.coeff, CNUM_ONE) ? arg :
        QExpr(arg.kind, CNUM_ONE, qexpr_copy_storage(arg))
    push!(factors, bare)
    return (coeff, false)
end

# Factors acting on disjoint spaces commute, since no operator kind is fermionic. Choosing the
# lexicographically least ordering among the commutation-equivalent ones (the trace-monoid
# normal form) makes `b * f(a)` and `f(a) * b` the same term. Factors sharing a space keep
# their relative order.
function qexpr_lex_order(factors::Vector{QExprArg})::Vector{Int}
    spaces = [acts_on(f) for f in factors]
    remaining = collect(eachindex(factors))
    order = Int[]
    sizehint!(order, length(factors))
    while !isempty(remaining)
        best = 0
        for (pos, j) in enumerate(remaining)
            movable = true
            for k in view(remaining, 1:(pos - 1))
                isdisjoint(spaces[k], spaces[j]) || (movable = false; break)
            end
            movable || continue
            if best == 0 || qexpr_arg_less(factors[j], factors[remaining[best]])
                best = pos
            end
        end
        push!(order, remaining[best])
        deleteat!(remaining, best)
    end
    return order
end

# Reordering can make two polynomial factors adjacent; they merge through `QAdd`
# multiplication, which may enable a further reordering, so iterate to a fixed point.
function qexpr_commuting_order(factors::Vector{QExprArg}, coeff::CNum)
    while length(factors) > 1
        order = qexpr_lex_order(factors)
        issorted(order) && break
        reordered = QExprArg[]
        sizehint!(reordered, length(factors))
        for i in order
            coeff, stopped = qexpr_push_product!(reordered, coeff, factors[i])
            stopped && return (factors, CNUM_ZERO, true)
        end
        factors = reordered
    end
    return (factors, coeff, false)
end

function qexpr_product(raw::Vector{QExprArg}, coefficient::CNum = CNUM_ONE)
    factors = QExprArg[]
    sizehint!(factors, length(raw))
    coeff = coefficient
    for arg in raw
        coeff, stopped = qexpr_push_product!(factors, coeff, arg)
        stopped && return qexpr_zero()
    end
    iszero_cnum(coeff) && return qexpr_zero()
    factors, coeff, stopped = qexpr_commuting_order(factors, coeff)
    stopped && return qexpr_zero()

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
    return QExpr(kind, CNUM_ONE, arg)
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
Base.inv(::QExpr) = throw(ArgumentError("Negative powers not supported"))

function Base.:^(a::QExpr, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return qexpr_one()
    n == 1 && return a
    n <= typemax(Int) || throw(ArgumentError("Power is too large for this platform"))
    factors = Vector{QExprArg}(undef, Int(n))
    fill!(factors, a)
    return qexpr_product(factors)
end

function Base.isequal(a::QExpr, b::QExpr)
    return a.kind == b.kind && isequal(a.coeff, b.coeff) &&
        isequal(getfield(a, :storage), getfield(b, :storage))
end
Base.:(==)(a::QExpr, b::QExpr) = isequal(a, b)
Base.hash(q::QExpr, h::UInt) =
    hash(getfield(q, :storage), hash(q.coeff, hash(q.kind, hash(:QExpr, h))))

function qexpr_arg_less(a::QExprArg, b::QExprArg)
    isequal(a, b) && return false
    a isa QAdd && b isa QExpr && return true
    a isa QExpr && b isa QAdd && return false
    if a isa QAdd
        return isless(a, b::QAdd)
    end
    return isless(a::QExpr, b::QExpr)
end

function qexpr_args_less(a::QExprArgs, b::QExprArgs)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        isequal(a[i], b[i]) && continue
        return qexpr_arg_less(a[i], b[i])
    end
    return length(a) < length(b)
end

function qexpr_storage_less(a::QExpr, b::QExpr)
    a_storage = getfield(a, :storage)
    b_storage = getfield(b, :storage)
    if a_storage isa Vector
        return qexpr_args_less(a_storage, b_storage::Vector{QExprArg})
    end
    return qexpr_arg_less(a_storage::QExprArg, b_storage::QExprArg)
end

function Base.isless(a::QExpr, b::QExpr)
    a.kind != b.kind && return Int(a.kind) < Int(b.kind)
    a_storage = getfield(a, :storage)
    b_storage = getfield(b, :storage)
    isequal(a_storage, b_storage) || return qexpr_storage_less(a, b)
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
        arg = qexpr_adjoint_arg(getfield(q, :storage)::QExprArg)
        return qexpr_scale(qexpr_call(q.kind, arg), c)
    end
    arg = qexpr_adjoint_arg(getfield(q, :storage)::QExprArg)
    return qexpr_scale(qexpr_call(QEXPR_EXPIM, -arg), c)
end
