@enum QExprKind::UInt8 begin
    QEXPR_ADD
    QEXPR_MUL
    QEXPR_SIN
    QEXPR_COS
    QEXPR_EXPIM
end

"""
    QExpr <: QField

Type of operator expressions that contain formal functions such as `sin`, `cos` or `expim`.

Polynomials in the operators are kept in canonical form, with products reordered and like
terms collected. An expression containing a formal function cannot be brought into that
form, so it is held as written: products keep the order of factors on a shared subsystem,
factors on different subsystems are sorted, and products are not multiplied out over sums.
Exact operations such as `adjoint`, [`substitute`](@ref) or [`simplify`](@ref) act on the
arguments of the functions. [`taylor`](@ref) turns a `QExpr` back into an ordinary
polynomial. The type is public but not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: QExpr

julia> a = Destroy(FockSpace(:cavity), :a);

julia> x = (a + a') * cos(a' * a);

julia> x isa QExpr, taylor(x, 0:0) isa QExpr
(true, false)

julia> x'
cos(a' * a) * (a + a')
```

See also [`sin`](@ref sin(::QField)), [`cos`](@ref cos(::QField)), [`expim`](@ref expim(::QField)), [`taylor`](@ref).
"""
struct QExpr <: QField
    kind::QExprKind
    coeff::CNum
    args::Vector{Union{QAdd, QExpr}}
    function QExpr(kind::QExprKind, coeff::CNum, args::Vector{Union{QAdd, QExpr}})
        iszero(coeff) && throw(
            ArgumentError("QExpr outer coefficient must be nonzero; use the canonical zero"),
        )
        if kind == QEXPR_ADD
            isequal(coeff, CNUM_ONE) || throw(
                ArgumentError("QExpr additive nodes require a unit outer coefficient"),
            )
        elseif kind == QEXPR_MUL
            length(args) == 1 && only(args) isa QAdd && !isequal(coeff, CNUM_ONE) &&
                throw(ArgumentError("a single-polynomial product requires a unit outer coefficient"))
        else
            length(args) == 1 || throw(
                ArgumentError("formal function nodes take exactly one argument"),
            )
        end
        return new(kind, coeff, args)
    end
end

const QExprArg = Union{QAdd, QExpr}

QExpr(kind::QExprKind, coeff::CNum, arg::QExprArg) = QExpr(kind, coeff, QExprArg[arg])

qexpr_zero() = QExpr(QEXPR_ADD, CNUM_ONE, QExprArg[])
qexpr_one() = QExpr(QEXPR_MUL, CNUM_ONE, QExprArg[])
qexpr_is_scalar(q::QExpr) = q.kind == QEXPR_MUL && isempty(q.args)

Base.zero(::Union{Type{QExpr}, QExpr}) = qexpr_zero()
Base.one(::Union{Type{QExpr}, QExpr}) = qexpr_one()
Base.iszero(q::QExpr) = q.kind == QEXPR_ADD && isempty(q.args)
Base.isone(q::QExpr) = qexpr_is_scalar(q) && isequal(q.coeff, CNUM_ONE)

qexpr_arg(q::QSym) = to_qadd(q)
qexpr_arg(q::QExprArg) = q

function is_scalar_qadd(q::QAdd)
    isempty(q.indices) || return false
    length(q.arguments) == 1 || return false
    term = first(keys(q.arguments))
    return isempty(term.ops) && isempty(term.ne)
end
scalar_qadd_coeff(q::QAdd) = first(values(q.arguments))

function qexpr_from_qadd(q::QAdd)
    iszero(q) && return qexpr_zero()
    is_scalar_qadd(q) && return QExpr(QEXPR_MUL, scalar_qadd_coeff(q), QExprArg[])
    return QExpr(QEXPR_MUL, CNUM_ONE, QExprArg[q])
end

function polynomial_body(q::QExpr)
    (q.kind == QEXPR_MUL && length(q.args) == 1) || return (false, ZERO_QADD)
    body = only(q.args)
    body isa QAdd || return (false, ZERO_QADD)
    return (true, body)
end

function polynomial_qadd(q::QExpr)
    qexpr_is_scalar(q) && return (true, single_qadd(q.coeff, EMPTY_OPS))
    iszero(q) && return (true, zero_qadd())
    return polynomial_body(q)
end

qexpr_same_body(a::QExpr, b::QExpr) = a.kind == b.kind && isequal(a.args, b.args)

function qexpr_scale(q::QExpr, c::CNum)
    iszero(c) && return qexpr_zero()
    isequal(c, CNUM_ONE) && return q
    polynomial, body = polynomial_body(q)
    polynomial && return qexpr_from_qadd(body * c)
    if q.kind == QEXPR_ADD
        args = QExprArg[]
        sizehint!(args, length(q.args))
        for arg in q.args
            push!(args, arg isa QAdd ? arg * c : qexpr_scale(arg, c))
        end
        return qexpr_sum(args)
    end
    new_c = q.coeff * c
    iszero(new_c) && return qexpr_zero()
    return QExpr(q.kind, new_c, copy(q.args))
end

function qexpr_collect_sum!(formal::Vector{QExpr}, polynomial::QAddBuilder, arg::QAdd)
    iszero(arg) || accumulate!(polynomial, arg)
    return polynomial
end

function qexpr_collect_sum!(formal::Vector{QExpr}, polynomial::QAddBuilder, arg::QExpr)
    iszero(arg) && return polynomial
    if arg.kind == QEXPR_ADD
        for child in arg.args
            qexpr_collect_sum!(formal, polynomial, child)
        end
    else
        is_polynomial, body = polynomial_qadd(arg)
        is_polynomial ? accumulate!(polynomial, body) : push!(formal, arg)
    end
    return polynomial
end

function qexpr_coalesce_formal(formal::Vector{QExpr})
    length(formal) <= 1 && return formal
    if length(formal) == 2
        a, b = formal
        if qexpr_same_body(a, b)
            new_c = a.coeff + b.coeff
            return iszero(new_c) ? QExpr[] : QExpr[QExpr(a.kind, new_c, copy(a.args))]
        end
        qexpr_body_less(b, a) && ((formal[1], formal[2]) = (b, a))
        return formal
    end
    sort!(formal; lt = qexpr_body_less)
    n = 0
    for arg in formal
        if n > 0 && qexpr_same_body(formal[n], arg)
            new_c = formal[n].coeff + arg.coeff
            iszero(new_c) ? (n -= 1) : (formal[n] = QExpr(arg.kind, new_c, copy(arg.args)))
        else
            n += 1
            formal[n] = arg
        end
    end
    return resize!(formal, n)
end

function qexpr_sum(raw::Vector{QExprArg})
    polynomial = QAddBuilder()
    formal = QExpr[]
    sizehint!(formal, length(raw))
    for arg in raw
        qexpr_collect_sum!(formal, polynomial, arg)
    end

    formal = qexpr_coalesce_formal(formal)
    polynomial_value = build(polynomial)
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

"""
    monomial_content(lead, coeffs)

Leading scalar of `lead` times each parameter's lowest power across `coeffs`, so dividing
every coefficient by it stays exact. A raw symbolic lead has no canonical scalar and
yields `CNUM_ONE`, which extracts nothing.
"""
function monomial_content(lead::CNum, coeffs)
    is_native(lead) && return lead
    tail = lead.tail
    tail isa Poly || return CNUM_ONE
    m = first(tail.terms)
    syms = copy(m.syms)
    exps = copy(m.exps)
    for c in coeffs
        ctail = c.tail
        ctail isa Poly || (empty!(syms); break)
        for term in ctail.terms
            common_powers!(syms, exps, term)
        end
        isempty(syms) && break
    end
    scalar = scalar_coeff(monomial_scalar(m))
    isempty(syms) && return scalar
    return scalar * tiered(monomial_terms, EXACT_ONE, syms, exps)
end

function common_powers!(
        syms::Vector{RawExpression}, exps::Vector{Rational{Int}}, term::Monomial,
    )
    for k in reverse(eachindex(syms))
        j = position_of(t -> t === syms[k], term.syms)
        if iszero(j)
            deleteat!(syms, k)
            deleteat!(exps, k)
        else
            exps[k] = min(exps[k], term.exps[j])
        end
    end
    return syms
end

function depends_on_any_index(c::CNum, indices::Vector{Index})
    for idx in indices
        depends_on_index_term(c, EMPTY_OPS, idx) && return true
    end
    return false
end

function qexpr_split_content(q::QAdd)
    uniform = allequal(values(q.arguments))
    content = if uniform
        first(values(q.arguments))
    else
        lead = argmin(term_order_key, keys(q.arguments))
        monomial_content(q.arguments[lead], values(q.arguments))
    end
    isequal(content, CNUM_ONE) && return (CNUM_ONE, q)
    depends_on_any_index(content, q.indices) && return (CNUM_ONE, q)
    scale = uniform ? CNUM_ONE : inv(content)
    if !uniform && is_native(scale)
        z = native_scalar(scale)
        z isa ComplexF64 && !isfinite(z) && return (CNUM_ONE, q)
    end
    d = QTermDict()
    sizehint!(d, length(q.arguments))
    for (term, c) in q.arguments
        normalized = uniform ? CNUM_ONE : c * scale
        isequal(content * normalized, c) || return (CNUM_ONE, q)
        d[copy_key(term)] = normalized
    end
    return (content, QAdd(d, q.indices))
end

function split_term(term::QTerm, partition::Vector{Int})
    inside = filter(op -> op.space_index in partition, term.ops)
    outside = filter(op -> op.space_index ∉ partition, term.ops)
    return (QTerm(inside, EMPTY_NE), QTerm(outside, EMPTY_NE))
end

function qexpr_tensor_factors(q::QAdd)::Vector{QAdd}
    spaces = acts_on(q)
    length(spaces) > 1 || return QAdd[q]
    isempty(q.indices) || return QAdd[q]
    any(term -> !isempty(term.ne), keys(q.arguments)) && return QAdd[q]
    pivot = argmin(term_order_key, keys(q.arguments))
    pivot_coeff = q.arguments[pivot]
    for n in 1:div(length(spaces), 2), partition in combinations(spaces, n)
        left_key, right_key = split_term(pivot, partition)
        left = QTermDict()
        right = QTermDict()
        for (term, c) in q.arguments
            l, r = split_term(term, partition)
            if isequal(r, right_key)
                left[l] = isequal(c, pivot_coeff) ? CNUM_ONE : c / pivot_coeff
            end
            isequal(l, left_key) && (right[r] = c)
        end
        lpoly = QAdd(left, EMPTY_INDICES)
        rpoly = QAdd(right, EMPTY_INDICES)
        if isequal(lpoly * rpoly, q)
            return vcat(qexpr_tensor_factors(lpoly), qexpr_tensor_factors(rpoly))
        end
    end
    return QAdd[q]
end

function qexpr_push_product!(factors::Vector{QExprArg}, coeff::CNum, arg::QAdd)
    iszero(arg) && return (CNUM_ZERO, true)
    if is_scalar_qadd(arg)
        coeff = coeff * scalar_qadd_coeff(arg)
        return (coeff, iszero(coeff))
    end
    parts = qexpr_tensor_factors(arg)
    if length(parts) > 1
        for part in parts
            coeff, stopped = qexpr_push_product!(factors, coeff, part)
            stopped && return (CNUM_ZERO, true)
        end
        return (coeff, false)
    end
    left = isempty(factors) ? arg : last(factors)
    if !isempty(factors) && left isa QAdd && !isdisjoint(acts_on(left), acts_on(arg))
        pop!(factors)
        merged = left * arg
        iszero(merged) && return (CNUM_ZERO, true)
        is_scalar_qadd(merged) || return qexpr_push_product!(factors, coeff, merged)
        coeff = coeff * scalar_qadd_coeff(merged)
    else
        content, normalized = qexpr_split_content(arg)
        coeff = coeff * content
        push!(factors, normalized)
    end
    return (coeff, iszero(coeff))
end

function qexpr_push_product!(factors::Vector{QExprArg}, coeff::CNum, arg::QExpr)
    iszero(arg) && return (CNUM_ZERO, true)
    coeff = coeff * arg.coeff
    iszero(coeff) && return (CNUM_ZERO, true)
    if arg.kind == QEXPR_MUL
        for factor in arg.args
            coeff, stopped = qexpr_push_product!(factors, coeff, factor)
            stopped && return (CNUM_ZERO, true)
        end
        return (coeff, false)
    end
    bare = isequal(arg.coeff, CNUM_ONE) ? arg :
        QExpr(arg.kind, CNUM_ONE, copy(arg.args))
    push!(factors, bare)
    return (coeff, false)
end

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

function qexpr_commuting_order(factors::Vector{QExprArg}, coeff::CNum)
    while length(factors) > 1
        order = qexpr_lex_order(factors)
        issorted(order) && break
        reordered = QExprArg[]
        sizehint!(reordered, length(factors))
        for i in order
            coeff, stopped = qexpr_push_reordered!(reordered, coeff, factors[i])
            stopped && return (factors, CNUM_ZERO, true)
        end
        factors = reordered
    end
    return (factors, coeff, false)
end

function qexpr_push_reordered!(factors::Vector{QExprArg}, coeff::CNum, arg::QExprArg)
    if arg isa QAdd && !isempty(factors) && last(factors) isa QAdd &&
            !isdisjoint(acts_on(last(factors)), acts_on(arg))
        return qexpr_push_product!(factors, coeff, arg)
    end
    push!(factors, arg)
    return (coeff, false)
end

function qexpr_product(raw::Vector{QExprArg}, coefficient::CNum = CNUM_ONE)
    factors = QExprArg[]
    sizehint!(factors, length(raw))
    coeff = coefficient
    for arg in raw
        coeff, stopped = qexpr_push_product!(factors, coeff, arg)
        stopped && return qexpr_zero()
    end
    iszero(coeff) && return qexpr_zero()
    factors, coeff, stopped = qexpr_commuting_order(factors, coeff)
    stopped && return qexpr_zero()

    isempty(factors) && return QExpr(QEXPR_MUL, coeff, QExprArg[])
    if all(factor -> factor isa QAdd, factors)
        polynomial = single_qadd(coeff, EMPTY_OPS)
        for factor in factors
            factor isa QAdd && (polynomial = polynomial * factor)
        end
        return qexpr_from_qadd(polynomial)
    end
    if length(factors) == 1
        factor = only(factors)
        factor isa QExpr && return qexpr_scale(factor, coeff)
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

function qexpr_rebuild(kind::QExprKind, coeff::CNum, args::Vector{QExprArg})
    if kind == QEXPR_ADD
        return qexpr_scale(qexpr_sum(args), coeff)
    elseif kind == QEXPR_MUL
        return qexpr_product(args, coeff)
    end
    return qexpr_scale(qexpr_call(kind, only(args)), coeff)
end

"""
    sin(A)

Sine of an operator expression, kept as a formal function.

No series is expanded. The function enters sums and products like any other operator, and
keeps its position relative to operators on the same subsystem. `A` may be any operator
expression. Use [`taylor`](@ref) to turn the function into a polynomial.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> sin(a + a') * a
sin(a + a') * a

julia> taylor(sin(a + a'), 0:3)
1//2 * a + 1//2 * a' - 1//6 * a * a * a - 1//2 * a' * a * a - 1//2 * a' * a' * a - 1//6 * a' * a' * a'
```

See also [`cos`](@ref cos(::QField)), [`expim`](@ref expim(::QField)), [`taylor`](@ref).
"""
Base.sin(q::QField) = qexpr_call(QEXPR_SIN, qexpr_arg(q))

"""
    cos(A)

Cosine of an operator expression, kept as a formal function.

No series is expanded. The function enters sums and products like any other operator, and
keeps its position relative to operators on the same subsystem; a commutator with it stays
formal. `A` may be any operator expression. Use [`taylor`](@ref) to turn the function into a
polynomial.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> commutator(cos(a + a'), a)
-a * cos(a + a') + cos(a + a') * a

julia> taylor(cos(a + a'), 0:2)
1//2 - 1//2 * a * a - a' * a - 1//2 * a' * a'
```

See also [`sin`](@ref sin(::QField)), [`expim`](@ref expim(::QField)), [`taylor`](@ref).
"""
Base.cos(q::QField) = qexpr_call(QEXPR_COS, qexpr_arg(q))

qexpr_provably_hermitian(q::Union{QSym, QAdd}) = iszero(simplify(q - adjoint(q)))
qexpr_provably_hermitian(q::QExpr) = isequal(q, adjoint(q))

"""
    expim(A)

Exponential of `im` times a Hermitian operator expression, kept as a formal function.

The result is the unitary ``e^{iA}``. `A` must be Hermitian as written, meaning that
`A - A'` simplifies to zero; parameters in it must therefore be real. Otherwise an
`ArgumentError` is thrown. No series is expanded until [`taylor`](@ref) is called, and
products such as ``e^{iA} e^{-iA}`` are not combined. To change the frame of an expression
by such a unitary, use a [`UnitaryTransform`](@ref), which applies ``e^{iA} X e^{-iA}``
exactly.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables θ;

julia> U = expim(θ * a' * a)
exp(im*(θ * a' * a))

julia> taylor(U, 0:1)
1 + θim * a' * a

julia> expim(a)
ERROR: ArgumentError: `expim(::QField)` requires a provably Hermitian argument
[...]
```

See also [`taylor`](@ref), [`UnitaryTransform`](@ref), [`expim(::Real)`](@ref expim(::Real)).
"""
function expim(q::QField)
    qexpr_provably_hermitian(q) || throw(
        ArgumentError("`expim(::QField)` requires a provably Hermitian argument"),
    )
    return qexpr_call(QEXPR_EXPIM, qexpr_arg(q))
end

Base.:+(a::QExpr) = a
Base.:+(a::QExpr, b::Union{QExpr, QAdd}) = qexpr_sum(QExprArg[a, b])
Base.:+(a::QAdd, b::QExpr) = qexpr_sum(QExprArg[a, b])
Base.:+(a::QExpr, b::QSym) = qexpr_sum(QExprArg[a, to_qadd(b)])
Base.:+(a::QSym, b::QExpr) = qexpr_sum(QExprArg[to_qadd(a), b])
Base.:+(a::QExpr, b::Coefficient) = qexpr_sum(QExprArg[a, single_qadd(to_cnum(b), EMPTY_OPS)])
Base.:+(a::Coefficient, b::QExpr) = b + a

Base.:-(a::QExpr) = qexpr_scale(a, CNUM_NEG1)

Base.:*(a::QExpr, b::Union{QExpr, QAdd}) = qexpr_product(QExprArg[a, b])
Base.:*(a::QAdd, b::QExpr) = qexpr_product(QExprArg[a, b])
Base.:*(a::QExpr, b::QSym) = qexpr_product(QExprArg[a, to_qadd(b)])
Base.:*(a::QSym, b::QExpr) = qexpr_product(QExprArg[to_qadd(a), b])
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
    return a.kind == b.kind && isequal(a.coeff, b.coeff) && isequal(a.args, b.args)
end
Base.:(==)(a::QExpr, b::QExpr) = isequal(a, b)
Base.hash(q::QExpr, h::UInt) = hash(q.args, hash(q.coeff, hash(q.kind, hash(:QExpr, h))))

function qexpr_arg_less(a::QExprArg, b::QExprArg)
    if a isa QAdd
        b isa QAdd || return true
        return !isequal(a, b) && isless(a, b)
    end
    b isa QAdd && return false
    return !isequal(a, b) && isless(a, b)
end

function qexpr_args_less(a::Vector{QExprArg}, b::Vector{QExprArg})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        isequal(a[i], b[i]) && continue
        return qexpr_arg_less(a[i], b[i])
    end
    return length(a) < length(b)
end

qexpr_body_less(a::QExpr, b::QExpr) =
    a.kind == b.kind ? qexpr_args_less(a.args, b.args) : a.kind < b.kind

function Base.isless(a::QExpr, b::QExpr)
    a.kind != b.kind && return a.kind < b.kind
    isequal(a.args, b.args) || return qexpr_args_less(a.args, b.args)
    return isless(coeff_key(a.coeff), coeff_key(b.coeff))
end

function Base.adjoint(q::QExpr)
    q.kind == QEXPR_EXPIM && return map_args(arg -> -adjoint(arg), q, conj(q.coeff))
    q.kind == QEXPR_MUL || return map_args(adjoint, q, conj(q.coeff))
    return qexpr_product(QExprArg[adjoint(arg) for arg in Iterators.reverse(q.args)], conj(q.coeff))
end
