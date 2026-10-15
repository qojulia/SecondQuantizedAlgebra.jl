@noinline function nested_taylor_error()
    throw(
        ArgumentError(
            "Taylor lowering of nested formal operator functions is not supported; " *
                "lower the inner function explicitly first",
        ),
    )
end

@noinline function bound_sum_taylor_error()
    throw(
        ArgumentError(
            "Taylor lowering of a formal function with a bound summation scope is not " *
                "supported; powers require fresh dummy indices",
        ),
    )
end

function taylor_prefix_order(ns::UnitRange{T})::Int where {T <: Integer}
    isempty(ns) && throw(ArgumentError("operator Taylor range must be non-empty"))
    first(ns) == 0 || throw(
        ArgumentError("operator Taylor range must start at 0, got $ns"),
    )
    last(ns) <= typemax(Int) || throw(
        ArgumentError("Taylor order is too large for this platform: $(last(ns))"),
    )
    return Int(last(ns))
end

function taylor_inverse_denominator(denominator::BigInt)::CNum
    denominator <= typemax(Int) && return to_cnum(1 // Int(denominator))
    return to_cnum(big(1) // denominator)
end

@inline taylor_scale(q::QAdd, c::CNum) = isequal(c, CNUM_ONE) ? q : q * c

lower_polynomial_arg(q::QAdd) = q

function lower_polynomial_arg(q::QExpr)::QAdd
    if q.kind == QEXPR_ADD
        builder = QAddBuilder()
        for arg in q.args
            accumulate!(builder, lower_polynomial_arg(arg))
        end
        return build(builder)
    elseif q.kind == QEXPR_MUL
        isempty(q.args) && return single_qadd(q.coeff, EMPTY_OPS)
        result = taylor_scale(lower_polynomial_arg(first(q.args)), q.coeff)
        @inbounds for i in 2:length(q.args)
            result = result * lower_polynomial_arg(q.args[i])
        end
        return result
    end
    nested_taylor_error()
end

function taylor_function_argument(arg::QExprArg)::QAdd
    polynomial = lower_polynomial_arg(arg)
    isempty(polynomial.indices) || bound_sum_taylor_error()
    return polynomial
end

function taylor_cos(argument::QAdd, order::Int)::QAdd
    builder = QAddBuilder()
    accumulate!(builder, CNUM_ONE)
    order < 2 && return build(builder)

    denominator = big(1)
    square = argument * argument
    power = square
    for n in 2:2:order
        denominator *= n - 1
        denominator *= n
        coefficient = taylor_inverse_denominator(denominator)
        isodd(n ÷ 2) && (coefficient = neg_cnum(coefficient))
        accumulate!(builder, power, coefficient)
        n + 2 <= order && (power = power * square)
    end
    return build(builder)
end

function taylor_sin(argument::QAdd, order::Int)::QAdd
    order < 1 && return zero_qadd()
    builder = QAddBuilder()
    accumulate!(builder, argument)
    order < 3 && return build(builder)

    denominator = big(1)
    square = argument * argument
    power = argument * square
    for n in 3:2:order
        denominator *= n - 1
        denominator *= n
        coefficient = taylor_inverse_denominator(denominator)
        isodd((n - 1) ÷ 2) && (coefficient = neg_cnum(coefficient))
        accumulate!(builder, power, coefficient)
        n + 2 <= order && (power = power * square)
    end
    return build(builder)
end

function taylor_expim(argument::QAdd, order::Int)::QAdd
    builder = QAddBuilder()
    accumulate!(builder, CNUM_ONE)
    order < 1 && return build(builder)

    denominator = big(1)
    power = argument
    for n in 1:order
        denominator *= n
        coefficient = taylor_inverse_denominator(denominator)
        phase = mod(n, 4)
        if phase == 1
            coefficient = mul_cnum(coefficient, CNUM_IM)
        elseif phase == 2
            coefficient = neg_cnum(coefficient)
        elseif phase == 3
            coefficient = mul_cnum(coefficient, CNUM_NEG_IM)
        end
        accumulate!(builder, power, coefficient)
        n < order && (power = power * argument)
    end
    return build(builder)
end

function taylor_function(kind::QExprKind, argument::QAdd, order::Int)::QAdd
    kind == QEXPR_COS && return taylor_cos(argument, order)
    kind == QEXPR_SIN && return taylor_sin(argument, order)
    return taylor_expim(argument, order)
end

lower_taylor_arg(q::QAdd, ::Int, ::Dict{QExpr, QAdd}) = q
lower_taylor_arg(q::QExpr, order::Int, cache::Dict{QExpr, QAdd}) =
    lower_taylor_qexpr(q, order, cache)

function lower_taylor_qexpr(q::QExpr, order::Int, cache::Dict{QExpr, QAdd})::QAdd
    haskey(cache, q) && return cache[q]

    result = if q.kind == QEXPR_ADD
        builder = QAddBuilder()
        for arg in q.args
            accumulate!(builder, lower_taylor_arg(arg, order, cache))
        end
        build(builder)
    elseif q.kind == QEXPR_MUL
        if isempty(q.args)
            single_qadd(q.coeff, EMPTY_OPS)
        else
            product = taylor_scale(lower_taylor_arg(first(q.args), order, cache), q.coeff)
            @inbounds for i in 2:length(q.args)
                product = product * lower_taylor_arg(q.args[i], order, cache)
            end
            product
        end
    else
        argument = taylor_function_argument(only(q.args))
        taylor_scale(taylor_function(q.kind, argument, order), q.coeff)
    end

    cache[q] = result
    return result
end

function lower_taylor_qexpr(q::QExpr, order::Int)::QAdd
    if q.kind == QEXPR_SIN || q.kind == QEXPR_COS || q.kind == QEXPR_EXPIM
        argument = taylor_function_argument(only(q.args))
        return taylor_scale(taylor_function(q.kind, argument, order), q.coeff)
    end
    return lower_taylor_qexpr(q, order, Dict{QExpr, QAdd}())
end

"""
    taylor(expr, 0:n)

Replace every formal operator function in `expr` by its Taylor polynomial up to order `n`.

The coefficients are exact and the result is a polynomial in canonical form; an expression
without formal functions is returned unchanged. Each function
is truncated on its own, so the order is not a total degree of the result: `sin(A)^2`
truncated at order `1` is `A^2`, and `sin(A)^2 + cos(A)^2` truncated at any order is not
`1`. A function inside another function throws an `ArgumentError`; lower the inner part
first. So does a function of an argument containing a [`Σ`](@ref), since a power such as
``(\\sum_i A_i)^2`` would need fresh summation indices. The range must start at `0`.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> A = a + a';

julia> taylor(sin(A)^2, 0:1)
1 + a * a + 2 * a' * a + a' * a'

julia> taylor(sin(A)^2, 0:1) == A^2
true
```

See also [`sin`](@ref sin(::QField)), [`cos`](@ref cos(::QField)), [`expim`](@ref expim(::QField)).
"""
function Symbolics.taylor(q::QExpr, ns::UnitRange{T})::QAdd where {T <: Integer}
    return lower_taylor_qexpr(q, taylor_prefix_order(ns))
end
function Symbolics.taylor(q::QAdd, ns::UnitRange{T})::QAdd where {T <: Integer}
    taylor_prefix_order(ns)
    return q
end
Symbolics.taylor(q::QSym, ns::UnitRange{T}) where {T <: Integer} = taylor(+q, ns)
