# Explicit Maclaurin lowering for the cold-path formal expression layer.

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

function taylor_inverse_factorial(n::Int)::CNum
    denominator = factorial(big(n))
    if denominator <= typemax(Int)
        return to_cnum(1 // Int(denominator))
    end
    return to_cnum(big(1) // denominator)
end

lower_polynomial_arg(q::QAdd) = q

function lower_polynomial_arg(q::QExpr)::QAdd
    if q.kind == QEXPR_ADD
        result = zero_qadd()
        for arg in q.args
            result = result + lower_polynomial_arg(arg)
        end
        return result * q.coeff
    elseif q.kind == QEXPR_MUL
        result = single_qadd(q.coeff, EMPTY_OPS)
        for arg in q.args
            result = result * lower_polynomial_arg(arg)
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
    result = single_qadd(CNUM_ONE, EMPTY_OPS)
    order < 2 && return result

    square = argument * argument
    power = square
    for n in 2:2:order
        coefficient = taylor_inverse_factorial(n)
        isodd(n ÷ 2) && (coefficient = neg_cnum(coefficient))
        result = result + power * coefficient
        n + 2 <= order && (power = power * square)
    end
    return result
end

function taylor_sin(argument::QAdd, order::Int)::QAdd
    order < 1 && return zero_qadd()
    result = argument
    order < 3 && return result

    square = argument * argument
    power = argument * square
    for n in 3:2:order
        coefficient = taylor_inverse_factorial(n)
        isodd((n - 1) ÷ 2) && (coefficient = neg_cnum(coefficient))
        result = result + power * coefficient
        n + 2 <= order && (power = power * square)
    end
    return result
end

function taylor_expim(argument::QAdd, order::Int)::QAdd
    result = single_qadd(CNUM_ONE, EMPTY_OPS)
    order < 1 && return result

    power = argument
    for n in 1:order
        coefficient = taylor_inverse_factorial(n)
        phase = mod(n, 4)
        if phase == 1
            coefficient = mul_cnum(coefficient, CNUM_IM)
        elseif phase == 2
            coefficient = neg_cnum(coefficient)
        elseif phase == 3
            coefficient = mul_cnum(coefficient, CNUM_NEG_IM)
        end
        result = result + power * coefficient
        n < order && (power = power * argument)
    end
    return result
end

function taylor_function(kind::QExprKind, argument::QAdd, order::Int)::QAdd
    kind == QEXPR_COS && return taylor_cos(argument, order)
    kind == QEXPR_SIN && return taylor_sin(argument, order)
    kind == QEXPR_EXPIM && return taylor_expim(argument, order)
    error("Taylor lowering requested for non-function QExpr kind $kind")
end

lower_taylor_arg(q::QAdd, ::Int) = q
lower_taylor_arg(q::QExpr, order::Int) = lower_taylor_qexpr(q, order)

function lower_taylor_qexpr(q::QExpr, order::Int)::QAdd
    if q.kind == QEXPR_ADD
        result = zero_qadd()
        for arg in q.args
            result = result + lower_taylor_arg(arg, order)
        end
        return result * q.coeff
    elseif q.kind == QEXPR_MUL
        result = single_qadd(q.coeff, EMPTY_OPS)
        for arg in q.args
            result = result * lower_taylor_arg(arg, order)
        end
        return result
    end

    argument = taylor_function_argument(only(q.args))
    return taylor_function(q.kind, argument, order) * q.coeff
end

"""
    taylor(expr::QExpr, ns::UnitRange{<:Integer}) -> QAdd

Explicitly lower formal `sin`, `cos`, and `expim` nodes to Maclaurin polynomials.
For the first operator-function API, `ns` must be a prefix range `0:n`.
Coefficients remain exact and generated products use the ordinary canonical `QAdd`
pipeline.

Each formal function must have a polynomial argument with no bound `QAdd` summation
scope. Nested formal functions and powers of bound sums are rejected explicitly rather
than assigned ambiguous implicit approximation or alpha-renaming semantics.
"""
function Symbolics.taylor(q::QExpr, ns::UnitRange{T})::QAdd where {T <: Integer}
    return lower_taylor_qexpr(q, taylor_prefix_order(ns))
end
