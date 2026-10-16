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

taylor_scale(q::QAdd, c::CNum) = isequal(c, CNUM_ONE) ? q : q * c

lower_qexpr(::F, q::QAdd) where {F} = q

function lower_qexpr(lower_function::F, q::QExpr)::QAdd where {F}
    if q.kind == QEXPR_ADD
        builder = QAddBuilder()
        for arg in q.args
            accumulate!(builder, lower_qexpr(lower_function, arg))
        end
        return build(builder)
    elseif q.kind == QEXPR_MUL
        isempty(q.args) && return single_qadd(q.coeff, EMPTY_OPS)
        result = taylor_scale(lower_qexpr(lower_function, first(q.args)), q.coeff)
        @inbounds for i in 2:length(q.args)
            result = result * lower_qexpr(lower_function, q.args[i])
        end
        return result
    end
    return lower_function(q)
end

function taylor_function_argument(arg::QExprArg)::QAdd
    polynomial = lower_qexpr(f -> nested_taylor_error(), arg)
    isempty(polynomial.indices) || bound_sum_taylor_error()
    return polynomial
end

function taylor_function(kind::QExprKind, argument::QAdd, order::Int)::QAdd
    builder = QAddBuilder()
    stride = kind == QEXPR_EXPIM ? 1 : 2
    lowest = kind == QEXPR_SIN ? 1 : 0
    lowest > order && return build(builder)
    step = stride == 1 || order < lowest + stride ? argument : argument * argument
    power = lowest == 1 ? argument : single_qadd(CNUM_ONE, EMPTY_OPS)
    denominator = big(1)
    for n in lowest:stride:order
        if n > lowest
            denominator *= stride == 1 ? n : n * (n - 1)
            power = n == stride ? step : power * step
        end
        coefficient = taylor_inverse_denominator(denominator)
        quarter = mod(n, 4)
        if quarter == 2 || (quarter == 3 && kind != QEXPR_EXPIM)
            coefficient = -coefficient
        elseif kind == QEXPR_EXPIM && isodd(quarter)
            coefficient = coefficient * (quarter == 1 ? CNUM_IM : CNUM_NEG_IM)
        end
        accumulate!(builder, power, coefficient)
    end
    return build(builder)
end

function lower_taylor_qexpr(q::QExpr, order::Int)::QAdd
    cache = Dict{QExpr, QAdd}()
    return lower_qexpr(q) do f
        get!(cache, f) do
            argument = taylor_function_argument(only(f.args))
            taylor_scale(taylor_function(f.kind, argument, order), f.coeff)
        end
    end
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
function Symbolics.taylor(q::Union{QSym, QAdd}, ns::UnitRange{T})::QAdd where {T <: Integer}
    taylor_prefix_order(ns)
    return to_qadd(q)
end
