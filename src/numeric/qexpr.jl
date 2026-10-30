@noinline qexpr_numeric_error() = throw(
    ArgumentError(
        "this operation does not support formal operator functions; " *
            "call `taylor(expr, 0:n)` first and convert the resulting `QAdd`",
    ),
)

to_numeric_lazy(::QExpr, ::NumericContext) = qexpr_numeric_error()

function to_numeric_translated(
        q::QExpr, ctx::NumericContext, parameter, time_parameter, op_type,
    )
    isempty(time_parameter) || throw(
        ArgumentError(
            "numeric conversion supports static QExpr only; nonempty time_parameter is unsupported, " *
                "including scalar time-dependent prefactors. Supply static values via parameter.",
        )
    )
    q_sub = substitute_with(q, parameter)
    vars = get_variables(q_sub)
    (has_numeric_sum(q_sub) || isempty(vars)) || throw(
        ArgumentError(
            "numeric conversion requires static QExpr with no free parameters; " *
                "supply values for $(join(vars, ", ")) via `parameter`.",
        )
    )
    return formal_static_operator(q_sub, ctx, op_type, parameter)
end

to_numeric(::QExpr, ::Basis, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::Integer, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::AbstractVector{<:Integer}, ::AbstractDict{<:QSym}) =
    qexpr_numeric_error()
to_numeric(::QExpr, ::Tuple, ::AbstractDict{<:QSym}) = qexpr_numeric_error()

has_numeric_sum(q::QAdd) = !isempty(q.indices)
has_numeric_sum(q::QExpr) = any(has_numeric_sum, q.args)

static_operator(q::QExpr, ctx::NumericContext, op_type) =
    formal_static_operator(q, ctx, op_type, SubstitutionRules())

function formal_numeric_argument(q::QAdd, ctx::NumericContext, rules::SubstitutionRules)
    isempty(q.indices) && return to_numeric_translated(q, ctx, rules, NO_TIME_PARAMETER, EagerOperator())
    return indexed_result(q, ctx, rules.scalars, NO_TIME_PARAMETER, EagerOperator())
end
formal_numeric_argument(q::QExpr, ctx::NumericContext, rules::SubstitutionRules) =
    formal_static_operator(q, ctx, EagerOperator(), rules)

function formal_static_operator(q::QExpr, ctx::NumericContext, op_type, rules::SubstitutionRules)
    if q.kind == QEXPR_ADD || q.kind == QEXPR_MUL
        if isempty(q.args)
            result = numeric_identity(ctx.backend, ctx.basis)
            q.kind == QEXPR_ADD && (result = 0.0im * result)
        else
            result = formal_numeric_argument(first(q.args), ctx, rules)
            for arg in Iterators.drop(q.args, 1)
                child = formal_numeric_argument(arg, ctx, rules)
                result = q.kind == QEXPR_ADD ? result + child : result * child
            end
        end
    else
        argument = formal_numeric_argument(only(q.args), ctx, rules)
        result = if q.kind == QEXPR_SIN
            numeric_function(ctx.backend, sin, argument)
        elseif q.kind == QEXPR_COS
            numeric_function(ctx.backend, cos, argument)
        else
            numeric_function(ctx.backend, expim, argument)
        end
    end
    return numeric_materialize(ctx.backend, to_complex(q.coeff) * result, op_type)
end
