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
    isempty(vars) || throw(
        ArgumentError(
            "numeric conversion requires static QExpr with no free parameters; " *
                "supply values for $(join(vars, ", ")) via `parameter`.",
        )
    )
    return static_operator(q_sub, ctx, op_type)
end

to_numeric(::QExpr, ::Basis, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::Integer, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::AbstractVector{<:Integer}, ::AbstractDict{<:QSym}) =
    qexpr_numeric_error()
to_numeric(::QExpr, ::Tuple, ::AbstractDict{<:QSym}) = qexpr_numeric_error()

function static_operator(q::QExpr, ctx::NumericContext, op_type)
    if q.kind == QEXPR_ADD || q.kind == QEXPR_MUL
        if isempty(q.args)
            result = numeric_identity(ctx.backend, ctx.basis)
            q.kind == QEXPR_ADD && (result = 0.0im * result)
        else
            result = static_operator(first(q.args), ctx, EagerOperator())
            for arg in Iterators.drop(q.args, 1)
                child = static_operator(arg, ctx, EagerOperator())
                result = q.kind == QEXPR_ADD ? result + child : result * child
            end
        end
    else
        argument = static_operator(only(q.args), ctx, EagerOperator())
        f = q.kind == QEXPR_SIN ? sin :
            q.kind == QEXPR_COS ? cos : A -> exp(im * A)
        result = numeric_function(ctx.backend, f, argument)
    end
    return numeric_materialize(ctx.backend, to_complex(q.coeff) * result, op_type)
end
