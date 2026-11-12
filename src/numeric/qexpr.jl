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
    q_sub = substitute_with(q, parameter)
    if !isempty(time_parameter)
        check_td_op_type(op_type)
        terms = formal_td_terms(q_sub, ctx, parameter)
        return formal_td_operator(terms, ctx, time_parameter)
    end
    vars = get_variables(q_sub)
    (has_numeric_sum(q_sub) || isempty(vars)) || throw(
        ArgumentError(
            "numeric conversion requires static QExpr with no free parameters; " *
                "supply values for $(join(vars, ", ")) via `parameter`.",
        )
    )
    return formal_static_operator(q_sub, ctx, op_type, parameter)
end

function formal_td_operator(terms::Vector{Tuple{CNum, O}}, ctx::NumericContext, time_parameter) where {O}
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    td_terms = Tuple{TDCoeff, Vector{O}}[]
    for (c, op) in terms
        push!(td_terms, (td_coeff(to_num(c), basevars, valuefuncs, p_aware), [op]))
    end
    return numeric_assemble_td(ctx.backend, ctx.basis, td_terms)
end

function formal_td_terms(q::QAdd, ctx::NumericContext, rules::SubstitutionRules)
    if isempty(q.indices)
        terms = coeff_terms(q)
    else
        ctx = sum_context(q, ctx, rules.scalars)
        validate_sites(ctx.backend, ctx.basis, ctx.sites)
        terms = unrolled_terms(q, ctx.sites, split_scalar_subs(rules.scalars))
    end
    identity_op = numeric_materialize(ctx.backend, numeric_identity(ctx.backend, ctx.basis), EagerOperator())
    result = Tuple{CNum, typeof(identity_op)}[]
    for (c, ops) in terms
        op = static_operator(single_qadd(CNUM_ONE, ops), ctx, EagerOperator())
        push!(result, (c, op))
    end
    isempty(result) && push!(result, (CNUM_ZERO, identity_op))
    return result
end

function formal_td_terms(q::QExpr, ctx::NumericContext, rules::SubstitutionRules)
    identity_op = numeric_materialize(ctx.backend, numeric_identity(ctx.backend, ctx.basis), EagerOperator())
    result = Tuple{CNum, typeof(identity_op)}[]
    if q.kind == QEXPR_ADD
        for arg in q.args
            append!(result, formal_td_terms(arg, ctx, rules))
        end
        isempty(result) && push!(result, (CNUM_ZERO, identity_op))
    elseif q.kind == QEXPR_MUL
        push!(result, (q.coeff, identity_op))
        for arg in q.args
            children = formal_td_terms(arg, ctx, rules)
            result = [(c * d, left * right) for (c, left) in result for (d, right) in children]
        end
    else
        # Only coefficients outside matrix functions may vary during a solve.
        argument = formal_numeric_argument(only(q.args), ctx, rules)
        op = numeric_materialize(ctx.backend, formal_function(q.kind, ctx, argument), EagerOperator())
        push!(result, (q.coeff, op))
    end
    return result
end

to_numeric(q::QExpr, b::Basis, d::AbstractDict{<:QSym}) =
    to_numeric(q, b; operators = d, adjoint_ops = false)

has_numeric_sum(q::QAdd) = !isempty(q.indices)
has_numeric_sum(q::QExpr) = any(has_numeric_sum, q.args)

static_operator(q::QExpr, ctx::NumericContext, op_type) =
    formal_static_operator(q, ctx, op_type, SubstitutionRules())

function formal_numeric_argument(q::QAdd, ctx::NumericContext, rules::SubstitutionRules)
    isempty(q.indices) && return static_operator(apply_parameter(q, rules), ctx, EagerOperator())
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
        result = formal_function(q.kind, ctx, formal_numeric_argument(only(q.args), ctx, rules))
    end
    return numeric_materialize(ctx.backend, to_complex(q.coeff) * result, op_type)
end

# One branch per function so each `numeric_function` call specializes on a concrete `f`.
function formal_function(kind::QExprKind, ctx::NumericContext, argument)
    kind == QEXPR_SIN && return numeric_function(ctx.backend, sin, argument)
    kind == QEXPR_COS && return numeric_function(ctx.backend, cos, argument)
    return numeric_function(ctx.backend, expim, argument)
end
