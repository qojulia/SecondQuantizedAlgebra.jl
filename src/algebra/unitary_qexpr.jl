# Exact unitary transformations for the cold-path formal expression layer.

qexpr_conjugate_arg(arg::QAdd, U::UnitaryTransform) = conjugate(arg, U)
qexpr_conjugate_arg(arg::QExpr, U::UnitaryTransform) = conjugate(arg, U)

function qexpr_reduce_unitary_coeff(c::CNum, U::UnitaryTransform)::CNum
    relations = U.action.relations
    isempty(relations) && return c
    return reduce_all(c, relations, true, ParamRelation[])
end

function conjugate(q::QExpr, U::UnitaryTransform)::QExpr
    args = QExprArg[]
    sizehint!(args, length(q.args))
    for arg in q.args
        push!(args, qexpr_conjugate_arg(arg, U))
    end
    coeff = qexpr_reduce_unitary_coeff(q.coeff, U)
    return qexpr_rebuild(q.kind, coeff, args)
end

transform(q::QExpr, U::UnitaryTransform{StaticTime}) = conjugate(q, U)
transform(q::QExpr, U::UnitaryTransform{DynamicTime}) = conjugate(q, U) + U.gauge
