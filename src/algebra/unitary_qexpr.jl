conjugate(q::QExpr, U::UnitaryTransform)::QExpr = map_args(
    arg -> conjugate(arg, U), q, reduce_affine(q.coeff, U.action.relations, ParamRelation[]),
)

transform(q::QExpr, U::UnitaryTransform{StaticTime}) = conjugate(q, U)
transform(q::QExpr, U::UnitaryTransform{DynamicTime}) = conjugate(q, U) + U.gauge
