# Deliberate numeric boundary for formal operator functions.

@noinline qexpr_numeric_error() = throw(
    ArgumentError(
        "direct numeric conversion of formal operator functions is not implemented; " *
            "call `taylor(expr, 0:n)` first and convert the resulting `QAdd`",
    ),
)

to_numeric_lazy(::QExpr, ::NumericContext) = qexpr_numeric_error()

to_numeric_translated(
    ::QExpr, ::NumericContext, parameter, time_parameter, op_type,
) = qexpr_numeric_error()

# Positional numeric routes bypass `to_numeric_translated`. Keep these strictly
# three-argument so the existing two-argument keyword methods remain the selected surface.
to_numeric(::QExpr, ::Basis, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::Integer, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::AbstractVector{<:Integer}, ::AbstractDict{<:QSym}) =
    qexpr_numeric_error()
# QExpr is rejected before backend dimension validation, so no element type parameter is
# needed here. Keeping the tuple unconstrained also avoids an Aqua unbound-typevar report.
to_numeric(::QExpr, ::Tuple, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
