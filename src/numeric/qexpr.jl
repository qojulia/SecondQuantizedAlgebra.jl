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

to_numeric(::QExpr, ::Basis, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::Integer, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
to_numeric(::QExpr, ::AbstractVector{<:Integer}, ::AbstractDict{<:QSym}) =
    qexpr_numeric_error()
to_numeric(::QExpr, ::Tuple, ::AbstractDict{<:QSym}) = qexpr_numeric_error()
