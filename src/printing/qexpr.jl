const QEXPR_PREC_ADD = 1
const QEXPR_PREC_MUL = 2
const QEXPR_PREC_FUNC = 3

struct TextStyle end
struct LatexStyle end

function qexpr_precedence(q::QExpr)
    q.kind == QEXPR_ADD && return QEXPR_PREC_ADD
    (q.kind == QEXPR_MUL || !is_unit(q.coeff)) && return QEXPR_PREC_MUL
    return QEXPR_PREC_FUNC
end

qexpr_qadd_needs_group(q::QAdd) = !isempty(q.indices) || length(q.arguments) > 1

group_delimiters(::TextStyle) = ("(", ")")
group_delimiters(::LatexStyle) = ("\\left( ", " \\right)")

product_separator(::TextStyle) = " * "
product_separator(::LatexStyle) = " "

function function_delimiters(::TextStyle, kind::QExprKind)
    kind == QEXPR_SIN && return ("sin(", ")")
    kind == QEXPR_COS && return ("cos(", ")")
    return ("exp(im*(", "))")
end

function function_delimiters(::LatexStyle, kind::QExprKind)
    kind == QEXPR_SIN && return ("\\sin\\left( ", " \\right)")
    kind == QEXPR_COS && return ("\\cos\\left( ", " \\right)")
    return ("e^{i\\left( ", " \\right)}")
end

show_polynomial(io::IO, ::TextStyle, q::QAdd) = show(io, q)
show_polynomial(io::IO, ::LatexStyle, q::QAdd) = write(io, latex_inline(q))

show_scalar(io::IO, ::TextStyle, c::CNum) = show_prefactor(io, c)
show_scalar(io::IO, ::LatexStyle, c::CNum) = write(io, latex_inline(c))

show_qexpr_prefactor(io::IO, ::TextStyle, c::CNum) = show_scaled_prefactor(io, c)

function show_qexpr_prefactor(io::IO, ::LatexStyle, c::CNum)
    is_neg_unit(c) && return write(io, "-")
    is_unit(c) && return
    write(io, latex_inline(c), " ")
    return
end

function show_qexpr(io::IO, style::S, q::QAdd, parent_precedence::Int) where {S}
    group = parent_precedence > QEXPR_PREC_ADD && qexpr_qadd_needs_group(q)
    open, close = group_delimiters(style)
    group && write(io, open)
    show_polynomial(io, style, q)
    group && write(io, close)
    return
end

function show_qexpr(io::IO, style::S, q::QExpr, parent_precedence::Int) where {S}
    group = qexpr_precedence(q) < parent_precedence
    open, close = group_delimiters(style)
    group && write(io, open)
    if q.kind == QEXPR_ADD
        show_qexpr_sum(io, style, q)
    elseif q.kind == QEXPR_MUL
        show_qexpr_product(io, style, q)
    else
        show_qexpr_prefactor(io, style, q.coeff)
        fopen, fclose = function_delimiters(style, q.kind)
        write(io, fopen)
        show_qexpr(io, style, only(q.args), 0)
        write(io, fclose)
    end
    group && write(io, close)
    return
end

function show_qexpr_sum(io::IO, style::S, q::QExpr) where {S}
    isempty(q.args) && return write(io, "0")
    show_qexpr(io, style, first(q.args), QEXPR_PREC_ADD)
    for arg in Iterators.drop(q.args, 1)
        if arg isa QExpr && is_real_negative(arg.coeff)
            write(io, " - ")
            show_qexpr(io, style, qexpr_scale(arg, CNUM_NEG1), QEXPR_PREC_ADD)
        else
            write(io, " + ")
            show_qexpr(io, style, arg, QEXPR_PREC_ADD)
        end
    end
    return
end

function show_qexpr_product(io::IO, style::S, q::QExpr) where {S}
    isempty(q.args) && return show_scalar(io, style, q.coeff)
    show_qexpr_prefactor(io, style, q.coeff)
    show_qexpr(io, style, first(q.args), QEXPR_PREC_MUL)
    for arg in Iterators.drop(q.args, 1)
        write(io, product_separator(style))
        show_qexpr(io, style, arg, QEXPR_PREC_MUL)
    end
    return
end

show_qexpr(io::IO, q::QExpr, parent_precedence::Int) =
    show_qexpr(io, TextStyle(), q, parent_precedence)
qexpr_latex(q::QExpr, parent_precedence::Int) =
    sprint(show_qexpr, LatexStyle(), q, parent_precedence)

Base.show(io::IO, q::QExpr) = show_qexpr(io, q, 0)

@latexrecipe function f(q::QExpr)
    return Expr(:latexifymerge, qexpr_latex(q, 0))
end
