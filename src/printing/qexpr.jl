# Precedence-aware rendering for the cold-path formal expression layer.

const QEXPR_PREC_ADD = 1
const QEXPR_PREC_MUL = 2
const QEXPR_PREC_FUNC = 3

function qexpr_precedence(q::QExpr)
    q.kind == QEXPR_ADD && return QEXPR_PREC_ADD
    (q.kind == QEXPR_MUL || !is_unit(q.coeff)) && return QEXPR_PREC_MUL
    return QEXPR_PREC_FUNC
end

qexpr_qadd_needs_group(q::QAdd) = !isempty(q.indices) || length(q.arguments) > 1

function show_qexpr_prefactor(io::IO, c::CNum)
    if is_neg_unit(c)
        write(io, "-")
        return
    elseif is_unit(c)
        return
    end

    tail = c.tail
    if tail isa RawSymbolicCoeff
        expr = tail.expr
        brace = SymbolicUtils.iscall(expr) &&
            (SymbolicUtils.operation(expr) === (+) || SymbolicUtils.operation(expr) === (/))
        brace && write(io, "(")
        print(io, expr)
        brace && write(io, ")")
    else
        d = to_num(c)
        brace = !is_native(c) && needs_pf_parens(d)
        brace && write(io, "(")
        show_display(io, d)
        brace && write(io, ")")
    end
    write(io, " * ")
    return
end

function show_qexpr_arg(io::IO, arg::QAdd, parent_precedence::Int)
    group = parent_precedence > QEXPR_PREC_ADD && qexpr_qadd_needs_group(arg)
    group && write(io, "(")
    show(io, arg)
    group && write(io, ")")
    return
end
show_qexpr_arg(io::IO, arg::QExpr, parent_precedence::Int) =
    show_qexpr(io, arg, parent_precedence)

function show_qexpr_add(io::IO, q::QExpr)
    isempty(q.args) && return write(io, "0")
    show_qexpr_arg(io, first(q.args), QEXPR_PREC_ADD)
    for arg in Iterators.drop(q.args, 1)
        if arg isa QExpr && is_real_negative(arg.coeff)
            write(io, " - ")
            show_qexpr(io, qexpr_scale(arg, CNUM_NEG1), QEXPR_PREC_ADD)
        else
            write(io, " + ")
            show_qexpr_arg(io, arg, QEXPR_PREC_ADD)
        end
    end
    return
end

function show_qexpr_mul(io::IO, q::QExpr)
    isempty(q.args) && return show_prefactor(io, q.coeff)
    show_qexpr_prefactor(io, q.coeff)
    show_qexpr_arg(io, first(q.args), QEXPR_PREC_MUL)
    for arg in Iterators.drop(q.args, 1)
        write(io, " * ")
        show_qexpr_arg(io, arg, QEXPR_PREC_MUL)
    end
    return
end

function show_qexpr_function(io::IO, q::QExpr)
    show_qexpr_prefactor(io, q.coeff)
    arg = only(q.args)
    if q.kind == QEXPR_SIN
        write(io, "sin(")
        show_qexpr_arg(io, arg, 0)
        write(io, ")")
    elseif q.kind == QEXPR_COS
        write(io, "cos(")
        show_qexpr_arg(io, arg, 0)
        write(io, ")")
    else
        write(io, "exp(im*(")
        show_qexpr_arg(io, arg, 0)
        write(io, "))")
    end
    return
end

function show_qexpr(io::IO, q::QExpr, parent_precedence::Int)
    precedence = qexpr_precedence(q)
    group = precedence < parent_precedence
    group && write(io, "(")
    if q.kind == QEXPR_ADD
        show_qexpr_add(io, q)
    elseif q.kind == QEXPR_MUL
        show_qexpr_mul(io, q)
    else
        show_qexpr_function(io, q)
    end
    group && write(io, ")")
    return
end

Base.show(io::IO, q::QExpr) = show_qexpr(io, q, 0)

latex_inline(x) = strip(String(latexify(x; env = :inline)), '$')

function qexpr_latex_prefactor(c::CNum)
    is_neg_unit(c) && return "-"
    is_unit(c) && return ""
    return string(latex_inline(latex_prefactor(c)), " ")
end

function qexpr_latex_arg(arg::QAdd, parent_precedence::Int)
    body = latex_inline(arg)
    if parent_precedence > QEXPR_PREC_ADD && qexpr_qadd_needs_group(arg)
        return string("\\left( ", body, " \\right)")
    end
    return body
end
qexpr_latex_arg(arg::QExpr, parent_precedence::Int) =
    qexpr_latex(arg, parent_precedence)

function qexpr_latex_add(q::QExpr)
    isempty(q.args) && return "0"
    io = IOBuffer()
    write(io, qexpr_latex_arg(first(q.args), QEXPR_PREC_ADD))
    for arg in Iterators.drop(q.args, 1)
        if arg isa QExpr && is_real_negative(arg.coeff)
            write(io, " - ")
            write(io, qexpr_latex(qexpr_scale(arg, CNUM_NEG1), QEXPR_PREC_ADD))
        else
            write(io, " + ")
            write(io, qexpr_latex_arg(arg, QEXPR_PREC_ADD))
        end
    end
    return String(take!(io))
end

function qexpr_latex_mul(q::QExpr)
    isempty(q.args) && return latex_inline(latex_prefactor(q.coeff))
    io = IOBuffer()
    write(io, qexpr_latex_prefactor(q.coeff))
    write(io, qexpr_latex_arg(first(q.args), QEXPR_PREC_MUL))
    for arg in Iterators.drop(q.args, 1)
        write(io, " ")
        write(io, qexpr_latex_arg(arg, QEXPR_PREC_MUL))
    end
    return String(take!(io))
end

function qexpr_latex_function(q::QExpr)
    prefix = qexpr_latex_prefactor(q.coeff)
    body = qexpr_latex_arg(only(q.args), 0)
    if q.kind == QEXPR_SIN
        return string(prefix, "\\sin\\left( ", body, " \\right)")
    elseif q.kind == QEXPR_COS
        return string(prefix, "\\cos\\left( ", body, " \\right)")
    end
    return string(prefix, "e^{i\\left( ", body, " \\right)}")
end

function qexpr_latex(q::QExpr, parent_precedence::Int)
    precedence = qexpr_precedence(q)
    body = if q.kind == QEXPR_ADD
        qexpr_latex_add(q)
    elseif q.kind == QEXPR_MUL
        qexpr_latex_mul(q)
    else
        qexpr_latex_function(q)
    end
    if precedence < parent_precedence
        return string("\\left( ", body, " \\right)")
    end
    return body
end

@latexrecipe function f(q::QExpr)
    return Expr(:latexifymerge, qexpr_latex(q, 0))
end
