function show_coeff_prefix(io::IO, c::CNum)
    is_unit(c) && return
    if is_neg_unit(c)
        write(io, "-")
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
        brace = !is_native(c) && (needs_pf_parens(d) || !iszero_num(imag(d)) && !iszero_num(real(d)))
        brace && write(io, "(")
        show_display(io, d)
        brace && write(io, ")")
    end
    write(io, " * ")
    return
end

function show_ops(io::IO, ops::Vector{Op})
    for (k, op) in enumerate(ops)
        k > 1 && write(io, " * ")
        show(io, op)
    end
    return
end

function show_sandwich(io::IO, left::Vector{Op}, right::Vector{Op})
    if !isempty(left)
        show_ops(io, left)
        write(io, " ⋅ ")
    end
    write(io, "(·)")
    if !isempty(right)
        write(io, " ⋅ ")
        show_ops(io, right)
    end
    return
end

function sorted_superterms(S::SuperOperator)
    ps = collect(S.arguments)
    sort!(ps; by = p -> superterm_order_key(p.first))
    return ps
end

function displays_negative(c::CNum)
    is_real_negative(c) && return true
    is_native(c) || return false
    z = native_scalar(c)
    return real(z) < 0 || (iszero(real(z)) && imag(z) < 0)
end

function write_sign!(io::IO, first_term::Bool, c::CNum)
    negative = !first_term && displays_negative(c)
    first_term || write(io, negative ? " - " : " + ")
    return negative ? neg_cnum(c) : c
end

function Base.show(io::IO, S::SuperOperator)
    isempty(S.arguments) && return write(io, "0")
    first_term = true
    for (t, c) in sorted_superterms(S)
        c = write_sign!(io, first_term, c)
        first_term = false
        used = bound_indices(t, c, S.indices)
        if !isempty(used)
            show_sum_prefix(io, used, t.ne)
            write(io, " ")
        end
        show_coeff_prefix(io, c)
        show_sandwich(io, t.left, t.right)
    end
    return
end

function show_side(io::IO, q::QField)
    s = sprint(show, q)
    print(io, s)
    return
end

function sorted_channels(L::Lindbladian)
    ps = collect(L.channels)
    sort!(ps; by = p -> channel_order_key(p.first))
    return ps
end

function show_scope_prefix(io::IO, ch::Channel)
    if !isempty(ch.scope)
        show_sum_prefix(io, ch.scope, ch.ne)
        write(io, " ")
    end
    return nothing
end

function show_dissipator(io::IO, ch::Channel)
    write(io, "𝒟[")
    show_side(io, ch.jump)
    if !is_diagonal(ch)
        write(io, ", ")
        show_side(io, ch.partner)
    end
    write(io, "]")
    return nothing
end

function sorted_formal(X::SuperExpr)
    ps = collect(X.formal)
    sort!(ps; by = p -> (sprint(show, p.first.left), sprint(show, p.first.right)))
    return ps
end

function Base.show(io::IO, L::Lindbladian)
    first_term = true
    if !iszero(L.hamiltonian)
        write(io, "-i[")
        show_side(io, L.hamiltonian)
        write(io, ", ·]")
        first_term = false
    end
    for (ch, γ) in sorted_channels(L)
        γ = write_sign!(io, first_term, γ)
        first_term = false
        show_scope_prefix(io, ch)
        show_coeff_prefix(io, γ)
        show_dissipator(io, ch)
    end
    first_term && write(io, "0")
    return
end

function Base.show(io::IO, X::SuperExpr)
    first_term = true
    if !iszero(X.poly)
        show(io, X.poly)
        first_term = false
    end
    for (t, c) in sorted_formal(X)
        c = write_sign!(io, first_term, c)
        first_term = false
        show_coeff_prefix(io, c)
        if !isone(t.left)
            show_qexpr(io, t.left, QEXPR_PREC_MUL)
            write(io, " ⋅ ")
        end
        write(io, "(·)")
        if !isone(t.right)
            write(io, " ⋅ ")
            show_qexpr(io, t.right, QEXPR_PREC_MUL)
        end
    end
    first_term && write(io, "0")
    return
end

Base.show(io::IO, t::SuperTerm) = show_sandwich(io, t.left, t.right)
function Base.show(io::IO, ch::Channel)
    show_scope_prefix(io, ch)
    show_dissipator(io, ch)
    return nothing
end

latex_ops(ops::Vector{Op}) = join((latex_inline(op) for op in ops), " ")

function latex_coeff_prefix(c::CNum)
    is_unit(c) && return ""
    is_neg_unit(c) && return "-"
    displays_negative(c) && return string("-", latex_coeff_prefix(neg_cnum(c)))
    body = latex_inline(c)
    if is_native(c)
        z = native_scalar(c)
        needs = !iszero(real(z)) && !iszero(imag(z))
    else
        d = to_num(c)
        needs = needs_pf_parens(d) || (!iszero_num(real(d)) && !iszero_num(imag(d)))
    end
    return needs ? string("\\left( ", body, " \\right) ") : string(body, " ")
end

function latex_sandwich(left::String, right::String)
    isempty(left) && isempty(right) && return "\\cdot"
    isempty(left) && return string("\\cdot\\,", right)
    isempty(right) && return string(left, "\\,\\cdot")
    return string(left, "\\,\\cdot\\,", right)
end

function super_latex(S::SuperOperator)
    isempty(S.arguments) && return "0"
    io = IOBuffer()
    first_term = true
    for (t, c) in sorted_superterms(S)
        c = write_sign!(io, first_term, c)
        first_term = false
        used = bound_indices(t, c, S.indices)
        isempty(used) || write(io, latex_sum_prefix(used, t.ne), " ")
        write(io, latex_coeff_prefix(c))
        write(io, latex_sandwich(latex_ops(t.left), latex_ops(t.right)))
    end
    return String(take!(io))
end

latex_part(q::QAdd) = latex_inline(q)
latex_part(q::QExpr) = qexpr_latex(q, 0)

function super_latex(L::Lindbladian)
    io = IOBuffer()
    first_term = true
    if !iszero(L.hamiltonian)
        write(io, "-i\\left[", latex_part(L.hamiltonian), ",\\cdot\\right]")
        first_term = false
    end
    for (ch, γ) in sorted_channels(L)
        γ = write_sign!(io, first_term, γ)
        first_term = false
        isempty(ch.scope) || write(io, latex_sum_prefix(ch.scope, ch.ne), " ")
        prefix = latex_coeff_prefix(γ)
        write(io, isempty(prefix) ? "" : string(rstrip(prefix), "\\,"))
        write(io, "\\mathcal{D}\\left[", latex_part(ch.jump))
        is_diagonal(ch) || write(io, ",", latex_part(ch.partner))
        write(io, "\\right]")
    end
    first_term && write(io, "0")
    return String(take!(io))
end

function super_latex(X::SuperExpr)
    io = IOBuffer()
    first_term = true
    if !iszero(X.poly)
        write(io, super_latex(X.poly))
        first_term = false
    end
    for (t, c) in sorted_formal(X)
        c = write_sign!(io, first_term, c)
        first_term = false
        write(io, latex_coeff_prefix(c))
        left = isone(t.left) ? "" : qexpr_latex(t.left, QEXPR_PREC_MUL)
        right = isone(t.right) ? "" : qexpr_latex(t.right, QEXPR_PREC_MUL)
        write(io, latex_sandwich(left, right))
    end
    first_term && write(io, "0")
    return String(take!(io))
end

@latexrecipe function f(S::SuperOperator)
    return Expr(:latexifymerge, super_latex(S))
end
@latexrecipe function f(L::Lindbladian)
    return Expr(:latexifymerge, super_latex(L))
end
@latexrecipe function f(X::SuperExpr)
    return Expr(:latexifymerge, super_latex(X))
end

const MapLaTeX = Union{SuperOperator, SuperExpr, Lindbladian}
Base.show(io::IO, ::MIME"text/latex", x::MapLaTeX) = write(io, latex_string(latexify(x)))
