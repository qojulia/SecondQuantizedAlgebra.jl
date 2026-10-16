const SuperMap = Union{SuperOperator, SuperExpr, Lindbladian}

function sorted_by_key(key::F, d::AbstractDict) where {F}
    ps = collect(d)
    return ps[sortperm(map(key, ps))]
end

sorted_superterms(S::SuperOperator) = sorted_by_key(p -> superterm_order_key(p.first), S.arguments)
sorted_formal(X::SuperExpr) =
    sorted_by_key(p -> (sprint(show, p.first.left), sprint(show, p.first.right)), X.formal)

function displays_negative(c::CNum)
    is_real_negative(c) && return true
    is_native(c) || return false
    z = native_scalar(c)
    return real(z) < 0 || (iszero(real(z)) && imag(z) < 0)
end

function write_sign!(io::IO, first_term::Bool, c::CNum)
    negative = !first_term && displays_negative(c)
    first_term || write(io, negative ? " - " : " + ")
    return negative ? -c : c
end

function latex_coeff_prefix(c::CNum)
    is_unit(c) && return ""
    is_neg_unit(c) && return "-"
    displays_negative(c) && return string("-", latex_coeff_prefix(-c))
    body = latex_inline(c)
    if is_native(c)
        needs = has_both_parts(native_scalar(c))
    else
        d = to_num(c)
        needs = needs_pf_parens(d) || has_both_parts(d)
    end
    return needs ? string("\\left( ", body, " \\right) ") : string(body, " ")
end

sandwich_tokens(::TextStyle) = (" ⋅ ", "(·)")
sandwich_tokens(::LatexStyle) = ("\\,", "\\cdot")

commutator_delimiters(::TextStyle) = ("-i[", ", ·]")
commutator_delimiters(::LatexStyle) = ("-i\\left[", ",\\cdot\\right]")

dissipator_delimiters(::TextStyle) = ("𝒟[", ", ", "]")
dissipator_delimiters(::LatexStyle) = ("\\mathcal{D}\\left[", ",", "\\right]")

render_ops(::TextStyle, ops::Vector{Op}) = sprint(show_ops, ops)
render_ops(style::LatexStyle, ops::Vector{Op}) = join((latex_inline(op) for op in ops), product_separator(style))

render_side(style::PrintStyle, q::QField) = isone(q) ? "" : sprint(show_qexpr, style, q, QEXPR_PREC_MUL)

show_rate(io::IO, ::TextStyle, c::CNum, ::Bool) = show_scaled_prefactor(io, c)
function show_rate(io::IO, ::LatexStyle, c::CNum, dissipator::Bool)
    prefix = latex_coeff_prefix(c)
    dissipator && !isempty(prefix) && return write(io, rstrip(prefix), "\\,")
    return write(io, prefix)
end

function show_sandwich(io::IO, style::PrintStyle, left::String, right::String)
    sep, dot = sandwich_tokens(style)
    isempty(left) || write(io, left, sep)
    write(io, dot)
    isempty(right) || write(io, sep, right)
    return
end

show_superterm(io::IO, style::PrintStyle, t::SuperTerm) =
    show_sandwich(io, style, render_ops(style, t.left), render_ops(style, t.right))

function show_channel(io::IO, style::PrintStyle, ch::Channel, rate::CNum)
    isempty(ch.scope) || show_sum_prefix(io, style, ch.scope, ch.ne)
    show_rate(io, style, rate, true)
    open, sep, close = dissipator_delimiters(style)
    write(io, open)
    show_qexpr(io, style, ch.jump, 0)
    if !is_diagonal(ch)
        write(io, sep)
        show_qexpr(io, style, ch.partner, 0)
    end
    write(io, close)
    return
end

function show_map(io::IO, style::PrintStyle, S::SuperOperator)
    isempty(S.arguments) && return write(io, "0")
    first_term = true
    for (t, c) in sorted_superterms(S)
        c = write_sign!(io, first_term, c)
        first_term = false
        used = bound_indices(t, c, S.indices)
        isempty(used) || show_sum_prefix(io, style, used, t.ne)
        show_rate(io, style, c, false)
        show_superterm(io, style, t)
    end
    return
end

function show_map(io::IO, style::PrintStyle, L::Lindbladian)
    first_term = iszero(L.hamiltonian)
    if !first_term
        open, close = commutator_delimiters(style)
        write(io, open)
        show_qexpr(io, style, L.hamiltonian, 0)
        write(io, close)
    end
    for (ch, γ) in channels(L)
        γ = write_sign!(io, first_term, γ)
        first_term = false
        show_channel(io, style, ch, γ)
    end
    first_term && write(io, "0")
    return
end

function show_map(io::IO, style::PrintStyle, X::SuperExpr)
    first_term = iszero(X.poly)
    first_term || show_map(io, style, X.poly)
    for (t, c) in sorted_formal(X)
        c = write_sign!(io, first_term, c)
        first_term = false
        show_rate(io, style, c, false)
        show_sandwich(io, style, render_side(style, t.left), render_side(style, t.right))
    end
    first_term && write(io, "0")
    return
end

Base.show(io::IO, x::SuperMap) = show_map(io, TextStyle(), x)
Base.show(io::IO, t::SuperTerm) = show_superterm(io, TextStyle(), t)
Base.show(io::IO, ch::Channel) = show_channel(io, TextStyle(), ch, CNUM_ONE)

@latexrecipe function f(x::SuperMap)
    return Expr(:latexifymerge, sprint(show_map, LatexStyle(), x))
end
