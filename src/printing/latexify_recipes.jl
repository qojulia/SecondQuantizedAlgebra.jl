const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool)

Choose whether transition levels are written as superscripts or subscripts in LaTeX.

With `true`, the default, ``\\sigma^{ij}`` is written; with `false`,
``\\sigma_{ij}``. The setting is global and applies to every later call of
`latexify`. The function is not exported.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: transition_superscript

julia> using Latexify

julia> σ = Transition(NLevelSpace(:atom, 2), :σ, 2, 1);

julia> transition_superscript(false); latexify(σ)
L"\${\\sigma}_{{21}}\$"

julia> transition_superscript(true); latexify(σ)
L"\${\\sigma}^{{21}}\$"
```

See also [`Transition`](@ref).
"""
function transition_superscript(x::Bool)
    transition_idx_script[] = x ? :^ : :_
    return x
end

function latex_index_suffix(idx::Index)
    has_index(idx) || return ""
    parts = split(string(index_name(idx)), '_')
    name = length(parts) == 1 ? parts[1] : string(parts[1], "_{", join(parts[2:end], ","), "}")
    return "_{$(name)}"
end

function latex_name(name)
    s = string(name)
    occursin('_', s) || return s
    head, rest = split(s, '_'; limit = 2)
    rest_escaped = replace(rest, '_' => "\\_")
    return string(head, "_{\\mathrm{", rest_escaped, "}}")
end

@latexrecipe function f(x::Op)
    suffix = latex_index_suffix(x.index)
    name = latex_name(operator_name(x))
    k = x.kind
    body = if k === OP_DESTROY
        "$(name)$(suffix)"
    elseif k === OP_CREATE
        "$(name)$(suffix)^{\\dagger}"
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        "{$(name)}$(suffix)$(transition_idx_script[]){{$(Int(x.l1))$(Int(x.l2))}}"
    elseif k === OP_PAULI || k === OP_SPIN
        "{$(name)}$(suffix)_{{$(xyz_sym[x.l1])}}"
    else
        "\\hat{$(name)}$(suffix)"
    end
    return Expr(:latexifymerge, body)
end

latex_real(x::Int) = x
latex_real(x::Float64) = x
latex_real(x::Rational{Int}) = isone(denominator(x)) ? numerator(x) : x

function latex_native_term(z::Union{ComplexF64, GaussianRational{Int}}, ops::Vector{Op})
    r = latex_real(real(z))
    i = latex_real(imag(z))
    iszero(imag(z)) && return latex_pf_term(r, ops)
    iszero(real(z)) && return latex_pf_term(complex(false, i), ops)
    return latex_pf_term(complex(r, i), ops)
end

@latexrecipe function f(c::Coeff)
    return latex_term(c, Op[])
end

const LATEX_FRAGMENT = Union{String, Symbol, QSym, Number, SymbolicUtils.BasicSymbolic}

needs_pf_brackets(pf::Number) = false
needs_pf_brackets(pf::Complex) = has_both_parts(pf)
has_both_parts(z) = !is_zero_part(real(z)) && !is_zero_part(imag(z))
is_zero_part(x::Num) = iszero_num(x)
is_zero_part(x::Real) = iszero(x)
needs_pf_brackets(pf::SymbolicUtils.BasicSymbolic) = is_call_to(pf, (/, +))

function latex_term(c::CNum, ops::Vector{Op})
    is_native(c) && return latex_native_term(native_scalar(c), ops)
    tail = c.tail
    tail isa Poly && has_complex_atom(tail) && return latex_pf_term(poly_to_raw(tail), ops)
    re, im = realimag(c)
    iszero_num(im) && return latex_pf_term(SymbolicUtils.unwrap(re), ops)
    return latex_pf_term(Complex(re, im), ops)
end

function latex_term(t::QAdd)
    ops, c, _ = term_signature(t)
    return latex_term(c, ops)
end

function latex_pf_term(pf, ops::Vector{Op})
    isempty(ops) && return Expr(:latexifymerge, pf)
    parts = LATEX_FRAGMENT[]
    if pf isa Number && pf == -1
        push!(parts, :(-))
    elseif pf isa Number && isone(pf)
    elseif needs_pf_brackets(pf)
        push!(parts, "\\left(")
        push!(parts, pf)
        push!(parts, "\\right) ")
    else
        push!(parts, pf)
        push!(parts, " ")
    end
    for op in ops
        push!(parts, op)
    end
    return Expr(:latexifymerge, parts...)
end

function latex_sum_prefix(indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    isempty(indices) && return ""
    last_name = index_name(indices[end])
    constraint = join(
        index_name(a) == last_name ? "{\\neq}$(index_name(b))" : ",$(index_name(a)){\\neq}$(index_name(b))"
            for (a, b) in ne_pairs
    )
    parts = map(enumerate(indices)) do (k, idx)
        r = Symbolics.value(SymbolicUtils.unwrap(index_range(idx)))
        label = string(index_name(idx), k == length(indices) ? constraint : "")
        "\\underset{$label}{\\overset{$r}{\\sum}}"
    end
    return join(parts, " ")
end

show_sum_prefix(io::IO, ::LatexStyle, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair}) =
    write(io, latex_sum_prefix(indices, ne_pairs), " ")

function latex_sum_group(indices::Vector{Index}, ne_pairs::Vector{NonEqualPair}, terms::Vector{QAdd})
    prefix = latex_sum_prefix(indices, ne_pairs)
    length(terms) == 1 && return Expr(:latexifymerge, prefix, latex_term(only(terms)))
    return Expr(:latexifymerge, prefix, Expr(:call, :+, map(latex_term, terms)...))
end

function latex_string(@nospecialize(tex))
    tex isa LaTeXString || throw(ArgumentError("Latexify returned a $(typeof(tex)), not a LaTeXString"))
    return tex
end

latex_body(@nospecialize(tex)) = String(strip(latex_string(tex), '$'))

latex_inline(x) = latex_body(latexify(x; env = :inline))

function latex_avg_expr(x)
    inner = SymbolicUtils.isconst(x) ? x.val : x
    return "\\langle $(latex_body(latexify(inner))) \\rangle"
end

function Symbolics._toexpr_op(::AvgFunc, args; kwargs...)
    return latex_avg_expr(only(args))
end

function Symbolics._toexpr_op(::typeof(expim), args; kwargs...)
    neg, body = split_phase_sign(only(args))
    return "e^{$(neg ? "-" : "")i $(latex_inline(body))}"
end

function Symbolics._toexpr_op(::SumFunc, args; kwargs...)
    scope = scope_of(args[2])
    body = args[1]
    body_tex = latex_inline(body)
    if is_call_to(body, (+,))
        body_tex = string("\\left( ", body_tex, " \\right)")
    end
    return string(latex_sum_prefix(scope.indices, scope.ne), body_tex)
end

function Symbolics._toexpr_metadata(
        x::SymbolicUtils.BasicSymbolic, ::Type{AverageOperator}, op;
        kwargs...,
    )
    iv_tex = latex_inline(only(SymbolicUtils.arguments(x)))
    return string(latex_avg_expr(op), "\\left( ", iv_tex, " \\right)")
end

@latexrecipe function f(x::QAdd)
    st = sorted_arguments(x)
    isempty(st) && return 0
    isempty(x.indices) && return Expr(:call, :+, map(latex_term, st)...)
    dep, indep = split_index_dependent(st, x.indices)
    terms = Expr[latex_term(t) for t in indep]
    for (used, ne_pairs, group) in group_dep_terms(dep, x.indices)
        push!(terms, latex_sum_group(used, ne_pairs, group))
    end
    return length(terms) == 1 ? only(terms) : Expr(:call, :+, terms...)
end

const QLaTeX = Union{QField, Coeff, SuperOperator, SuperExpr, Lindbladian}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latex_string(latexify(x)))
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{<:QLaTeX}) =
    print(io, "\$\$ ", latex_string(latexify(x; env = :raw)), " \$\$")
