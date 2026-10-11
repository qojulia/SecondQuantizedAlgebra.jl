const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool) -> Bool

Set whether [`Transition`](@ref) and [`CollectiveTransition`](@ref) level indices
are rendered as superscripts (`true`, default) or subscripts (`false`) in LaTeX
output via Latexify.jl.

- `true`: ``{\\sigma}^{{ij}}``
- `false`: ``{\\sigma}_{{ij}}``

# Examples
```jldoctest
julia> SecondQuantizedAlgebra.transition_superscript(false)
false

julia> SecondQuantizedAlgebra.transition_superscript(true)
true
```
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

function needs_pf_brackets(pf::Number)
    return false
end
function needs_pf_brackets(pf::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(pf) || return false
    op = SymbolicUtils.operation(pf)
    return op === (/) || op === (+)
end

function latex_term(c::CNum, ops::Vector{Op})
    is_native(c) && return latex_native_term(native_scalar(c), ops)
    re, im = realimag(c)
    iszero_num(im) && return latex_pf_term(SymbolicUtils.unwrap(re), ops)
    return latex_pf_term(Complex(re, im), ops)
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
    idx_parts = String[]
    last_name = index_name(indices[end])
    for (k, idx) in enumerate(indices)
        r = Symbolics.value(SymbolicUtils.unwrap(index_range(idx)))
        if k == length(indices) && !isempty(ne_pairs)
            chain = string(index_name(idx))
            for (a, b) in ne_pairs
                if index_name(a) == last_name
                    chain *= "{\\neq}$(index_name(b))"
                else
                    chain *= ",$(index_name(a)){\\neq}$(index_name(b))"
                end
            end
            push!(idx_parts, "\\underset{$chain}{\\overset{$r}{\\sum}}")
        else
            push!(idx_parts, "\\underset{$(index_name(idx))}{\\overset{$r}{\\sum}}")
        end
    end
    return join(idx_parts, " ")
end

function latex_sum_group(indices::Vector{Index}, ne_pairs::Vector{NonEqualPair}, terms::Vector{QAdd})
    prefix = latex_sum_prefix(indices, ne_pairs)
    term_exprs = Expr[
        let
            ops, c, _ = term_signature(t)
            latex_term(c, ops)
        end for t in terms
    ]
    if length(term_exprs) == 1
        return Expr(:latexifymerge, prefix, term_exprs[1])
    end
    return Expr(:latexifymerge, prefix, Expr(:call, :+, term_exprs...))
end

function latex_fragment(x)
    inner = SymbolicUtils.isconst(x) ? x.val : x
    return strip(String(latexify(inner)), '$')
end

function latex_avg_expr(inner)
    return "\\langle $(latex_fragment(inner)) \\rangle"
end

function Symbolics._toexpr_op(::AvgFunc, args; kwargs...)
    return latex_avg_expr(only(args))
end

function Symbolics._toexpr_op(::typeof(expim), args; kwargs...)
    arg = only(args)
    neg = leading_sign(arg) < 0
    body = neg ? SymbolicUtils.unwrap(expand(-arg)) : arg
    return "e^{$(neg ? "-" : "")i $(strip(String(latexify(body; env = :inline)), '$'))}"
end

function Symbolics._toexpr_op(::SumFunc, args; kwargs...)
    scope = scope_of(args[2])
    body = args[1]
    body_tex = strip(String(latexify(body; env = :inline)), '$')
    if SymbolicUtils.iscall(body) && SymbolicUtils.operation(body) === (+)
        body_tex = string("\\left( ", body_tex, " \\right)")
    end
    return string(latex_sum_prefix(scope.indices, scope.ne), body_tex)
end

function Symbolics._toexpr_metadata(
        x::SymbolicUtils.BasicSymbolic, ::Type{AverageOperator}, op;
        kwargs...,
    )
    iv_tex = strip(String(latexify(only(SymbolicUtils.arguments(x)); env = :inline)), '$')
    return string(latex_avg_expr(op), "\\left( ", iv_tex, " \\right)")
end

@latexrecipe function f(x::QAdd)
    st = sorted_arguments(x)
    if !isempty(x.indices)
        dep_qadds = QAdd[]
        terms_out = Expr[]
        for t in st
            ops, c, _ = term_signature(t)
            term_expr = latex_term(c, ops)
            if any(idx -> depends_on_index_term(c, ops, idx), x.indices)
                push!(dep_qadds, t)
            else
                push!(terms_out, term_expr)
            end
        end
        if !isempty(dep_qadds)
            for (used, ne_pairs, terms) in group_dep_terms(dep_qadds, x.indices)
                push!(terms_out, latex_sum_group(used, ne_pairs, terms))
            end
        end
        if isempty(terms_out)
            return 0
        elseif length(terms_out) == 1
            return only(terms_out)
        end
        return Expr(:call, :+, terms_out...)
    end
    terms = Expr[
        let
            ops, c, _ = term_signature(t)
            latex_term(c, ops)
        end for t in st
    ]
    return Expr(:call, :+, terms...)
end

const QLaTeX = Union{QField, Coeff}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{<:QLaTeX}) =
    print(io, "\$\$ ", latexify(x; env = :raw), " \$\$")
