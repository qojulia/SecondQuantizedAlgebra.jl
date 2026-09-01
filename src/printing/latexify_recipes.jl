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

# Join accumulated per-slot suffixes into one comma subscript (`i_2_1` -> `i_{2,1}`);
# a bare `_{i_2_1}` is a double subscript that MathJax rejects.
function latex_index_suffix(idx::Index)
    has_index(idx) || return ""
    parts = split(string(index_name(idx)), '_')
    name = length(parts) == 1 ? parts[1] : string(parts[1], "_{", join(parts[2:end], ","), "}")
    return "_{$(name)}"
end

# Render an operator name for LaTeX. A bare name (`:a`) passes through, but a
# Julia-style compound name (`:a_pol`, `:c_bog`) is split at the first `_` and
# rendered as `a_{\mathrm{pol}}` — otherwise KaTeX reads the `_` as a subscript
# operator and the trailing characters render as stray italic letters.
function latex_name(name)
    s = string(name)
    idx = findfirst('_', s)
    idx === nothing && return s
    head = s[1:prevind(s, idx)]
    rest = s[nextind(s, idx):end]
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
    else   # OP_POSITION, OP_MOMENTUM
        "\\hat{$(name)}$(suffix)"
    end
    return Expr(:latexifymerge, body)
end

# Lower through the same public coefficient representation used by terminal display.
function latex_prefactor(c::CNum)
    d = to_num(c)
    r_unwrap = SymbolicUtils.unwrap(real(d))
    i_unwrap = SymbolicUtils.unwrap(imag(d))
    r_val = Symbolics.value(r_unwrap)
    i_val = Symbolics.value(i_unwrap)
    # `iszero` on a BasicSymbolic returns a symbolic expression, not Bool, which
    # blows up the `if` below on older Symbolics (Julia 1.10 CI). Use structural
    # equality on the unwrapped form to get a Bool either way.
    i_is_zero = isequal(i_unwrap, 0) || (i_val isa Number && iszero(i_val))
    r_is_zero = isequal(r_unwrap, 0) || (r_val isa Number && iszero(r_val))
    if i_is_zero
        return r_val
    elseif r_is_zero
        # Pure imaginary: `complex(false, x)` only works for `x <: Real`, so on
        # symbolic prefactors we fall through to the full `Complex{Num}` form below.
        if i_val isa Real
            return complex(false, i_val)
        end
        return d
    elseif r_val isa Real && i_val isa Real
        return complex(r_val, i_val)
    else
        return d
    end
end
latex_prefactor(c::Number) = c

const LATEX_TERM = Union{Expr, Number, SymbolicUtils.BasicSymbolic}
const LATEX_FRAGMENT = Union{String, Symbol, QSym, Number, SymbolicUtils.BasicSymbolic}

# Check if a symbolic prefactor needs \left( \right) brackets when followed by operators.
# Fractions (/) and sums (+) are visually ambiguous without grouping.
function needs_pf_brackets(pf::Number)
    return false
end
function needs_pf_brackets(pf::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(pf) || return false
    op = SymbolicUtils.operation(pf)
    return op === (/) || op === (+)
end

# Helper: render a single term (prefactor * operators) as LaTeX
function latex_term(c::CNum, ops::Vector{Op})
    pf = latex_prefactor(c)
    if isempty(ops)
        return pf
    end
    parts = LATEX_FRAGMENT[]
    if pf isa Number && pf == -1
        push!(parts, :(-))
    elseif pf isa Number && isone(pf)
        # skip prefactor
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
    term_exprs = LATEX_TERM[
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

# The LaTeX twin of the `show_call(::typeof(expim))` override, so a phase renders
# as an exponential rather than as the raw head name.
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
    # Wrap a multi-term (`+`) body so the Σ binds the whole sum, not just its
    # first summand (mirrors the Unicode `show_call(::SumFunc)`).
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
        # Split terms into index-dependent and index-independent
        dep_qadds = QAdd[]
        terms_out = LATEX_TERM[]
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
    terms = LATEX_TERM[
        let
            ops, c, _ = term_signature(t)
            latex_term(c, ops)
        end for t in st
    ]
    return Expr(:call, :+, terms...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
