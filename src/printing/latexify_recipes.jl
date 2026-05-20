const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool) -> Bool

Set whether [`Transition`](@ref) level indices are rendered as superscripts (`true`, default)
or subscripts (`false`) in LaTeX output via Latexify.jl.

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

function _latex_index_suffix(idx::Index)
    has_index(idx) || return ""
    return "_{$(idx.name)}"
end

# Render an operator name for LaTeX. A bare name (`:a`) passes through, but a
# Julia-style compound name (`:a_pol`, `:c_bog`) is split at the first `_` and
# rendered as `a_{\mathrm{pol}}` — otherwise KaTeX reads the `_` as a subscript
# operator and the trailing characters render as stray italic letters.
function _latex_name(name)
    s = string(name)
    idx = findfirst('_', s)
    idx === nothing && return s
    head = s[1:prevind(s, idx)]
    rest = s[nextind(s, idx):end]
    rest_escaped = replace(rest, '_' => "\\_")
    return string(head, "_{\\mathrm{", rest_escaped, "}}")
end

@latexrecipe function f(x::Destroy)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "$(_latex_name(x.name))$(suffix)")
end

@latexrecipe function f(x::Create)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "$(_latex_name(x.name))$(suffix)^{\\dagger}")
end

@latexrecipe function f(x::Transition)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(_latex_name(x.name))}$(suffix)$(transition_idx_script[]){{$(x.i)$(x.j)}}")
end

@latexrecipe function f(x::Pauli)
    ax = _xyz_sym[x.axis]
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(_latex_name(x.name))}$(suffix)_{{$ax}}")
end

@latexrecipe function f(x::Spin)
    ax = _xyz_sym[x.axis]
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(_latex_name(x.name))}$(suffix)_{{$ax}}")
end

@latexrecipe function f(x::Position)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "\\hat{$(_latex_name(x.name))}$(suffix)")
end

@latexrecipe function f(x::Momentum)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "\\hat{$(_latex_name(x.name))}$(suffix)")
end

# Extract a plain Julia number from CNum for LaTeX rendering
function _latex_prefactor(c::CNum)
    r_unwrap = SymbolicUtils.unwrap(real(c))
    i_unwrap = SymbolicUtils.unwrap(imag(c))
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
        # symbolic prefactors we fall through to the full `c` form below.
        if i_val isa Real
            return complex(false, i_val)
        end
        return c
    elseif r_val isa Real && i_val isa Real
        return complex(r_val, i_val)
    else
        return c
    end
end
_latex_prefactor(c::Number) = c

const _LATEX_TERM = Union{Expr, Number, SymbolicUtils.BasicSymbolic}
const _LATEX_FRAGMENT = Union{String, Symbol, QSym, Number, SymbolicUtils.BasicSymbolic}

# Check if a symbolic prefactor needs \left( \right) brackets when followed by operators.
# Fractions (/) and sums (+) are visually ambiguous without grouping.
function _needs_pf_brackets(pf::Number)
    return false
end
function _needs_pf_brackets(pf::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(pf) || return false
    op = SymbolicUtils.operation(pf)
    return op === (/) || op === (+)
end

# Helper: render a single term (prefactor * operators) as LaTeX
function _latex_term(c::CNum, ops::Vector{QSym})
    pf = _latex_prefactor(c)
    if isempty(ops)
        return pf
    end
    parts = _LATEX_FRAGMENT[]
    if pf isa Number && pf == -1
        push!(parts, :(-))
    elseif pf isa Number && isone(pf)
        # skip prefactor
    elseif _needs_pf_brackets(pf)
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

function _latex_sum_prefix(indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    isempty(indices) && return ""
    idx_parts = String[]
    last_name = indices[end].name
    for (k, idx) in enumerate(indices)
        r = Symbolics.value(SymbolicUtils.unwrap(idx.range))
        if k == length(indices) && !isempty(ne_pairs)
            chain = string(idx.name)
            for (a, b) in ne_pairs
                if a.name == last_name
                    chain *= "{\\neq}$(b.name)"
                else
                    chain *= ",$(a.name){\\neq}$(b.name)"
                end
            end
            push!(idx_parts, "\\underset{$chain}{\\overset{$r}{\\sum}}")
        else
            push!(idx_parts, "\\underset{$(idx.name)}{\\overset{$r}{\\sum}}")
        end
    end
    return join(idx_parts, " ")
end

function _latex_sum_group(indices::Vector{Index}, ne_pairs::Vector{NonEqualPair}, terms::Vector{QAdd})
    prefix = _latex_sum_prefix(indices, ne_pairs)
    term_exprs = _LATEX_TERM[
        let
                ops, c, _ = _term_signature(t)
                _latex_term(c, ops)
        end for t in terms
    ]
    if length(term_exprs) == 1
        return Expr(:latexifymerge, prefix, term_exprs[1])
    end
    return Expr(:latexifymerge, prefix, Expr(:call, :+, term_exprs...))
end

@latexrecipe function f(x::QAdd)
    st = sorted_arguments(x)
    if !isempty(x.indices)
        # Split terms into index-dependent and index-independent
        dep_qadds = QAdd[]
        terms_out = _LATEX_TERM[]
        for t in st
            ops, c, _ = _term_signature(t)
            term_expr = _latex_term(c, ops)
            if any(idx -> _depends_on_index_term(c, ops, idx), x.indices)
                push!(dep_qadds, t)
            else
                push!(terms_out, term_expr)
            end
        end
        if !isempty(dep_qadds)
            for (used, ne_pairs, terms) in _group_dep_terms(dep_qadds, x.indices)
                push!(terms_out, _latex_sum_group(used, ne_pairs, terms))
            end
        end
        if isempty(terms_out)
            return 0
        elseif length(terms_out) == 1
            return only(terms_out)
        end
        return Expr(:call, :+, terms_out...)
    end
    terms = _LATEX_TERM[
        let
                ops, c, _ = _term_signature(t)
                _latex_term(c, ops)
        end for t in st
    ]
    return Expr(:call, :+, terms...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
