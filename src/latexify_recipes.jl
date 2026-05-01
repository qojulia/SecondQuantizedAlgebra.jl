# Transition superscript toggle
const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool) -> Bool

Set whether [`Transition`](@ref) level indices are rendered as superscripts (`true`, default)
or subscripts (`false`) in LaTeX output via Latexify.jl.

- `true`: ``{\\sigma}^{{ij}}``
- `false`: ``{\\sigma}_{{ij}}``

# Examples
```julia
transition_superscript(false)   # use subscript notation
transition_superscript(true)    # restore default superscript
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

@latexrecipe function f(x::Destroy)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "$(x.name)$(suffix)")
end

@latexrecipe function f(x::Create)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "$(x.name)$(suffix)^{\\dagger}")
end

@latexrecipe function f(x::Transition)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(x.name)}$(suffix)$(transition_idx_script[]){{$(x.i)$(x.j)}}")
end

@latexrecipe function f(x::Pauli)
    ax = _xyz_sym[x.axis]
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(x.name)}$(suffix)_{{$ax}}")
end

@latexrecipe function f(x::Spin)
    ax = _xyz_sym[x.axis]
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "{$(x.name)}$(suffix)_{{$ax}}")
end

@latexrecipe function f(x::Position)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "\\hat{$(x.name)}$(suffix)")
end

@latexrecipe function f(x::Momentum)
    suffix = _latex_index_suffix(x.index)
    return Expr(:latexifymerge, "\\hat{$(x.name)}$(suffix)")
end

# Extract a plain Julia number from CNum for LaTeX rendering
function _latex_prefactor(c::CNum)
    r_val = Symbolics.value(SymbolicUtils.unwrap(real(c)))
    i_val = Symbolics.value(SymbolicUtils.unwrap(imag(c)))
    if iszero(i_val)
        return r_val
    elseif iszero(r_val)
        return complex(zero(r_val), i_val)
    else
        return complex(r_val, i_val)
    end
end
_latex_prefactor(c::Number) = c

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
    parts = []
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

@latexrecipe function f(x::QAdd)
    st = sorted_arguments(x)
    if !isempty(x.indices)
        idx_parts = []
        for idx in x.indices
            r = Symbolics.value(SymbolicUtils.unwrap(idx.range))
            push!(idx_parts, "\\underset{$(idx.name)}{\\overset{$r}{\\sum}}")
        end
        if !isempty(x.non_equal)
            ne_str = join(["$(a.name){\\neq}$(b.name)" for (a, b) in x.non_equal], ",")
            idx_parts[end] = replace(idx_parts[end], "$(x.indices[end].name)" => "$(x.indices[end].name){\\neq}$(ne_str)")
        end
        prefix = join(idx_parts, " ")
        # Split terms into index-dependent and index-independent
        dep_terms = []
        indep_terms = []
        for t in st
            ops = first(keys(t.arguments))
            c = first(values(t.arguments))
            term_expr = _latex_term(c, ops)
            if any(idx -> _depends_on_index_term(c, ops, idx), x.indices)
                push!(dep_terms, term_expr)
            else
                push!(indep_terms, term_expr)
            end
        end
        sum_expr = if length(dep_terms) == 1
            Expr(:latexifymerge, prefix, dep_terms[1])
        else
            Expr(:latexifymerge, prefix, Expr(:call, :+, dep_terms...))
        end
        if isempty(indep_terms)
            return sum_expr
        else
            return Expr(:call, :+, indep_terms..., sum_expr)
        end
    end
    terms = [_latex_term(first(values(t.arguments)), first(keys(t.arguments))) for t in st]
    return Expr(:call, :+, terms...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
