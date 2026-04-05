# Transition superscript toggle
const transition_idx_script = Ref(:^)

"""
    transition_superscript(x::Bool)

Toggle whether Transition level indices are printed as superscript (true)
or subscript (false) in LaTeX. Default is superscript.
"""
function transition_superscript(x::Bool)
    transition_idx_script[] = x ? :^ : :_
    return x
end

# Helper: format index suffix for LaTeX
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

# Extract a plain Julia number from CNum for comparison in LaTeX recipes
function _latex_prefactor(c::CNum)
    r_val = Symbolics.value(Symbolics.unwrap(real(c)))
    i_val = Symbolics.value(Symbolics.unwrap(imag(c)))
    if iszero(i_val)
        return r_val
    elseif iszero(r_val)
        return complex(zero(r_val), i_val)
    else
        return complex(r_val, i_val)
    end
end
_latex_prefactor(c) = c

@latexrecipe function f(x::QMul)
    c = _latex_prefactor(x.arg_c)
    if isempty(x.args_nc)
        return c
    end
    parts = []
    if c isa Number && c == -1
        push!(parts, :(-))
    elseif c isa Number && isone(c)
        # skip prefactor
    else
        push!(parts, c)
    end
    for op in x.args_nc
        push!(parts, op)
    end
    return Expr(:latexifymerge, parts...)
end

@latexrecipe function f(x::QAdd)
    if !isempty(x.indices)
        # Sum notation: Σ_{i=1}^{N} (terms)
        idx_parts = []
        for idx in x.indices
            r = Symbolics.value(Symbolics.unwrap(idx.range))
            push!(idx_parts, "\\underset{$(idx.name)}{\\overset{$r}{\\sum}}")
        end
        if !isempty(x.non_equal)
            ne_str = join(["$(a.name){\\ne}$(b.name)" for (a, b) in x.non_equal], ",")
            idx_parts[end] = replace(idx_parts[end], "$(x.indices[end].name)" => "$(x.indices[end].name){\\ne}$(ne_str)")
        end
        prefix = join(idx_parts, " ")
        # Wrap sum contents
        return Expr(:latexifymerge, prefix, " ", Expr(:call, :+, x.arguments...))
    end
    return Expr(:call, :+, x.arguments...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
