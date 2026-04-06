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
    else
        push!(parts, pf)
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
            r = Symbolics.value(Symbolics.unwrap(idx.range))
            push!(idx_parts, "\\underset{$(idx.name)}{\\overset{$r}{\\sum}}")
        end
        if !isempty(x.non_equal)
            ne_str = join(["$(a.name){\\neq}$(b.name)" for (a, b) in x.non_equal], ",")
            idx_parts[end] = replace(idx_parts[end], "$(x.indices[end].name)" => "$(x.indices[end].name){\\neq}$(ne_str)")
        end
        prefix = join(idx_parts, " ")
        terms = [_latex_term(first(values(t.arguments)), first(keys(t.arguments))) for t in st]
        return Expr(:latexifymerge, prefix, " ", Expr(:call, :+, terms...))
    end
    terms = [_latex_term(first(values(t.arguments)), first(keys(t.arguments))) for t in st]
    return Expr(:call, :+, terms...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
