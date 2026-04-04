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

@latexrecipe function f(x::Destroy)
    return Expr(:latexifymerge, x.name)
end

@latexrecipe function f(x::Create)
    return Expr(:latexifymerge, "$(x.name)^{\\dagger}")
end

@latexrecipe function f(x::Transition)
    return Expr(:latexifymerge, "{$(x.name)}$(transition_idx_script[]){{$(x.i)$(x.j)}}")
end

@latexrecipe function f(x::Pauli)
    ax = _xyz_sym[x.axis]
    return Expr(:latexifymerge, "{$(x.name)}_{{$ax}}")
end

@latexrecipe function f(x::Spin)
    ax = _xyz_sym[x.axis]
    return Expr(:latexifymerge, "{$(x.name)}_{{$ax}}")
end

@latexrecipe function f(x::Position)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end

@latexrecipe function f(x::Momentum)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end

# Strip zero imaginary part from Complex for clean display
_latex_prefactor(c::Complex) = iszero(imag(c)) ? real(c) : c
_latex_prefactor(c) = c

@latexrecipe function f(x::QMul)
    c = _latex_prefactor(x.arg_c)
    if isempty(x.args_nc)
        return c
    end
    parts = []
    if c == -1
        push!(parts, :(-))
    elseif !isone(c)
        push!(parts, c)
    end
    for op in x.args_nc
        push!(parts, op)
    end
    return Expr(:latexifymerge, parts...)
end

@latexrecipe function f(x::QAdd)
    return Expr(:call, :+, x.arguments...)
end

const QLaTeX = Union{<:QField}
Base.show(io::IO, ::MIME"text/latex", x::QLaTeX) = write(io, latexify(x))
