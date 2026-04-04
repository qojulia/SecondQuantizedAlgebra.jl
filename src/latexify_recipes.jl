@latexrecipe function f(x::Destroy)
    return Expr(:latexifymerge, x.name)
end

@latexrecipe function f(x::Create)
    return Expr(:latexifymerge, "$(x.name)^{\\dagger}")
end

@latexrecipe function f(x::QMul)
    if isempty(x.args_nc)
        return x.arg_c
    end
    parts = []
    if x.arg_c == -1
        push!(parts, :(-))
    elseif !isone(x.arg_c)
        push!(parts, x.arg_c)
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
