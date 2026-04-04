# HilbertSpace
Base.show(io::IO, h::FockSpace) = write(io, "ℋ(", string(h.name), ")")
function Base.show(io::IO, h::ProductSpace)
    show(io, h.spaces[1])
    for i in 2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
    return
end

# Operators
Base.show(io::IO, x::Destroy) = write(io, string(x.name))
Base.show(io::IO, x::Create) = write(io, string(x.name), "†")

# QMul
function Base.show(io::IO, x::QMul)
    if isempty(x.args_nc)
        show(io, x.arg_c)
        return
    end
    c = x.arg_c
    if c isa Union{Integer,AbstractFloat,Rational} && c == -1
        write(io, "-")
    elseif !(c isa Union{Integer,AbstractFloat,Rational} && isone(c))
        _needs_parens = !(c isa Union{Integer,AbstractFloat,Rational,Complex})
        _needs_parens && write(io, "(")
        show(io, c)
        _needs_parens && write(io, ")")
        write(io, " * ")
    end
    show(io, x.args_nc[1])
    for i in 2:length(x.args_nc)
        write(io, " * ")
        show(io, x.args_nc[i])
    end
    return
end

# QAdd
function Base.show(io::IO, x::QAdd)
    isempty(x.arguments) && return write(io, "0")
    show(io, x.arguments[1])
    for i in 2:length(x.arguments)
        term = x.arguments[i]
        if term.arg_c isa Union{Integer,AbstractFloat,Rational} && term.arg_c < 0
            write(io, " - ")
            show(io, QMul(-term.arg_c, term.args_nc))
        else
            write(io, " + ")
            show(io, term)
        end
    end
    return
end
