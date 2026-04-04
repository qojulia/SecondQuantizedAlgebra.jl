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

# Hilbert spaces — new types
Base.show(io::IO, h::NLevelSpace) = write(io, "ℋ(", string(h.name), ")")
Base.show(io::IO, h::PauliSpace) = write(io, "ℋ(", string(h.name), ")")
Base.show(io::IO, h::SpinSpace) = write(io, "ℋ(", string(h.name), ")")

# Operators — Fock
Base.show(io::IO, x::Destroy) = write(io, string(x.name))
Base.show(io::IO, x::Create) = write(io, string(x.name), "†")

# Operators — Transition: σ₁₂
const _subscript_digits = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
function _to_subscript(n::Int)
    return join(_subscript_digits[d + 1] for d in reverse(digits(n)))
end
function Base.show(io::IO, x::Transition)
    write(io, string(x.name), _to_subscript(x.i), _to_subscript(x.j))
end

# Operators — Pauli/Spin: σx, Sz
const _xyz_sym = [:x, :y, :z]
Base.show(io::IO, x::Pauli) = write(io, string(x.name), string(_xyz_sym[x.axis]))
Base.show(io::IO, x::Spin) = write(io, string(x.name), string(_xyz_sym[x.axis]))

# Helper: clean display of a numeric prefactor
function _show_prefactor(io::IO, c)
    if c isa Complex
        r, i = reim(c)
        if iszero(i)
            _show_prefactor(io, r)
        elseif iszero(r)
            if isone(i)
                write(io, "im")
            elseif i == -1
                write(io, "-im")
            else
                show(io, i)
                write(io, "im")
            end
        else
            write(io, "(")
            show(io, c)
            write(io, ")")
        end
    elseif !(c isa Union{Integer,AbstractFloat,Rational})
        write(io, "(")
        show(io, c)
        write(io, ")")
    else
        show(io, c)
    end
end

# Check if a prefactor is effectively 1
_is_unit(c::Union{Integer,AbstractFloat,Rational}) = isone(c)
_is_unit(c::Complex) = isone(real(c)) && iszero(imag(c))
_is_unit(::Any) = false

# Check if a prefactor is effectively -1
_is_neg_unit(c::Union{Integer,AbstractFloat,Rational}) = c == -1
_is_neg_unit(c::Complex) = real(c) == -1 && iszero(imag(c))
_is_neg_unit(::Any) = false

# Check if a prefactor is a real negative number (safe for < comparison)
_is_real_negative(c::Union{Integer,AbstractFloat,Rational}) = c < 0
_is_real_negative(::Any) = false

# QMul
function Base.show(io::IO, x::QMul)
    if isempty(x.args_nc)
        show(io, x.arg_c)
        return
    end
    c = x.arg_c
    if _is_neg_unit(c)
        write(io, "-")
    elseif !_is_unit(c)
        _show_prefactor(io, c)
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
        if _is_real_negative(term.arg_c)
            write(io, " - ")
            show(io, QMul(-term.arg_c, term.args_nc))
        else
            write(io, " + ")
            show(io, term)
        end
    end
    return
end
