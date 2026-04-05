# HilbertSpace
function Base.show(io::IO, h::FockSpace)
    write(io, "ℋ(")
    print(io, h.name)
    return write(io, ")")
end
function Base.show(io::IO, h::ProductSpace)
    show(io, h.spaces[1])
    for i in 2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
    return
end

# Hilbert spaces — other types
Base.show(io::IO, h::NLevelSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PauliSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::SpinSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PhaseSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
function Base.show(io::IO, h::ClusterSpace)
    write(io, "Cluster(")
    show(io, h.original_space)
    write(io, ", N=")
    print(io, h.N)
    write(io, ", order=")
    print(io, h.order)
    return write(io, ")")
end

# Index suffix helpers
function _show_copy_suffix(io::IO, ci::Int)
    if ci > 1
        write(io, "_")
        print(io, ci)
    end
    return
end

function _show_index_suffix(io::IO, idx::Index)
    if has_index(idx)
        write(io, "_")
        print(io, idx.name)
    end
    return
end

# Operators — Fock
function Base.show(io::IO, x::Destroy)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Create)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    _show_index_suffix(io, x.index)
    return write(io, "†")
end

# Operators — Transition: σ₁₂ or σ_i₁₂
const _subscript_digits = ('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')
function _write_subscript(io::IO, n::Int)
    return if 0 <= n <= 9
        write(io, _subscript_digits[n + 1])
    else
        for d in reverse(digits(n))
            write(io, _subscript_digits[d + 1])
        end
    end
end
function Base.show(io::IO, x::Transition)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    _show_index_suffix(io, x.index)
    _write_subscript(io, x.i)
    return _write_subscript(io, x.j)
end

# Operators — Pauli/Spin: σx, Sz
const _xyz_sym = (:x, :y, :z)
function Base.show(io::IO, x::Pauli)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end
function Base.show(io::IO, x::Spin)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end

# Operators — Position/Momentum
function Base.show(io::IO, x::Position)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Momentum)
    print(io, x.name)
    _show_copy_suffix(io, x.copy_index)
    return _show_index_suffix(io, x.index)
end

# Helper: clean display of a numeric prefactor
function _show_prefactor(io::IO, c::CNum)
    return if iszero(imag(c))
        print(io, real(c))
    elseif iszero(real(c))
        i = imag(c)
        if isone(i)
            write(io, "im")
        elseif isequal(i, Num(-1))
            write(io, "-im")
        else
            print(io, i)
            write(io, "im")
        end
    else
        write(io, "(")
        print(io, real(c))
        write(io, " + ")
        print(io, imag(c))
        write(io, "im)")
    end
end
_show_prefactor(io::IO, c) = print(io, c)

# Check if a CNum prefactor is effectively 1 (real part 1, imag part 0)
function _is_unit(c::CNum)
    return iszero(imag(c)) && isone(real(c))
end
_is_unit(c) = isone(c)

# Check if a CNum prefactor is effectively -1
function _is_neg_unit(c::CNum)
    return iszero(imag(c)) && isequal(real(c), Num(-1))
end
_is_neg_unit(c) = isequal(c, -1)

# Check if a prefactor is a real negative number
function _is_real_negative(c::CNum)
    if iszero(imag(c))
        r = Symbolics.unwrap(real(c))
        return r isa Real && r < 0
    end
    return false
end
_is_real_negative(::Any) = false

# QMul
function Base.show(io::IO, x::QMul)
    if isempty(x.args_nc)
        _show_prefactor(io, x.arg_c)
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
    # Show Σ prefix if this is a sum
    if !isempty(x.indices)
        for idx in x.indices
            write(io, "Σ($(idx.name)=1:")
            print(io, idx.range)
            write(io, ")")
        end
        if !isempty(x.non_equal)
            write(io, "(")
            for (k, (a, b)) in enumerate(x.non_equal)
                k > 1 && write(io, ",")
                write(io, "$(a.name)≠$(b.name)")
            end
            write(io, ")")
        end
        write(io, " ")
    end

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
