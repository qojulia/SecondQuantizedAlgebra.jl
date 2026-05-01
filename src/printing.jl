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

Base.show(io::IO, h::NLevelSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PauliSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::SpinSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PhaseSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))

# Index suffix helper
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
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Create)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return write(io, "'")
end

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
    _show_index_suffix(io, x.index)
    _write_subscript(io, x.i)
    return _write_subscript(io, x.j)
end

const _xyz_sym = (:x, :y, :z)
function Base.show(io::IO, x::Pauli)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end
function Base.show(io::IO, x::Spin)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end

function Base.show(io::IO, x::Position)
    print(io, x.name)
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Momentum)
    print(io, x.name)
    return _show_index_suffix(io, x.index)
end

# Prefactor display helpers
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
_show_prefactor(io::IO, c::Number) = print(io, c)

function _is_unit(c::CNum)
    return iszero(imag(c)) && isone(real(c))
end
_is_unit(c::Number) = isone(c)

function _is_neg_unit(c::CNum)
    return iszero(imag(c)) && isequal(real(c), Num(-1))
end
_is_neg_unit(c::Number) = isequal(c, -1)

function _is_real_negative(c::CNum)
    if iszero(imag(c))
        r = SymbolicUtils.unwrap(real(c))
        return r isa Real && r < 0
    end
    return false
end
_is_real_negative(::Number) = false

# Check if a symbolic prefactor needs parentheses when followed by operators.
function _needs_pf_parens(c::CNum)
    iszero(imag(c)) || return false
    r = Symbolics.value(SymbolicUtils.unwrap(real(c)))
    r isa Number && return false
    SymbolicUtils.iscall(r) || return false
    op = SymbolicUtils.operation(r)
    return op === (/) || op === (+)
end

# Show a single term: prefactor * op1 * op2 * ...
function _show_term(io::IO, c::CNum, ops::Vector{QSym})
    if isempty(ops)
        _show_prefactor(io, c)
        return
    end
    if _is_neg_unit(c)
        write(io, "-")
    elseif !_is_unit(c)
        if _needs_pf_parens(c)
            write(io, "(")
            _show_prefactor(io, c)
            write(io, ")")
        else
            _show_prefactor(io, c)
        end
        write(io, " * ")
    end
    show(io, ops[1])
    for i in 2:length(ops)
        write(io, " * ")
        show(io, ops[i])
    end
    return
end

# Show a list of terms joined by + / -
function _show_terms(io::IO, st::Vector{QAdd})
    isempty(st) && return write(io, "0")
    ops1, c1 = first(keys(st[1].arguments)), first(values(st[1].arguments))
    _show_term(io, c1, ops1)
    for i in 2:length(st)
        ops_i = first(keys(st[i].arguments))
        c_i = first(values(st[i].arguments))
        if _is_real_negative(c_i)
            write(io, " - ")
            _show_term(io, -c_i, ops_i)
        else
            write(io, " + ")
            _show_term(io, c_i, ops_i)
        end
    end
    return
end

# QAdd
function Base.show(io::IO, x::QAdd)
    st = sorted_arguments(x)
    isempty(st) && return write(io, "0")

    if !isempty(x.indices)
        # Split into index-dependent and index-independent terms
        dep = eltype(st)[]
        indep = eltype(st)[]
        for t in st
            ops = first(keys(t.arguments))
            c = first(values(t.arguments))
            if any(idx -> _depends_on_index_term(c, ops, idx), x.indices)
                push!(dep, t)
            else
                push!(indep, t)
            end
        end
        # Print index-independent terms first
        if !isempty(indep)
            _show_terms(io, indep)
            if !isempty(dep)
                write(io, " + ")
            end
        end
        # Print summation with dependent terms
        if !isempty(dep)
            for idx in x.indices
                write(io, "Σ(")
                print(io, idx.name)
                write(io, "=1:")
                print(io, idx.range)
                write(io, ")")
            end
            if !isempty(x.non_equal)
                write(io, "(")
                for (k, (a, b)) in enumerate(x.non_equal)
                    k > 1 && write(io, ",")
                    print(io, a.name)
                    write(io, "≠")
                    print(io, b.name)
                end
                write(io, ")")
            end
            write(io, " ")
            if length(dep) > 1
                write(io, "(")
                _show_terms(io, dep)
                write(io, ")")
            else
                _show_terms(io, dep)
            end
        end
    else
        _show_terms(io, st)
    end
    return
end
