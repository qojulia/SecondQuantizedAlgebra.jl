"""
    QMul <: QTerm

Lazy product of quantum operators with a commutative prefactor.
`*` never applies commutation relations — it just collects operators
in canonical order and multiplies prefactors.

Fields:
- `arg_c::Num` — commutative (c-number) prefactor
- `args_nc::Vector{QSym}` — non-commutative operator factors, canonically sorted
"""

# Cached constants — avoids repeated Num() allocations on the hot path
const _NUM_ZERO = Num(0)
const _NUM_ONE = Num(1)
const _CNUM_ZERO = Complex(_NUM_ONE - _NUM_ONE, _NUM_ONE - _NUM_ONE)
const _CNUM_ONE = Complex(_NUM_ONE, _NUM_ZERO)
const _CNUM_NEG1 = Complex(-_NUM_ONE, _NUM_ZERO)
const _CNUM_IM = Complex(_NUM_ZERO, _NUM_ONE)
const _CNUM_NEG_IM = Complex(_NUM_ZERO, -_NUM_ONE)

_to_cnum(x::Complex{Num}) = x
_to_cnum(x::Num) = Complex(x, _NUM_ZERO)
function _to_cnum(x::Real)
    x == 0 && return _CNUM_ZERO
    x == 1 && return _CNUM_ONE
    x == -1 && return _CNUM_NEG1
    return Complex(Num(x), _NUM_ZERO)
end
function _to_cnum(x::Complex)
    x == im && return _CNUM_IM
    x == -im && return _CNUM_NEG_IM
    return Complex(Num(real(x)), Num(imag(x)))
end
function _to_cnum(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === complex
        args = SymbolicUtils.arguments(x)
        return Complex(Num(args[1]), Num(args[2]))
    end
    return Complex(Num(x), _NUM_ZERO)
end

function _iszero_cnum(c::CNum)
    r = Symbolics.unwrap(real(c))
    i = Symbolics.unwrap(imag(c))
    # Fast path: plain numbers avoid SymbolicUtils structural comparison
    r isa Number && i isa Number && return iszero(r) && iszero(i)
    return isequal(real(c), _NUM_ZERO) && isequal(imag(c), _NUM_ZERO)
end


_sort_key(op::QSym) = (op.space_index, op.copy_index, op.index.name)

# Group operators by site. Must be stable: same-site operators preserve
# their left-to-right multiplication order (they don't commute).
_site_sort!(v::Vector{QSym}) = sort!(v; by = _sort_key, alg = Base.MergeSort)

struct QMul <: QTerm
    arg_c::CNum
    args_nc::Vector{QSym}
    function QMul(arg_c::CNum, args_nc::Vector{QSym})
        return new(arg_c, args_nc)
    end
end
QMul(arg_c::Number, args_nc::Vector{QSym}) = QMul(_to_cnum(arg_c), args_nc)
QMul(args_nc::Vector{QSym}) = QMul(_to_cnum(1), args_nc)

"""
    prefactor(m::QMul) -> CNum

The commutative (c-number) prefactor of the product.
"""
prefactor(m::QMul) = m.arg_c

"""
    operators(m::QMul) -> Vector{QSym}

The non-commutative operator factors, in canonical (site-sorted) order.
"""
operators(m::QMul) = m.args_nc

Base.length(a::QMul) = length(a.args_nc)
Base.iszero(a::QMul) = iszero(a.arg_c)
Base.zero(::QMul) = QMul(_to_cnum(0), QSym[])
Base.zero(::Type{QMul}) = QMul(_to_cnum(0), QSym[])

# Equality and hashing
function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc) == length(b.args_nc) || return false
    for (x, y) in zip(a.args_nc, b.args_nc)
        isequal(x, y) || return false
    end
    return true
end
Base.:(==)(a::QMul, b::QMul) = isequal(a, b)
Base.hash(q::QMul, h::UInt) = hash(:QMul, hash(q.arg_c, hash(q.args_nc, h)))

# Adjoint: (ABC)† = C†B†A† — reverse and adjoint each factor.
function Base.adjoint(q::QMul)
    args_nc = QSym[adjoint(op) for op in reverse(q.args_nc)]
    _site_sort!(args_nc)
    return QMul(conj(q.arg_c), args_nc)
end

## Multiplication — always returns QMul

# QSym * QSym → QMul
function Base.:*(a::QSym, b::QSym)
    args = QSym[a, b]
    _site_sort!(args)
    return QMul(_to_cnum(1), args)
end

# QSym * Number → QMul
Base.:*(a::QSym, b::Number) = QMul(_to_cnum(b), QSym[a])
Base.:*(b::Number, a::QSym) = QMul(_to_cnum(b), QSym[a])

# QMul * Number → QMul
Base.:*(a::QMul, b::Number) = QMul(a.arg_c * _to_cnum(b), a.args_nc)
Base.:*(b::Number, a::QMul) = a * b

# QSym * QMul → QMul
function Base.:*(a::QSym, b::QMul)
    na = length(b.args_nc)
    args_nc = Vector{QSym}(undef, 1 + na)
    args_nc[1] = a
    copyto!(args_nc, 2, b.args_nc, 1, na)
    _site_sort!(args_nc)
    return QMul(b.arg_c, args_nc)
end
function Base.:*(a::QMul, b::QSym)
    na = length(a.args_nc)
    args_nc = Vector{QSym}(undef, na + 1)
    copyto!(args_nc, 1, a.args_nc, 1, na)
    args_nc[na + 1] = b
    _site_sort!(args_nc)
    return QMul(a.arg_c, args_nc)
end

# QMul * QMul → QMul
function Base.:*(a::QMul, b::QMul)
    na, nb = length(a.args_nc), length(b.args_nc)
    args_nc = Vector{QSym}(undef, na + nb)
    copyto!(args_nc, 1, a.args_nc, 1, na)
    copyto!(args_nc, na + 1, b.args_nc, 1, nb)
    _site_sort!(args_nc)
    return QMul(a.arg_c * b.arg_c, args_nc)
end

# Division
Base.:/(a::QSym, b::Number) = QMul(_to_cnum(1) / _to_cnum(b), QSym[a])
Base.:/(a::QMul, b::Number) = QMul(a.arg_c / _to_cnum(b), a.args_nc)

# Power
function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return QMul(Num(1), QSym[])
    args_nc = QSym[a for _ in 1:n]
    return QMul(_to_cnum(1), args_nc)
end
function Base.:^(a::QMul, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return QMul(Num(1), QSym[])
    args_nc = repeat(a.args_nc, n)
    _site_sort!(args_nc)
    return QMul(a.arg_c^n, args_nc)
end

# Negation
Base.:-(a::QSym) = QMul(_to_cnum(-1), QSym[a])
Base.:-(a::QMul) = QMul(-a.arg_c, a.args_nc)
