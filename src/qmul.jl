"""
    QMul <: QTerm

Lazy product of quantum operators with a commutative prefactor.
`*` never applies commutation relations — it just collects operators
in canonical order and multiplies prefactors.

Fields:
- `arg_c::Num` — commutative (c-number) prefactor
- `args_nc::Vector{QSym}` — non-commutative operator factors, canonically sorted
"""
const CNum = Complex{Num}

_to_cnum(x::Complex{Num}) = x
_to_cnum(x::Num) = Complex(x, Num(0))
_to_cnum(x::Real) = Complex(Num(x), Num(0))
_to_cnum(x::Complex) = Complex(Num(real(x)), Num(imag(x)))

_sort_key(op::QSym) = (op.space_index, op.copy_index, op.index.name)

struct QMul <: QTerm
    arg_c::CNum
    args_nc::Vector{QSym}
    function QMul(arg_c::CNum, args_nc::Vector{QSym})
        return new(arg_c, args_nc)
    end
end
QMul(arg_c::Number, args_nc::Vector{QSym}) = QMul(_to_cnum(arg_c), args_nc)
QMul(args_nc::Vector{QSym}) = QMul(_to_cnum(1), args_nc)

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
    sort!(args_nc; by = _sort_key)
    return QMul(conj(q.arg_c), args_nc)
end

## Multiplication — always returns QMul

# QSym * QSym → QMul
function Base.:*(a::QSym, b::QSym)
    args = QSym[a, b]
    sort!(args; by = _sort_key)
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
    sort!(args_nc; by = _sort_key)
    return QMul(b.arg_c, args_nc)
end
function Base.:*(a::QMul, b::QSym)
    na = length(a.args_nc)
    args_nc = Vector{QSym}(undef, na + 1)
    copyto!(args_nc, 1, a.args_nc, 1, na)
    args_nc[na + 1] = b
    sort!(args_nc; by = _sort_key)
    return QMul(a.arg_c, args_nc)
end

# QMul * QMul → QMul
function Base.:*(a::QMul, b::QMul)
    na, nb = length(a.args_nc), length(b.args_nc)
    args_nc = Vector{QSym}(undef, na + nb)
    copyto!(args_nc, 1, a.args_nc, 1, na)
    copyto!(args_nc, na + 1, b.args_nc, 1, nb)
    sort!(args_nc; by = _sort_key)
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
    sort!(args_nc; by = _sort_key)
    return QMul(a.arg_c^n, args_nc)
end

# Negation
Base.:-(a::QSym) = QMul(_to_cnum(-1), QSym[a])
Base.:-(a::QMul) = QMul(-a.arg_c, a.args_nc)
