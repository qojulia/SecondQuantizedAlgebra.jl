"""
CNum helpers — commutative number type and canonical sorting for operator sequences.
"""

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

@inline function _iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, _NUM_ZERO)
end

@inline function _iszero_cnum(c::CNum)
    return _iszero_num(c.re) && _iszero_num(c.im)
end

# Short-circuit CNum arithmetic: most prefactors have zero imaginary part,
# so we avoid 4 Num multiplications and 2 Num subtractions in the common case.
@inline function _mul_cnum(a::CNum, b::CNum)::CNum
    ar, ai = a.re, a.im
    br, bi = b.re, b.im
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return Complex(ar * br, _NUM_ZERO)
    elseif ai_zero
        return Complex(ar * br, ar * bi)
    elseif bi_zero
        return Complex(ar * br, ai * br)
    else
        return Complex(ar * br - ai * bi, ar * bi + ai * br)
    end
end

@inline function _neg_cnum(a::CNum)::CNum
    return _mul_cnum(_CNUM_NEG1, a)
end

@inline function _add_cnum(a::CNum, b::CNum)::CNum
    ar, ai = a.re, a.im
    br, bi = b.re, b.im
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return Complex(ar + br, _NUM_ZERO)
    elseif ai_zero
        return Complex(ar + br, bi)
    elseif bi_zero
        return Complex(ar + br, ai)
    else
        return Complex(ar + br, ai + bi)
    end
end

_sort_key(op::QSym) = (op.space_index, op.index.name)
