"""
Coefficient representation and canonical sorting for operator sequences.

A `Coeff` is the prefactor stored against each [`QTerm`](@ref). It carries a
native `ComplexF64` fast path for concrete numbers and falls back to a
`Complex{Num}` tail only for genuinely symbolic coefficients, so the common
numeric arithmetic never touches SymbolicUtils hashconsing.
"""

struct Coeff
    z::ComplexF64
    tail::Union{Nothing, Complex{Num}}
end
const CNum = Coeff

const _NUM_ZERO = Num(0)
const _NUM_ONE = Num(1)

# Adding 0.0+0.0im normalizes any signed zero (`-0.0 -> 0.0`) so that structurally
# equal coefficients (e.g. `conj(2)` vs `2`) stay `isequal` and hash identically.
@inline _native(z::ComplexF64) = Coeff(z + complex(0.0, 0.0), nothing)
@inline _symbolic(c::Complex{Num}) = Coeff(zero(ComplexF64), c)
@inline _is_native(c::Coeff) = c.tail === nothing

const _CNUM_ZERO = _native(zero(ComplexF64))
const _CNUM_ONE = _native(one(ComplexF64))
const _CNUM_NEG1 = _native(-one(ComplexF64))
const _CNUM_IM = _native(ComplexF64(0, 1))
const _CNUM_NEG_IM = _native(ComplexF64(0, -1))

# Only take the native path when the value round-trips through ComplexF64 with no
# loss, so exact rationals / bignums stay symbolic instead of silently truncating.
@inline function _native_complex(x::Number)
    z = ComplexF64(x)
    return z == x ? z : nothing
end

# Concrete numeric content of an (unwrapped) symbolic value, else `nothing`.
# `substitute`/`simplify` yield `BasicSymbolic` *constants*, not plain numbers,
# so both must be recognized to keep numeric coefficients on the native path.
@inline function _const_value(v)
    v isa Number && return v
    (v isa SymbolicUtils.BasicSymbolic && SymbolicUtils.isconst(v)) && return v.val
    return nothing
end
@inline _numeric_value(x::Num) = _const_value(SymbolicUtils.unwrap(x))

# Combine numeric real/imag parts as `rv + i*iv`. Uses the `Bool` imaginary unit
# so exact types (Int, Rational) are preserved, and tolerates an already-complex
# `rv`/`iv` (a complex value substituted into the real slot, e.g. `g => 2+3im`).
@inline _combine_reim(rv::Number, iv::Number) = rv + Complex(false, true) * iv

# Recover the smallest faithful `Num` for display/materialization: integer-valued
# floats print as integers, matching the pre-native coefficient behaviour.
@inline function _num_from_float(x::Float64)
    return (isinteger(x) && abs(x) <= 9.007199254740992e15) ? Num(Int(x)) : Num(x)
end

_to_cnum(x::Coeff) = x
_to_cnum(x::Num) = _to_cnum(SymbolicUtils.unwrap(x))
function _to_cnum(x::Real)
    z = _native_complex(x)
    return z === nothing ? _symbolic(Complex(Num(x), _NUM_ZERO)) : _native(z)
end
function _to_cnum(x::Complex)
    z = _native_complex(x)
    return z === nothing ? _symbolic(Complex(Num(real(x)), Num(imag(x)))) : _native(z)
end
function _to_cnum(x::Complex{Num})
    rv = _numeric_value(real(x))
    iv = _numeric_value(imag(x))
    if rv !== nothing && iv !== nothing
        z = _native_complex(_combine_reim(rv, iv))
        z === nothing || return _native(z)
    end
    return _symbolic(x)
end
function _to_cnum(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === complex
        args = SymbolicUtils.arguments(x)
        return _to_cnum(Complex(Num(args[1]), Num(args[2])))
    end
    SymbolicUtils.isconst(x) && return _to_cnum(x.val)
    return _symbolic(Complex(Num(x), _NUM_ZERO))
end

# Rebuild a coefficient from real/imag `Num` parts, folding back to native when
# both parts are concrete numbers. Used by the symbolic boundary functions.
function _cnum(re::Num, im::Num)
    rv = _numeric_value(re)
    iv = _numeric_value(im)
    if rv !== nothing && iv !== nothing
        z = _native_complex(_combine_reim(rv, iv))
        z === nothing || return _native(z)
    end
    return _symbolic(Complex(re, im))
end

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = _to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = _to_cnum(x)

Base.real(c::Coeff) = c.tail === nothing ? _num_from_float(real(c.z)) : real(c.tail)
Base.imag(c::Coeff) = c.tail === nothing ? _num_from_float(imag(c.z)) : imag(c.tail)

# Materialize the full `Complex{Num}`; the only place a coefficient lowers back to
# Symbolics for symbolic boundaries (substitute / simplify / expand / printing).
to_num(c::Coeff) = c.tail === nothing ? Complex(real(c), imag(c)) : c.tail

Base.show(io::IO, c::Coeff) = show(io, to_num(c))

function Base.isequal(a::Coeff, b::Coeff)
    an = a.tail === nothing
    bn = b.tail === nothing
    an && bn && return isequal(a.z, b.z)
    (an || bn) && return false
    return isequal(a.tail, b.tail)
end
Base.:(==)(a::Coeff, b::Coeff) = isequal(a, b)
function Base.hash(c::Coeff, h::UInt)
    return c.tail === nothing ? hash(c.z, hash(:CoeffNative, h)) : hash(c.tail, hash(:CoeffSym, h))
end

# Coefficients are routinely compared against plain numbers / `Complex{Num}`
# (e.g. `prefactor(x) == 2`, `q[key] == _CNUM_ONE`); promote the number side.
Base.isequal(a::Coeff, b::Number) = isequal(a, _to_cnum(b))
Base.isequal(a::Number, b::Coeff) = isequal(_to_cnum(a), b)
Base.:(==)(a::Coeff, b::Number) = isequal(a, _to_cnum(b))
Base.:(==)(a::Number, b::Coeff) = isequal(_to_cnum(a), b)

Base.iszero(c::Coeff) = _iszero_cnum(c)
Base.conj(c::Coeff) = _conj_cnum(c)

# Coefficients flow through downstream code (and tests) as numbers; support the
# usual scalar arithmetic, promoting any `Number` operand into a `Coeff` first.
Base.:-(c::Coeff) = _neg_cnum(c)
Base.:+(a::Coeff, b::Coeff) = _add_cnum(a, b)
Base.:+(a::Coeff, b::Number) = _add_cnum(a, _to_cnum(b))
Base.:+(a::Number, b::Coeff) = _add_cnum(_to_cnum(a), b)
Base.:-(a::Coeff, b::Coeff) = _add_cnum(a, _neg_cnum(b))
Base.:-(a::Coeff, b::Number) = _add_cnum(a, _neg_cnum(_to_cnum(b)))
Base.:-(a::Number, b::Coeff) = _add_cnum(_to_cnum(a), _neg_cnum(b))
Base.:*(a::Coeff, b::Coeff) = _mul_cnum(a, b)
Base.:*(a::Coeff, b::Number) = _mul_cnum(a, _to_cnum(b))
Base.:*(a::Number, b::Coeff) = _mul_cnum(_to_cnum(a), b)
function Base.:/(a::Coeff, b::Coeff)
    (a.tail === nothing && b.tail === nothing) && return _native(a.z / b.z)
    return _to_cnum(to_num(a) / to_num(b))
end
Base.:/(a::Coeff, b::Number) = a / _to_cnum(b)
Base.:/(a::Number, b::Coeff) = _to_cnum(a) / b

_sym_conj(x::Num) = SymbolicUtils.symtype(x) <: Real ? x : Num(conj(SymbolicUtils.unwrap(x)))

@inline function _conj_cnum(c::Coeff)
    c.tail === nothing && return _native(conj(c.z))
    # conj(re + i*im) = conj(re) - i*conj(im); each part may be a Number-symtype symbol.
    return _cnum(_sym_conj(real(c.tail)), -_sym_conj(imag(c.tail)))
end

@inline function _iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, _NUM_ZERO)
end

@inline function _iszero_cnum(c::Coeff)
    c.tail === nothing && return iszero(c.z)
    return _iszero_num(real(c.tail)) && _iszero_num(imag(c.tail))
end

# `unwrap` returns a `BasicSymbolic` even for numeric constants, so test the
# node kind rather than `isa Number`.
@inline function _is_symbolic_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    return SymbolicUtils.issym(v) || SymbolicUtils.iscall(v)
end

@inline function _is_symbolic_cnum(c::Coeff)
    c.tail === nothing && return false
    return _is_symbolic_num(real(c.tail)) || _is_symbolic_num(imag(c.tail))
end

# Structural `a == -b`; see `_addto_key!` for why this is needed.
@inline function _isneg_cnum(a::Coeff, b::Coeff)
    an = a.tail === nothing
    bn = b.tail === nothing
    an && bn && return isequal(a.z, -b.z)
    (an || bn) && return false
    return isequal(real(a.tail), -real(b.tail)) && isequal(imag(a.tail), -imag(b.tail))
end

@inline function _mul_cnum(a::Coeff, b::Coeff)
    (a.tail === nothing && b.tail === nothing) && return _native(a.z * b.z)
    return _mul_cnum_slow(a, b)
end

# Short-circuit the symbolic multiply: most prefactors have zero imaginary part,
# so we avoid 4 Num multiplications and 2 Num subtractions in the common case.
@noinline function _mul_cnum_slow(a::Coeff, b::Coeff)
    ca, cb = to_num(a), to_num(b)
    ar, ai = real(ca), imag(ca)
    br, bi = real(cb), imag(cb)
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return _cnum(ar * br, _NUM_ZERO)
    elseif ai_zero
        return _cnum(ar * br, ar * bi)
    elseif bi_zero
        return _cnum(ar * br, ai * br)
    else
        return _cnum(ar * br - ai * bi, ar * bi + ai * br)
    end
end

@inline function _neg_cnum(a::Coeff)
    a.tail === nothing && return _native(-a.z)
    return _cnum(-real(a.tail), -imag(a.tail))
end

@inline function _add_cnum(a::Coeff, b::Coeff)
    (a.tail === nothing && b.tail === nothing) && return _native(a.z + b.z)
    return _add_cnum_slow(a, b)
end

@noinline function _add_cnum_slow(a::Coeff, b::Coeff)
    ca, cb = to_num(a), to_num(b)
    ar, ai = real(ca), imag(ca)
    br, bi = real(cb), imag(cb)
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return _cnum(ar + br, _NUM_ZERO)
    elseif ai_zero
        return _cnum(ar + br, bi)
    elseif bi_zero
        return _cnum(ar + br, ai)
    else
        return _cnum(ar + br, ai + bi)
    end
end

_sort_key(op::QSym) = (op.space_index, op.index.name)
