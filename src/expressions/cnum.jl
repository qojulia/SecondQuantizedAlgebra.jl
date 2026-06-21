"""
Coefficient representation and canonical sorting for operator sequences.

A `Coeff` is the prefactor stored against each [`QTerm`](@ref). It has three forms:
a native `ComplexF64` fast path for concrete numbers, a [`Poly`](@ref) sparse
parameter polynomial for products and sums of named parameters, and a
`Complex{Num}` fallback for genuinely non-polynomial coefficients (`exp`, `sin`,
divisions by sums, ...). Numeric and parameter-polynomial arithmetic stay native,
so the common coefficient never routes through SymbolicUtils hashconsing; the
coefficient lowers to `Complex{Num}` only at the symbolic boundaries (`to_num`).
"""

# Zero-size tag marking the native fast path (the value lives inline in `z`).
# A named singleton rather than `nothing`: same layout and cost, but the type
# says "native" instead of "absent".
struct Native end
const NATIVE = Native()

struct Coeff
    z::ComplexF64
    tail::Union{Native, Poly, Complex{Num}}
end
const CNum = Coeff

const _NUM_ZERO = Num(0)
const _NUM_ONE = Num(1)
const _ONE_C = ComplexF64(1)
const _EMPTY_SYMS = SymbolicUtils.BasicSymbolic[]
const _EMPTY_EXPS = Int[]

# Adding 0.0+0.0im normalizes any signed zero (`-0.0 -> 0.0`) so that structurally
# equal coefficients (e.g. `conj(2)` vs `2`) stay `isequal` and hash identically.
@inline _native(z::ComplexF64) = Coeff(z + complex(0.0, 0.0), NATIVE)
@inline _symbolic(c::Complex{Num}) = Coeff(zero(ComplexF64), c)
@inline _poly_coeff(p::Poly) = Coeff(zero(ComplexF64), p)
@inline _is_native(c::Coeff) = c.tail isa Native
@inline _is_poly(c::Coeff) = c.tail isa Poly

const _CNUM_ZERO = _native(zero(ComplexF64))
const _CNUM_ONE = _native(one(ComplexF64))
const _CNUM_NEG1 = _native(-one(ComplexF64))
const _CNUM_IM = _native(ComplexF64(0, 1))
const _CNUM_NEG_IM = _native(ComplexF64(0, -1))

# Concrete numeric content of an (unwrapped) symbolic value, else `nothing`.
# `substitute`/`simplify` yield `BasicSymbolic` *constants*, not plain numbers,
# so both must be recognized to keep numeric coefficients on the native path.
# (`nothing` is consumed by `isa` checks at the call sites, never threaded.)
@inline function _const_value(v)
    v isa Number && return v
    (v isa SymbolicUtils.BasicSymbolic && SymbolicUtils.isconst(v)) && return v.val
    return nothing
end
@inline _numeric_value(x::Num) = _const_value(SymbolicUtils.unwrap(x))

# Recover the smallest faithful `Num` for display/materialization: integer-valued
# floats print as integers, matching the pre-native coefficient behaviour.
@inline function _num_from_float(x::Float64)
    return (isinteger(x) && abs(x) <= 9.007199254740992e15) ? Num(Int(x)) : Num(x)
end

_to_cnum(x::Coeff) = x
_to_cnum(x::Num) = _rec(SymbolicUtils.unwrap(x))
# Native only when the value round-trips through ComplexF64 with no loss, so exact
# rationals / bignums stay symbolic instead of silently truncating.
function _to_cnum(x::Real)
    z = ComplexF64(x)
    return z == x ? _native(z) : _symbolic(Complex(Num(x), _NUM_ZERO))
end
function _to_cnum(x::Complex)
    z = ComplexF64(x)
    return z == x ? _native(z) : _symbolic(Complex(Num(real(x)), Num(imag(x))))
end
_to_cnum(x::Complex{Num}) = _cnum(real(x), imag(x))
_to_cnum(x::SymbolicUtils.BasicSymbolic) = _rec(x)

# Build a coefficient from real/imag `Num` parts by interpreting `re + im*i` in the
# Coeff algebra (see `_rec`). The canonicalizing constructor used by symbolic
# boundaries that may yield a polynomial (substitute / conj / change_index), so
# native / polynomial / symbolic forms are recovered uniformly.
_cnum(re::Num, im::Num) =
    _add_cnum(_rec(SymbolicUtils.unwrap(re)), _mul_cnum(_CNUM_IM, _rec(SymbolicUtils.unwrap(im))))

# Rebuild after a *materialized symbolic arithmetic* step (the `Complex{Num}`
# fallback of `_mul_cnum`/`_add_cnum`/`_neg_cnum`). One operand was already a
# genuinely non-polynomial symbolic tail, so the result is native (when it folds to
# a constant) or symbolic, never a fresh polynomial. It therefore must NOT re-enter
# `_rec`, which would re-decompose the product and recurse without end; it only
# folds a numeric constant back to the native path.
function _cnum_sym(re::Num, im::Num)
    rv = _numeric_value(re)
    iv = _numeric_value(im)
    if rv isa Number && iv isa Number
        w = rv + Complex(false, true) * iv
        z = ComplexF64(w)
        z == w && return _native(z)
    end
    return _symbolic(Complex(re, im))
end

# === Parameter-polynomial tier: folding, materialization, recognition ===

# Fold a canonical term list into a Coeff: empty -> zero, a single constant term
# -> native, otherwise a Poly. Producers (`_poly_add`/`_poly_mul`/`_poly_scale`)
# return canonical term lists, so no re-sort is needed here.
function _from_poly(terms::Vector{Monomial})
    isempty(terms) && return _CNUM_ZERO
    if length(terms) == 1 && isempty(terms[1].syms)
        return _to_cnum(terms[1].scalar)
    end
    return _poly_coeff(Poly(terms))
end

# One monomial term -> Complex{Num}. Powers are built by repeated `Num`
# multiplication/division (both infer to `Num`) rather than `Num ^ Int`, which
# SymbolicUtils infers as `Any` and would degrade `to_num` / `real` everywhere.
function _term_to_num(m::Monomial)
    prod = _NUM_ONE
    @inbounds for i in eachindex(m.syms)
        base = Num(m.syms[i])
        e = m.exps[i]
        if e >= 0
            for _ in 1:e
                prod = prod * base
            end
        else
            for _ in 1:(-e)
                prod = prod / base
            end
        end
    end
    sr = _num_from_float(real(m.scalar))
    si = _num_from_float(imag(m.scalar))
    return Complex(sr * prod, si * prod)
end

# Sum the terms; the only place a polynomial lowers to SymbolicUtils.
function _poly_to_num(p::Poly)
    isempty(p.terms) && return Complex(_NUM_ZERO, _NUM_ZERO)
    acc = _term_to_num(p.terms[1])
    @inbounds for i in 2:length(p.terms)
        acc = acc + _term_to_num(p.terms[i])
    end
    return acc
end

# An "atom" factor is a bare symbol, an array-variable index (`ω[i]`, a
# `getindex` call), or `conj(atom)`.
@inline function _is_atom(b)
    b isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.issym(b) && return true
    SymbolicUtils.iscall(b) || return false
    op = SymbolicUtils.operation(b)
    op === getindex && return true
    return op === conj && length(SymbolicUtils.arguments(b)) == 1 &&
        _is_atom(only(SymbolicUtils.arguments(b)))
end

# A bare atom (symbol / array index / `conj(atom)`) as a single-monomial Coeff.
@inline _atom_coeff(x::SymbolicUtils.BasicSymbolic) =
    _poly_coeff(Poly(Monomial[Monomial(_ONE_C, SymbolicUtils.BasicSymbolic[x], [1])]))
# An unrecognized symbolic value, kept on the `Complex{Num}` symbolic path.
@inline _sym_leaf(x::SymbolicUtils.BasicSymbolic) = _symbolic(Complex(Num(x), _NUM_ZERO))

# Total recognizer: interpret a symbolic value as a `Coeff` by evaluating its
# expression tree in the Coeff algebra. Numbers / atoms map to native / 1-term
# polynomial coefficients; `+`, `*`, integer `^`, and monomial `/` compose through
# the coefficient arithmetic (which keeps polynomials native and escalates to the
# `Complex{Num}` tail only when an operand is non-polynomial); everything else
# (`exp`, `sin`, symbolic powers, divisions by sums, ...) is a symbolic leaf.
# There is no failure sentinel: an unrecognized subterm is just a symbolic-tail
# `Coeff` that the arithmetic absorbs, so recognition is total and never `nothing`.
_rec(x::Number)::Coeff = _to_cnum(x)
_rec(x::Num)::Coeff = _rec(SymbolicUtils.unwrap(x))
_rec(x)::Coeff = _to_cnum(x)
function _rec(x::SymbolicUtils.BasicSymbolic)::Coeff
    SymbolicUtils.isconst(x) && return _rec(x.val)
    SymbolicUtils.issym(x) && return _atom_coeff(x)
    SymbolicUtils.iscall(x) || return _sym_leaf(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        c = _CNUM_ZERO
        for a in args
            c = _add_cnum(c, _rec(a))
        end
        return c
    elseif op === (*)
        c = _CNUM_ONE
        for a in args
            c = _mul_cnum(c, _rec(a))
        end
        return c
    elseif op === (^)
        length(args) == 2 || return _sym_leaf(x)
        pv = _const_value(args[2])
        pv isa Integer || return _sym_leaf(x)
        n = Int(pv)
        base = _rec(args[1])
        if n >= 0
            c = _CNUM_ONE
            for _ in 1:n
                c = _mul_cnum(c, base)
            end
            return c
        elseif _is_native(base)
            return _native(base.z^n)
        elseif base.tail isa Poly && length(base.tail.terms) == 1
            m = base.tail.terms[1]
            s = inv(m.scalar)
            s * m.scalar == _ONE_C &&
                return _poly_coeff(Poly(Monomial[Monomial(s^(-n), m.syms, Int[e * n for e in m.exps])]))
        end
        return _sym_leaf(x)
    elseif op === getindex
        return _atom_coeff(x)
    elseif op === conj
        return _is_atom(x) ? _atom_coeff(x) : _sym_leaf(x)
    elseif op === (/)
        length(args) == 2 || return _sym_leaf(x)
        den = _rec(args[2])
        if _is_native(den)
            return _mul_cnum(_rec(args[1]), _native(inv(den.z)))
        elseif den.tail isa Poly && length(den.tail.terms) == 1
            m = den.tail.terms[1]
            s = inv(m.scalar)
            s * m.scalar == _ONE_C &&
                return _mul_cnum(
                _rec(args[1]),
                _poly_coeff(Poly(Monomial[Monomial(s, m.syms, Int[-e for e in m.exps])]))
            )
        end
        return _sym_leaf(x)
    elseif op === complex
        length(args) == 2 && return _cnum(Num(args[1]), Num(args[2]))
        return _sym_leaf(x)
    end
    return _sym_leaf(x)
end

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = _to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = _to_cnum(x)

Base.real(c::Coeff) = _is_native(c) ? _num_from_float(real(c.z)) : real(to_num(c))
Base.imag(c::Coeff) = _is_native(c) ? _num_from_float(imag(c.z)) : imag(to_num(c))

# Materialize the full `Complex{Num}`; the only place a coefficient lowers back to
# Symbolics for symbolic boundaries (substitute / simplify / expand / printing).
function to_num(c::Coeff)
    t = c.tail
    t isa Native && return Complex(_num_from_float(real(c.z)), _num_from_float(imag(c.z)))
    t isa Poly && return _poly_to_num(t)
    return t::Complex{Num}
end

Base.show(io::IO, c::Coeff) = show(io, to_num(c))

# Branch on the tail type so each `isequal` / `hash` call sees a concrete operand
# (no `Union{Native,Poly,Complex{Num}}` split dispatch).
function Base.isequal(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Native
        return tb isa Native && isequal(a.z, b.z)
    elseif ta isa Poly
        return tb isa Poly && isequal(ta, tb)
    else
        return tb isa Complex{Num} && isequal(ta, tb)
    end
end
Base.:(==)(a::Coeff, b::Coeff) = isequal(a, b)
function Base.hash(c::Coeff, h::UInt)
    t = c.tail
    t isa Native && return hash(c.z, hash(:CoeffNative, h))
    t isa Poly && return hash(t, hash(:CoeffSym, h))
    return hash(t::Complex{Num}, hash(:CoeffSym, h))
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
    (_is_native(a) && _is_native(b)) && return _native(a.z / b.z)
    return _to_cnum(to_num(a) / to_num(b))
end
Base.:/(a::Coeff, b::Number) = a / _to_cnum(b)
Base.:/(a::Number, b::Coeff) = _to_cnum(a) / b

_sym_conj(x::Num) = SymbolicUtils.symtype(x) <: Real ? x : Num(conj(SymbolicUtils.unwrap(x)))

@inline function _conj_cnum(c::Coeff)
    _is_native(c) && return _native(conj(c.z))
    # conj(re + i*im) = conj(re) - i*conj(im); each part may be a Number-symtype
    # symbol. Materialized parts serve the polynomial and symbolic tails alike;
    # `_cnum` re-recognizes the polynomial after conjugation.
    return _cnum(_sym_conj(real(c)), -_sym_conj(imag(c)))
end

@inline function _iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, _NUM_ZERO)
end

@inline function _iszero_cnum(c::Coeff)
    _is_native(c) && return iszero(c.z)
    c.tail isa Poly && return false   # a canonical Poly never sums to zero
    return _iszero_num(real(c.tail)) && _iszero_num(imag(c.tail))
end

# `unwrap` returns a `BasicSymbolic` even for numeric constants, so test the
# node kind rather than `isa Number`.
@inline function _is_symbolic_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    return SymbolicUtils.issym(v) || SymbolicUtils.iscall(v)
end

@inline function _is_symbolic_cnum(c::Coeff)
    _is_native(c) && return false
    c.tail isa Poly && return true
    return _is_symbolic_num(real(c.tail)) || _is_symbolic_num(imag(c.tail))
end

# Structural `a == -b`; see `_addto_key!` for why this is needed.
@inline function _isneg_cnum(a::Coeff, b::Coeff)
    an = _is_native(a)
    bn = _is_native(b)
    an && bn && return isequal(a.z, -b.z)
    (an != bn) && return false   # one native, one symbolic: never exact negatives
    if a.tail isa Poly && b.tail isa Poly
        return isempty(_poly_add(a.tail.terms, b.tail.terms))
    end
    return isequal(real(a), -real(b)) && isequal(imag(a), -imag(b))
end

@inline function _mul_cnum(a::Coeff, b::Coeff)
    (_is_native(a) && _is_native(b)) && return _native(a.z * b.z)
    return _mul_cnum_slow(a, b)
end

# Polynomial fast paths (native scalar mul + integer-exponent merge, no CAS), then
# the materialized symbolic multiply. Most prefactors have zero imaginary part, so
# the symbolic path skips the extra Num mul/sub in that case.
@noinline function _mul_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return _from_poly(_poly_mul(ta.terms, tb.terms))
    elseif ta isa Poly && tb isa Native
        return _from_poly(_poly_scale(ta.terms, b.z))
    elseif tb isa Poly && ta isa Native
        return _from_poly(_poly_scale(tb.terms, a.z))
    end
    ca, cb = to_num(a), to_num(b)
    ar, ai = real(ca), imag(ca)
    br, bi = real(cb), imag(cb)
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return _cnum_sym(ar * br, _NUM_ZERO)
    elseif ai_zero
        return _cnum_sym(ar * br, ar * bi)
    elseif bi_zero
        return _cnum_sym(ar * br, ai * br)
    else
        return _cnum_sym(ar * br - ai * bi, ar * bi + ai * br)
    end
end

@inline function _neg_cnum(a::Coeff)
    t = a.tail
    t isa Native && return _native(-a.z)
    t isa Poly && return _from_poly(_poly_scale(t.terms, -_ONE_C))
    return _cnum_sym(-real(t), -imag(t))
end

@inline function _add_cnum(a::Coeff, b::Coeff)
    (_is_native(a) && _is_native(b)) && return _native(a.z + b.z)
    return _add_cnum_slow(a, b)
end

# Polynomial addition is a native merge (no escalation): this is what makes the
# tier pay off on sum-heavy workloads, where the single-monomial design regressed.
@noinline function _add_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return _from_poly(_poly_add(ta.terms, tb.terms))
    elseif ta isa Poly && tb isa Native
        return _from_poly(_poly_add(ta.terms, Monomial[Monomial(b.z, _EMPTY_SYMS, _EMPTY_EXPS)]))
    elseif tb isa Poly && ta isa Native
        return _from_poly(_poly_add(tb.terms, Monomial[Monomial(a.z, _EMPTY_SYMS, _EMPTY_EXPS)]))
    end
    ca, cb = to_num(a), to_num(b)
    ar, ai = real(ca), imag(ca)
    br, bi = real(cb), imag(cb)
    ai_zero = _iszero_num(ai)
    bi_zero = _iszero_num(bi)
    if ai_zero && bi_zero
        return _cnum_sym(ar + br, _NUM_ZERO)
    elseif ai_zero
        return _cnum_sym(ar + br, bi)
    elseif bi_zero
        return _cnum_sym(ar + br, ai)
    else
        return _cnum_sym(ar + br, ai + bi)
    end
end

_sort_key(op::QSym) = (op.space_index, op.index.name)
