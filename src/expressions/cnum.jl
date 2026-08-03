# Zero-size tag marking the native fast path (the value lives inline in `z`).
struct Native end
const NATIVE = Native()

"""
Coefficient representation for operator prefactors. A `Coeff` has three forms: a
native `ComplexF64` (concrete numbers), a `Poly` parameter polynomial
(products/sums of named parameters), and a `Complex{Num}` fallback. Native and
polynomial arithmetic stay off SymbolicUtils, lowering to `Complex{Num}` only at
the symbolic boundaries (`to_num`).
"""
struct Coeff
    z::ComplexF64
    tail::Union{Native, Poly, Complex{Num}}
end
const CNum = Coeff

const _NUM_ZERO = Num(0)
const _NUM_ONE = Num(1)
const _ONE_C = ComplexF64(1)
const _EMPTY_SYMS = SymbolicUtils.BasicSymbolic[]
const _EMPTY_EXPS = Rational{Int}[]

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
const _CNUM_HALF = _native(ComplexF64(0.5))

"""
    expim(x)

The unit phase `exp(im*x)`, as a coefficient a symbolic argument keeps intact.

The polynomial tier holds this as one indivisible factor, so phases multiply by adding
exponents and conjugate by negating them, and `expim(x)*expim(-x)` cancels to `1` on its
own. Splitting a phase into `cos`/`sin` instead leaves two unrelated factors whose
relationship has to be supplied separately, which is why rotating frames are built from
`expim` and why a time-dependent coefficient written this way composes with them.

Displays as `exp(im*x)`, and a conjugate pair folds back to `cos`/`sin` for display.

Returns a coefficient, not a `Num`: `conj` and multiplication by a complex scalar have to
act on the phase as one factor, and a `Num` is declared real, so both would silently split
it into unrelated `real`/`imag` halves.

```jldoctest
julia> using SecondQuantizedAlgebra: expim

julia> @variables ω t;

julia> h = FockSpace(:f); a = Destroy(h, :a);

julia> expim(ω * t) * expim(-ω * t) * a
a

julia> conj(expim(ω * t)) * expim(ω * t)
1
```
"""
expim(x::Number) = exp(im * x)
expim(x::Num) = _phase_coeff(x)
expim(x::SymbolicUtils.BasicSymbolic) = _phase_coeff(x)

@inline _is_phase(b) =
    b isa SymbolicUtils.BasicSymbolic &&
    SymbolicUtils.iscall(b) &&
    SymbolicUtils.operation(b) === expim

# The explicit scalar `shape` matters: without it the term is shaped `Unknown` and every
# later `Complex{Num}` addition it takes part in fails on a shape mismatch.
const _SCALAR_SHAPE = UnitRange{Int}[]
_expim_expanded(x) = SymbolicUtils.term(
    expim, SymbolicUtils.unwrap(x); type = Complex{Real}, shape = _SCALAR_SHAPE,
)
# `expand` first: `(ω + 2J)*t` and `ω*t + 2J*t` are one phase and must intern to one atom.
_expim(x) = _expim_expanded(expand(x))

# `type` and `shape` take part in hash-consing, so without these every `maketerm` rebuild
# (each `substitute`, each `Postwalk`) recomputes them from the generic fallbacks and mints
# a *different* atom, one whose symtype is `Real` and whose `conj` is therefore the identity.
SymbolicUtils.promote_symtype(::typeof(expim), ::SymbolicUtils.TypeT) = Complex{Real}
SymbolicUtils.promote_shape(::typeof(expim), ::SymbolicUtils.ShapeT) =
    SymbolicUtils.ShapeVecT()

# Without a rule `expand_derivatives` hits the global `nothing` fallback and leaves an inert
# `Differential` node. The body must yield a bare `BasicSymbolic`, hence `_expim`.
Symbolics.@register_derivative expim(x) 1 im * _expim(x)

function _leading_sign(v)
    u = SymbolicUtils.unwrap(v)
    if u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(u)
        op = SymbolicUtils.operation(u)
        args = SymbolicUtils.arguments(u)
        (op === (+) && !isempty(args)) && return _leading_sign(args[1])
        if op === (*)
            for f in args
                n = _const_value(SymbolicUtils.unwrap(f))
                n isa Real && !iszero(n) && return sign(n)
            end
            return 0.0
        end
    end
    n = _const_value(u)
    return n isa Real ? sign(n) : 0.0
end

# Never route this through `Symbolics.wrap`: `Num` cannot hold a complex symtype, so wrapping
# splits the phase into `real`/`imag` halves and loses the single factor everything below
# depends on. The argument is oriented so `expim(-x)` interns as `expim(x)` at exponent `-1`;
# without that, two modes rotating at opposite rates carry unrelated atoms and never cancel.
function _phase_coeff(x)
    a = expand(x)
    v = _const_value(SymbolicUtils.unwrap(a))
    # A phase over a literal is its value. Interning it instead would keep every coefficient
    # it touches off the native tier, so numerically cancelling terms would never fold.
    v isa Number && return _native(ComplexF64(exp(im * v)))
    neg = _leading_sign(a) < 0
    c = _atom_coeff(_expim_expanded(neg ? expand(-a) : a))
    return neg ? _conj_cnum(c) : c
end

# Exact value of an elementary function at argument `0`, or `nothing`. Folds the `exp(0)`
# Symbolics leaves after Euler-expanding `exp(im*ω*t)`. Spelled as `===` (not `in` a tuple)
# to stay statically resolved; user-registered functions (`pulse(t)`) return `nothing`.
@inline function _unary_at_zero(op)::Union{ComplexF64, Nothing}
    (op === exp || op === cos || op === cosh) && return _ONE_C
    (op === sin || op === tan || op === sinh || op === tanh) && return zero(ComplexF64)
    return nothing
end

# Concrete numeric content of an (unwrapped) symbolic value, else `nothing`
# (`BasicSymbolic` constants from substitute/simplify count, keeping them native).
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
_to_cnum(x::Num) = _recognize(SymbolicUtils.unwrap(x))
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
_to_cnum(x::SymbolicUtils.BasicSymbolic) = _recognize(x)

# Canonicalizing constructor from real/imag `Num` parts (`re + im*i`), used by the
# symbolic boundaries (substitute / conj / change_index) that may yield a polynomial.
_cnum(re::Num, im::Num) =
    _add_cnum(_recognize(SymbolicUtils.unwrap(re)), _mul_cnum(_CNUM_IM, _recognize(SymbolicUtils.unwrap(im))))

# Rebuild after a materialized symbolic arithmetic step: folds a numeric constant
# back to native, else stays symbolic. Must NOT re-enter `_recognize` (would recurse).
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

# Fold a canonical term list into a Coeff: empty -> zero, one constant term ->
# native, else a Poly. Inputs are already canonical, so no re-sort here.
function _from_poly(terms::Vector{Monomial})
    isempty(terms) && return _CNUM_ZERO
    if length(terms) == 1 && isempty(terms[1].syms)
        return _to_cnum(terms[1].scalar)
    end
    return _poly_coeff(Poly(terms))
end

# One monomial term -> Complex{Num}. The integer part of each exponent is built by
# repeated multiply/divide (both infer `Num`, unlike `Num ^ Int` which infers `Any`).
function _term_to_num(m::Monomial)
    prod = _NUM_ONE
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        e = m.exps[i]
        # `exp(-im*x)` rather than `1 / exp(im*x)`: same value, and it keeps an inverse
        # phase readable. Re-recognising it lands back on this atom, since `_recognize`
        # re-orients every `expim`.
        if e < 0 && _is_phase(s)
            s = _expim(-only(SymbolicUtils.arguments(s)))
            e = -e
        end
        base = Num(s)
        q = div(numerator(e), denominator(e))   # integer part, toward zero
        if q >= 0
            for _ in 1:q
                prod = prod * base
            end
        else
            for _ in 1:(-q)
                prod = prod / base
            end
        end
        f = e - q
        if f == 1 // 2
            prod = prod * sqrt(base)
        elseif f == -1 // 2
            prod = prod / sqrt(base)
        elseif f != 0
            prod = prod * (base^f)::Num   # `::Num`: `Num ^ Rational` infers `Any`
        end
    end
    # Guard the zero halves: `0 * x` does not always fold for a complex-symtype factor,
    # and an unfolded `0*expim(...)` would survive all the way to display.
    rz, iz = real(m.scalar), imag(m.scalar)
    re = iszero(rz) ? _NUM_ZERO : _num_from_float(rz) * prod
    imag_ = iszero(iz) ? _NUM_ZERO : _num_from_float(iz) * prod
    return Complex(re, imag_)
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

# An "atom" is an irreducible scalar the polynomial tier treats as one opaque
# variable: a symbol, an array index (`ω[i]`), or a non-algebraic one-arg call on an
# atom (`real(g)`, `imag(g)`, `sqrt`, `exp`, `conj`, ...). Algebraic ops (`+ * ^ /`,
# `complex`) are decomposed by `_recognize` instead, keeping their structure native.
@inline function _is_atom(b)
    b isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.issym(b) && return true
    SymbolicUtils.iscall(b) || return false
    op = SymbolicUtils.operation(b)
    op === getindex && return true
    # Atomic whatever its argument. The "argument must itself be an atom" rule below keeps
    # `cos(ω*t)` on the symbolic tail where the CAS can fold it; an `expim` has nothing for
    # the CAS to fold, and holding it off the polynomial tier breaks phase cancellation.
    op === expim && return true
    (op === (+) || op === (*) || op === (^) || op === (/) || op === complex) &&
        return false
    args = SymbolicUtils.arguments(b)
    return length(args) == 1 && _is_atom(only(args))
end

# A bare atom (symbol / array index / `conj(atom)`) as a single-monomial Coeff.
@inline _atom_coeff(x::SymbolicUtils.BasicSymbolic) =
    _poly_coeff(Poly(Monomial[Monomial(_ONE_C, SymbolicUtils.BasicSymbolic[x], Rational{Int}[1])]))
# An unrecognized symbolic value, kept on the `Complex{Num}` symbolic path.
@inline _sym_leaf(x::SymbolicUtils.BasicSymbolic) = _symbolic(Complex(Num(x), _NUM_ZERO))

# A fractional power `base^r`. Native only for a numeric base or a single-atom
# unit-scalar monomial (giving that atom a rational exponent); any other base would
# need to distribute the radical (unsound), so it becomes a symbolic leaf.
function _rational_power(basearg, r::Rational{Int}, x)
    base = _recognize(basearg)
    _is_native(base) && return _native(base.z^r)
    if base.tail isa Poly && length(base.tail.terms) == 1
        m = base.tail.terms[1]
        if length(m.syms) == 1 && m.scalar == _ONE_C
            return _poly_coeff(Poly(Monomial[Monomial(_ONE_C, m.syms, Rational{Int}[m.exps[1] * r])]))
        end
    end
    return _sym_leaf(x)
end

# Total recognizer: evaluate a symbolic expression tree in the Coeff algebra.
# Numbers/atoms map to native/Poly; `+ * ^(int) /` compose through the coefficient
# arithmetic; radicals fold to rational exponents; everything else is a symbolic
# leaf. Always returns a `Coeff` (no `nothing` sentinel).
_recognize(x::Number)::Coeff = _to_cnum(x)
_recognize(x::Num)::Coeff = _recognize(SymbolicUtils.unwrap(x))
_recognize(x)::Coeff = _to_cnum(x)
function _recognize(x::SymbolicUtils.BasicSymbolic)::Coeff
    SymbolicUtils.isconst(x) && return _recognize(x.val)
    SymbolicUtils.issym(x) && return _atom_coeff(x)
    SymbolicUtils.iscall(x) || return _sym_leaf(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        return _recognize_sum(args)
    elseif op === (*)
        return _recognize_prod(args)
    elseif op === (^)
        length(args) == 2 || return _sym_leaf(x)
        # `isa` guards narrow the `Any` exponent before arithmetic (a helper call on
        # the `Any` value would force a runtime dispatch).
        pv = _const_value(args[2])
        if pv isa Integer
            n = Int(pv)
            base = _recognize(args[1])
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
                    return _poly_coeff(Poly(Monomial[Monomial(s^(-n), m.syms, Rational{Int}[e * n for e in m.exps])]))
            end
        elseif pv isa Rational{Int}
            return _rational_power(args[1], pv, x)
        end
        return _sym_leaf(x)
    elseif op === getindex
        return _atom_coeff(x)
    elseif op === expim
        # Through `_phase_coeff`, not `_atom_coeff`: a phase recovered from a lowered
        # expression has to orient the same way as one built directly, or the two spellings
        # become unrelated atoms and stop cancelling.
        length(args) == 1 && return _phase_coeff(only(args))
        return _sym_leaf(x)
    elseif op === conj
        # Fold conj of a constant (e.g. conj(0.0) left by substituting a complex
        # parameter to a real value) so the coefficient can collapse to zero.
        inner = _recognize(args[1])
        _is_native(inner) && return _native(conj(inner.z))
        # Without this, `conj(expim(x))` interns as a fresh atom that never cancels against
        # the phase it came from.
        _is_phase(args[1]) && return _conj_cnum(inner)
        return _is_atom(x) ? _atom_coeff(x) : _sym_leaf(x)
    elseif op === (/)
        length(args) == 2 || return _sym_leaf(x)
        den = _recognize(args[2])
        if _is_native(den)
            return _mul_cnum(_recognize(args[1]), _native(inv(den.z)))
        elseif den.tail isa Poly && length(den.tail.terms) == 1
            m = den.tail.terms[1]
            s = inv(m.scalar)
            s * m.scalar == _ONE_C &&
                return _mul_cnum(
                _recognize(args[1]),
                _poly_coeff(Poly(Monomial[Monomial(s, m.syms, Rational{Int}[-e for e in m.exps])]))
            )
        end
        return _sym_leaf(x)
    elseif op === complex
        length(args) == 2 && return _cnum(Num(args[1]), Num(args[2]))
        return _sym_leaf(x)
    elseif op === sqrt
        length(args) == 1 && return _rational_power(only(args), 1 // 2, x)
        return _sym_leaf(x)
    elseif op === cbrt
        length(args) == 1 && return _rational_power(only(args), 1 // 3, x)
        return _sym_leaf(x)
    end
    # Fold an elementary function of a literal zero (`exp(0) -> 1`, `sin(0) -> 0`, ...).
    if length(args) == 1
        v = _unary_at_zero(op)
        if v !== nothing
            a1 = _recognize(only(args))
            (_is_native(a1) && iszero(a1.z)) && return _native(v)
        end
    end
    # Irreducible one-arg call on an atom (`exp`, `sin`, `real`, `imag`, ...): keep it
    # native as an opaque integer-exponent atom. Radicals are handled above instead.
    return _is_atom(x) ? _atom_coeff(x) : _sym_leaf(x)
end

function _recognize_sum(args)::Coeff
    monos = Monomial[]
    znative = zero(ComplexF64)
    sym = _CNUM_ZERO
    have_sym = false
    for a in args
        ca = _recognize(a)
        t = ca.tail
        if t isa Native
            znative += ca.z
        elseif t isa Poly
            append!(monos, t.terms)
        else
            sym = _add_cnum(sym, ca)
            have_sym = true
        end
    end
    znative == 0 || push!(monos, Monomial(znative, _EMPTY_SYMS, _EMPTY_EXPS))
    poly = _from_poly(_canonical_terms!(monos))
    return have_sym ? _add_cnum(poly, sym) : poly
end

# Insertion-sort a factor list (syms + matching exps) in place by objectid key.
# Avoids Base's `sortperm` machinery for the abstract `BasicSymbolic` eltype, whose
# inference is a nontrivial slice of first-call latency; these lists are tiny.
function _sort_factors!(syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}})
    @inbounds for i in 2:length(syms)
        s = syms[i]; e = exps[i]; k = _fkey(s)
        j = i - 1
        while j >= 1 && _fkey(syms[j]) > k
            syms[j + 1] = syms[j]; exps[j + 1] = exps[j]
            j -= 1
        end
        syms[j + 1] = s; exps[j + 1] = e
    end
    return nothing
end

function _merge_factor_list(syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}})
    n = length(syms)
    n <= 1 && return (syms, exps)
    _sort_factors!(syms, exps)
    osyms = SymbolicUtils.BasicSymbolic[]
    oexps = Rational{Int}[]
    sizehint!(osyms, n); sizehint!(oexps, n)
    i = 1
    @inbounds while i <= n
        s = syms[i]
        e = exps[i]
        j = i + 1
        while j <= n && syms[j] === s
            e += exps[j]
            j += 1
        end
        e != 0 && (push!(osyms, s); push!(oexps, e))
        i = j
    end
    return (osyms, oexps)
end

function _recognize_prod(args)::Coeff
    scalar = _ONE_C
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    other = _CNUM_ONE
    have_other = false
    for a in args
        ca = _recognize(a)
        t = ca.tail
        if t isa Native
            scalar *= ca.z
        elseif t isa Poly && length(t.terms) == 1
            m = t.terms[1]
            scalar *= m.scalar
            append!(syms, m.syms)
            append!(exps, m.exps)
        else
            other = _mul_cnum(other, ca)
            have_other = true
        end
    end
    iszero(scalar) && return _CNUM_ZERO
    ms, me = _merge_factor_list(syms, exps)
    mono = _from_poly(Monomial[Monomial(scalar + complex(0.0, 0.0), ms, me)])
    return have_other ? _mul_cnum(mono, other) : mono
end

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = _to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = _to_cnum(x)

Base.real(c::Coeff) = _is_native(c) ? _num_from_float(real(c.z)) : real(to_num(c))
Base.imag(c::Coeff) = _is_native(c) ? _num_from_float(imag(c.z)) : imag(to_num(c))

@inline function _realimag(c::Coeff)
    _is_native(c) && return (_num_from_float(real(c.z)), _num_from_float(imag(c.z)))
    cn = to_num(c)
    return (real(cn), imag(cn))
end

"""
    to_num(c::Coeff) -> Complex{Num}

Lower a stored [`Coeff`](@ref) to the public Symbolics representation used at
package boundaries. This is the supported way to read the coefficient returned
when iterating a [`QAdd`](@ref).
"""
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

# `conj(conj(x)) == x`, so unwrap an existing `conj(...)` rather than nesting a
# second one (which never folds and survives downstream).
_is_conj_call(x) =
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === conj
function _sym_conj(x::Num)
    SymbolicUtils.symtype(x) <: Real && return x
    u = SymbolicUtils.unwrap(x)
    _is_conj_call(u) && return Num(SymbolicUtils.arguments(u)[1])
    # Distribute over sums and products, so a real factor (a radical, a real-symtype
    # parameter) is left alone and a phase flips on its own. Wrapping the whole expression
    # leaves a `conj(...)` that never folds, and the residual built from it cannot cancel.
    if SymbolicUtils.iscall(u)
        op = SymbolicUtils.operation(u)
        args = SymbolicUtils.arguments(u)
        if op === (+)
            return sum(_sym_conj(Num(a)) for a in args)
        elseif op === (*)
            return prod(_sym_conj(Num(a)) for a in args)
        end
    end
    return Num(conj(u))
end

# Conjugate an atom factor: real-symtype atoms are self-conjugate, an existing
# `conj(...)` unwraps (involution), else wrap in `conj(...)` (still an atom to
# `_is_atom`).
@inline function _conj_atom(s::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.symtype(s) <: Real && return s
    _is_conj_call(s) && return SymbolicUtils.arguments(s)[1]
    return SymbolicUtils.unwrap(conj(s))
end

# Native conjugation of a Poly: conjugate each scalar and atom, re-sort the rekeyed
# factors by `objectid`, re-canonicalize. Avoids the `to_num` round-trip per term.
function _conj_poly(p::Poly)
    terms = Vector{Monomial}(undef, length(p.terms))
    @inbounds for k in eachindex(p.terms)
        m = p.terms[k]
        n = length(m.syms)
        if n == 0
            terms[k] = Monomial(conj(m.scalar), m.syms, m.exps)
            continue
        end
        nsyms = Vector{SymbolicUtils.BasicSymbolic}(undef, n)
        nexps = copy(m.exps)
        for i in 1:n
            s = m.syms[i]
            # Unimodular: conjugation flips the exponent and keeps the atom, which is what
            # lets `p * conj(p)` cancel instead of growing a second unrelated factor.
            if _is_phase(s)
                nsyms[i] = s
                nexps[i] = -nexps[i]
            else
                nsyms[i] = _conj_atom(s)
            end
        end
        _sort_factors!(nsyms, nexps)
        terms[k] = Monomial(conj(m.scalar), nsyms, nexps)
    end
    return _from_poly(_canonical_terms!(terms))
end

@inline function _conj_cnum(c::Coeff)
    t = c.tail
    t isa Native && return _native(conj(c.z))
    t isa Poly && return _conj_poly(t)
    # Symbolic tail: conj(re + i*im) = conj(re) - i*conj(im); each part may be a
    # Number-symtype symbol. `_cnum` re-recognizes the result.
    re, im = _realimag(c)
    return _cnum(_sym_conj(re), -_sym_conj(im))
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
    # `iszero(a.z + b.z)`, not `isequal(a.z, -b.z)`: negating a normalized `+0.0im`
    # gives `-0.0im`, and `isequal(0.0, -0.0)` is false, so the latter misses real
    # negatives (e.g. `1` vs `-1`).
    an && bn && return iszero(a.z + b.z)
    (an != bn) && return false   # one native, one symbolic: never exact negatives
    if a.tail isa Poly && b.tail isa Poly
        return isempty(_poly_add(a.tail.terms, b.tail.terms))
    end
    ar, ai = _realimag(a)
    br, bi = _realimag(b)
    return isequal(ar, -br) && isequal(ai, -bi)
end

# Materialize both operands once and return their real/imag parts plus the
# zero-imaginary flags the symbolic mul/add fast paths share.
@inline function _cnum_parts(a::Coeff, b::Coeff)
    ca, cb = to_num(a), to_num(b)
    ar, ai = real(ca), imag(ca)
    br, bi = real(cb), imag(cb)
    return (ar, ai, br, bi, _iszero_num(ai), _iszero_num(bi))
end

@inline function _mul_cnum(a::Coeff, b::Coeff)
    (_is_native(a) && _is_native(b)) && return _native(a.z * b.z)
    return _mul_cnum_slow(a, b)
end

# Native Poly fast paths first, then the materialized symbolic multiply (which skips
# the extra Num mul/sub when an operand has zero imaginary part).
@noinline function _mul_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return _from_poly(_poly_mul(ta.terms, tb.terms))
    elseif ta isa Poly && tb isa Native
        return _from_poly(_poly_scale(ta.terms, b.z))
    elseif tb isa Poly && ta isa Native
        return _from_poly(_poly_scale(tb.terms, a.z))
    end
    ar, ai, br, bi, ai_zero, bi_zero = _cnum_parts(a, b)
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
    # Skip add-by-zero: `_recognize` folds every sum from `_CNUM_ZERO`, so without this each
    # fold would splice a throwaway zero `Monomial` into the Poly and merge it away.
    _is_native(a) && iszero(a.z) && return b
    _is_native(b) && iszero(b.z) && return a
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
    ar, ai, br, bi, ai_zero, bi_zero = _cnum_parts(a, b)
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

_sort_key(op::QSym) = (op.space_index, _name_rank(op.index.name_id))
