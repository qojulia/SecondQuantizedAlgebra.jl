# Zero-size tag marking the native fast path (the value lives inline in `z`).
struct Native end
const NATIVE = Native()

struct RawSymbolicCoeff
    expr::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    real_slot::Bool

    function RawSymbolicCoeff(
            expr::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}, real_slot::Bool = false,
        )
        SymbolicUtils.symtype(expr) <: Number ||
            throw(ArgumentError("a symbolic coefficient must have numeric symtype"))
        isempty(SymbolicUtils.shape(expr)) ||
            throw(ArgumentError("a symbolic coefficient must be scalar"))
        return new(expr, real_slot)
    end
end

"""
Coefficient representation for operator prefactors. A `Coeff` has three forms: a
native `ComplexF64` (concrete numbers), a `Poly` parameter polynomial, and a raw
SymbolicUtils fallback. The latter preserves a single complex expression tree and
is lowered to `Complex{Num}` only at public boundaries (`to_num`).
"""
struct Coeff
    z::ComplexF64
    tail::Union{Native, Poly, RawSymbolicCoeff}
end
const CNum = Coeff

const _NUM_ZERO = Num(0)
const _NUM_ONE = Num(1)
const _ONE_C = ComplexF64(1)
const _ONE_R = ExactComplex(1 // 1, 0 // 1)
const _EMPTY_SYMS = SymbolicUtils.BasicSymbolic[]
const _EMPTY_EXPS = Rational{Int}[]

# Adding 0.0+0.0im normalizes any signed zero (`-0.0 -> 0.0`) so that structurally
# equal coefficients (e.g. `conj(2)` vs `2`) stay `isequal` and hash identically.
@inline _native(z::ComplexF64) = Coeff(z + complex(0.0, 0.0), NATIVE)
@inline _symbolic(
    x::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}; real_slot::Bool = false,
) = Coeff(zero(ComplexF64), RawSymbolicCoeff(x, real_slot))
@inline _poly_coeff(p::Poly) = Coeff(zero(ComplexF64), p)
@inline _is_native(c::Coeff) = c.tail isa Native
@inline _is_poly(c::Coeff) = c.tail isa Poly

const _CNUM_ZERO = _native(zero(ComplexF64))
const _CNUM_ONE = _native(one(ComplexF64))
const _CNUM_NEG1 = _native(-one(ComplexF64))
const _CNUM_IM = _native(ComplexF64(0, 1))
const _CNUM_NEG_IM = _native(ComplexF64(0, -1))
# This is an internal trigonometric constant, not user input. Keep it native so the
# Euler expansion remains on the hot path; exact user-provided rationals still use Poly.
const _CNUM_HALF = _native(ComplexF64(0.5))

"""
    expim(x)

Return the unit phase `exp(im*x)` for a provably real argument `x`.

Symbolic phases form a canonical multiplicative group: arguments add under multiplication,
integer powers scale the argument, and opposite phases cancel exactly. They remain compact
under conjugation, substitution, differentiation, and numerical evaluation. For example,
`expim(x) * expim(y) == expim(x + y)` and `expim(x) * expim(-x) == 1`. Phases always display
as exponentials; use
[`trigonometric_form`](@ref) for an explicit change of representation.

```jldoctest
julia> using SecondQuantizedAlgebra

julia> import SecondQuantizedAlgebra: expim

julia> @variables ω t;

julia> h = FockSpace(:f); a = Destroy(h, :a);

julia> expim(ω * t) * expim(-ω * t) * a
a

julia> conj(expim(ω * t)) * expim(ω * t)
1
```
"""
expim(x::Real) = exp(im * x)
expim(x::Num) = _phase_coeff(x)
expim(x::SymbolicUtils.BasicSymbolic) = _phase_coeff(x)

@noinline function _nonreal_phase_argument(x)
    throw(ArgumentError("`expim` requires a provably real argument; got `$x`"))
end

# A complex number is never accepted merely because its current imaginary part happens to
# be zero. The atom promises unit modulus structurally, so its domain has to be real by type,
# not by a value-dependent test.
expim(x::Number) = _nonreal_phase_argument(x)

@inline _is_phase(b) =
    b isa SymbolicUtils.BasicSymbolic &&
    SymbolicUtils.iscall(b) &&
    SymbolicUtils.operation(b) === expim

@inline _is_imaginary_unit(x) = isequal(x, Symbolics.IM)

@inline function _phase_factor_index(syms)::Int
    @inbounds for i in eachindex(syms)
        _is_phase(syms[i]) && return i
    end
    return 0
end

# The explicit scalar `shape` matters: without it the term is shaped `Unknown` and every
# later `Complex{Num}` addition it takes part in fails on a shape mismatch.
const _SCALAR_SHAPE = UnitRange{Int}[]
_raw_complex(re, im) = SymbolicUtils.term(
    complex, SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im);
    type = Complex{Real}, shape = _SCALAR_SHAPE,
)
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
function _negate_expanded(a::Num)::Num
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = _NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result -= Num(term)
        end
        return result
    end
    return -a
end

function _scale_expanded(a::Num, n::Int)::Num
    isone(n) && return a
    n == -1 && return _negate_expanded(a)
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = _NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result += n * Num(term)
        end
        return result
    end
    return n * a
end

function _phase_coeff_expanded(a::Num)
    u = SymbolicUtils.unwrap(a)
    v = _const_value(u)
    # A phase over a literal is its value. Interning it instead would keep every coefficient
    # it touches off the native tier, so numerically cancelling terms would never fold.
    v isa Real && return _native(ComplexF64(exp(im * v)))
    v isa Number && return _nonreal_phase_argument(a)
    (u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.symtype(u) <: Real) ||
        return _nonreal_phase_argument(a)
    neg = _leading_sign(a) < 0
    c = _atom_coeff(_expim_expanded(neg ? _negate_expanded(a) : a))
    return neg ? _conj_cnum(c) : c
end

_phase_coeff(x) = _phase_coeff_expanded(Num(expand(x)))

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

@inline function _phase_power(x)
    _is_phase(x) && return (only(SymbolicUtils.arguments(x)), 1)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    SymbolicUtils.operation(x) === (^) || return nothing
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return nothing
    _is_phase(args[1]) || return nothing
    exponent = _const_value(args[2])
    exponent isa Integer || return nothing
    return (only(SymbolicUtils.arguments(args[1])), Int(exponent))
end

# One bottom-up rewrite step for the identities SymbolicUtils cannot infer for the
# package-local `expim` operation. Multiplication is handled as one AC node so phase
# collection is independent of factor order and tree grouping.
function _rewrite_phase(x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if length(args) == 1 && _is_phase(args[1])
        argument = only(SymbolicUtils.arguments(args[1]))
        op === conj && return _expim_expanded(-argument)
        op === real && return cos(argument)
        op === imag && return sin(argument)
        (op === abs || op === abs2) && return 1
    elseif op === (/) && length(args) == 2
        angle = 0
        found_phase = false
        denominator = args[2]
        denominator_phase = _phase_power(denominator)
        if denominator_phase !== nothing
            argument, exponent = denominator_phase
            return args[1] * _expim_expanded(SymbolicUtils.expand(-exponent * argument))
        elseif SymbolicUtils.iscall(denominator) && SymbolicUtils.operation(denominator) === (*)
            factors = SymbolicUtils.arguments(denominator)
            ordinary_count = 0
            for factor in factors
                phase = _phase_power(factor)
                if phase === nothing
                    ordinary_count += 1
                else
                    argument, exponent = phase
                    angle += exponent * argument
                    found_phase = true
                end
            end
            found_phase || return nothing
            if ordinary_count == 0
                return args[1] * _expim_expanded(SymbolicUtils.expand(-angle))
            end
            ordinary_denominator = Any[]
            sizehint!(ordinary_denominator, ordinary_count)
            for factor in factors
                _phase_power(factor) === nothing && push!(ordinary_denominator, factor)
            end
            denominator = foldl(*, ordinary_denominator; init = 1)
            return (args[1] / denominator) * _expim_expanded(SymbolicUtils.expand(-angle))
        end
    elseif op === (*)
        ordinary = Any[]
        angle = 0
        phase_factors = 0
        changed = false
        for factor in args
            phase = _phase_power(factor)
            if phase === nothing
                push!(ordinary, factor)
                continue
            end
            argument, exponent = phase
            angle += exponent * argument
            phase_factors += 1
            changed |= exponent != 1
        end
        (phase_factors >= 2 || changed) || return nothing
        expanded_angle = SymbolicUtils.expand(angle)
        value = _const_value(expanded_angle)
        phase = value isa Number && iszero(value) ? 1 : _expim_expanded(expanded_angle)
        # Build the real amplitude before attaching the complex phase.  Multiplying the
        # phase by an exact real rational first makes SymbolicUtils promote `1//2` to a
        # floating-point complex constant.
        amplitude = foldl(*, ordinary; init = 1)
        phase === 1 && return amplitude
        # Keep the real amplitude in an explicit complex slot while rebuilding sums: a
        # real term and a complex phase term otherwise get promoted together as Float64.
        return _raw_complex(amplitude, 0 // 1) * phase
    elseif op === (^) && length(args) == 2 && _is_phase(args[1])
        exponent = _const_value(args[2])
        exponent isa Integer || return nothing
        argument = only(SymbolicUtils.arguments(args[1]))
        iszero(exponent) && return 1
        return _expim_expanded(SymbolicUtils.expand(exponent * argument))
    end
    return nothing
end

const _PHASE_NORMALIZER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Fixpoint(SymbolicUtils.Rewriters.Postwalk(_rewrite_phase)),
)

_normalize_phase(x) = _PHASE_NORMALIZER(x)
function _simplify_raw(x; kwargs...)
    normalized = _normalize_phase(x)
    simplified = SymbolicUtils.simplify(normalized; kwargs...)
    return _normalize_phase(simplified)
end

function _from_raw(x; normalize::Bool = true, real_slot::Bool = false)::Coeff
    value = _const_value(x)
    value isa Number && return _to_cnum(value)
    x isa SymbolicUtils.BasicSymbolic || return _to_cnum(x)
    expr = normalize ? _normalize_phase(x) : x
    value = _const_value(expr)
    value isa Number && return _to_cnum(value)
    # Complex-valued intermediate expressions can acquire an explicit zero imaginary
    # slot while expanding exact rational coefficients.  Canonicalize that case back to
    # the real symbolic tree so equivalent expressions do not differ only by `complex(..., 0)`.
    real_part, imaginary_part = _raw_realimag(expr)
    imaginary_value = _const_value(imaginary_part)
    if imaginary_value isa Number && iszero(imaginary_value)
        real_expr = real_part isa SymbolicUtils.BasicSymbolic ? real_part :
            SymbolicUtils.unwrap(Num(real_part))
        real_value = _const_value(real_expr)
        real_value isa Number && return _to_cnum(real_value)
        return _symbolic(real_expr; real_slot = true)
    end
    _is_phase(expr) && return _phase_coeff(only(SymbolicUtils.arguments(expr)))
    return _symbolic(expr; real_slot)
end

@inline function _from_raw_arithmetic(x, real_slot::Bool)::Coeff
    value = _const_value(x)
    value isa Number && return _to_cnum(value)
    x isa SymbolicUtils.BasicSymbolic || return _to_cnum(x)
    return _symbolic(x; real_slot)
end

@inline function _num_from_float(x::Float64)
    return (isinteger(x) && abs(x) <= 9.007199254740992e15) ? Num(Int(x)) : Num(x)
end

@inline _num_from_scalar(x::Rational{Int}) =
    denominator(x) == 1 ? Num(Int(numerator(x))) : Num(x)
@inline _num_from_scalar(x::Float64) = _num_from_float(x)

@inline function _exact_coeff(x::ExactComplex)::Coeff
    re, im = real(x), imag(x)
    if denominator(re) == 1 && denominator(im) == 1
        z = ComplexF64(Int(numerator(re)), Int(numerator(im)))
        z == x && return _native(z)
    end
    return _poly_coeff(Poly(Monomial[Monomial(x, _EMPTY_SYMS, _EMPTY_EXPS)]))
end

_to_cnum(x::Coeff) = x
_to_cnum(x::Num) = _recognize(SymbolicUtils.unwrap(x))
_to_cnum(x::Rational{Int}) = _exact_coeff(ExactComplex(x, 0 // 1))
_to_cnum(x::ExactComplex) = _exact_coeff(x)
# Native only when the value round-trips through ComplexF64 with no loss; non-rational
# values that cannot be represented faithfully (for example, large bignums) stay symbolic.
function _to_cnum(x::Real)
    z = ComplexF64(x)
    return z == x ? _native(z) : _symbolic(SymbolicUtils.unwrap(Num(x)))
end
function _to_cnum(x::Complex)
    z = ComplexF64(x)
    z == x && return _native(z)
    raw_re = SymbolicUtils.unwrap(Num(real(x)))
    raw_im = SymbolicUtils.unwrap(Num(imag(x)))
    return _symbolic(raw_re + im * raw_im)
end
function _to_cnum(x::Complex{Num})
    re, im = real(x), imag(x)
    raw_re, raw_im = SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im)
    im_value = _const_value(raw_im)
    if _is_phase(raw_re) && im_value isa Number && iszero(im_value)
        return _phase_coeff(only(SymbolicUtils.arguments(raw_re)))
    end
    # The fields of `Complex{Num}` are explicit real and imaginary coefficient slots,
    # even when a contained symbolic atom has the conservative `Number` symtype.
    # Canonicalize both through the coefficient algebra so materializing and then
    # rebuilding a coefficient does not switch between Raw and Poly tiers.
    return _cnum(re, im)
end
_to_cnum(x::SymbolicUtils.BasicSymbolic) = _recognize(x)

# Keep symbolic Fourier amplitudes exact without changing the native numeric fast path used
# by ordinary internal trigonometric calculations.
const _CNUM_HALF_EXACT = _exact_coeff(ExactComplex(1 // 2, 0 // 1))

# Canonicalizing constructor from real/imag `Num` parts (`re + im*i`), used by the
# symbolic boundaries (substitute / conj / change_index) that may yield a polynomial.
function _mul_by_im(c::Coeff)::Coeff
    tail = c.tail
    tail isa Native && return _native(im * c.z)
    tail isa Poly && return _from_poly(_poly_scale(tail.terms, ComplexF64(im)))
    return _from_raw(im * tail.expr)
end

_cnum(re::Num, im::Num) =
    _add_cnum(_recognize(SymbolicUtils.unwrap(re)), _mul_by_im(_recognize(SymbolicUtils.unwrap(im))))

# Rebuild after a materialized symbolic arithmetic step: folds a numeric constant
# back to native, else stays symbolic. Must NOT re-enter `_recognize` (would recurse).
function _cnum_sym(re::Num, im::Num)
    rv = _numeric_value(re)
    iv = _numeric_value(im)
    if rv isa Number && iv isa Number
        w = rv + Complex(false, true) * iv
        w isa ExactComplex && return _exact_coeff(w)
        z = ComplexF64(w)
        z == w && return _native(z)
    end
    raw_re = SymbolicUtils.unwrap(re)
    raw_im = SymbolicUtils.unwrap(im)
    return _from_raw(raw_re + Symbolics.IM * raw_im)
end

# === Parameter-polynomial tier: folding, materialization, recognition ===

# Polynomial monomial scalars retain exact Gaussian rationals; only an explicitly
# floating-point input or an operation involving one enters the ComplexF64 path.

# Fold a canonical term list into a Coeff: empty -> zero, one constant term ->
# native, else a Poly. Inputs are already canonical, so no re-sort here.
function _from_poly(terms::Vector{Monomial})
    isempty(terms) && return _CNUM_ZERO
    if length(terms) == 1 && isempty(terms[1].syms)
        scalar = terms[1].scalar
        return scalar isa ExactComplex ? _exact_coeff(scalar) : _to_cnum(scalar)
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
    real_scalar, imag_scalar = real(m.scalar), imag(m.scalar)
    re = iszero(real_scalar) ? _NUM_ZERO : _num_from_scalar(real_scalar) * prod
    imag_ = iszero(imag_scalar) ? _NUM_ZERO : _num_from_scalar(imag_scalar) * prod
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

function _term_to_raw(m::Monomial)
    result = if m.scalar isa ExactComplex
        re, im = real(m.scalar), imag(m.scalar)
        re_num = _num_from_scalar(re)
        im_num = _num_from_scalar(im)
        iszero(im) ? SymbolicUtils.unwrap(re_num) : _raw_complex(re_num, im_num)
    else
        m.scalar
    end
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        exponent = m.exps[i]
        factor = if denominator(exponent) == 1
            symbol^numerator(exponent)
        else
            symbol^exponent
        end
        # The polynomial convention for a non-real-symtype atom is carried by the
        # `real_slot` bit on the raw coefficient. Do not wrap the factor here: SymbolicUtils
        # simplifies `complex(z, 0)` back to `z` before the bit can be observed.
        result *= factor
    end
    return result
end

function _poly_to_raw(p::Poly)
    result = _term_to_raw(first(p.terms))
    @inbounds for i in 2:length(p.terms)
        result += _term_to_raw(p.terms[i])
    end
    return result
end

@inline function _raw_expression(c::Coeff)
    tail = c.tail
    tail isa Native && return c.z
    tail isa Poly && return _poly_to_raw(tail)
    return tail.expr
end

function _raw_realimag(x)
    x isa Number && return (real(x), imag(x))
    value = _const_value(x)
    value isa Number && return (real(value), imag(value))
    _is_imaginary_unit(x) && return (0, 1)
    SymbolicUtils.symtype(x) <: Real && return (x, 0)
    _is_phase(x) && begin
        argument = only(SymbolicUtils.arguments(x))
        return (cos(argument), sin(argument))
    end
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        op === complex && return (args[1], args[2])
        if op === (+)
            re, im = 0, 0
            for argument in args
                ar, ai = _raw_realimag(argument)
                re += ar
                im += ai
            end
            return (re, im)
        elseif op === (*)
            re, im = 1, 0
            for argument in args
                ar, ai = _raw_realimag(argument)
                re, im = re * ar - im * ai, re * ai + im * ar
            end
            return (re, im)
        elseif op === (/) && length(args) == 2 &&
                SymbolicUtils.symtype(args[2]) <: Real
            numerator = args[1]
            if SymbolicUtils.iscall(numerator) &&
                    SymbolicUtils.operation(numerator) === complex
                slots = SymbolicUtils.arguments(numerator)
                length(slots) == 2 || return (real(x), imag(x))
                return (slots[1] / args[2], slots[2] / args[2])
            end
        elseif op === conj
            re, im = _raw_realimag(only(args))
            return (re, -im)
        end
    end
    # An opaque Number-symtype expression may itself be complex-valued. Preserve the
    # public `real`/`imag` semantics here; the polynomial tier makes its real-slot
    # convention explicit in `_term_to_raw` before this fallback is reached.
    return (real(x), imag(x))
end

function _poly_is_real(p::Poly)
    for monomial in p.terms
        iszero(imag(monomial.scalar)) || return false
        for symbol in monomial.syms
            (_is_phase(symbol) || _is_imaginary_unit(symbol)) && return false
        end
    end
    return true
end

@inline _raw_is_real(t::RawSymbolicCoeff) =
    t.real_slot || SymbolicUtils.symtype(t.expr) <: Real

@inline function _cnum_is_real(c::Coeff)
    t = c.tail
    t isa Native && return iszero(imag(c.z))
    t isa Poly && return _poly_is_real(t)
    return _raw_is_real(t)
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

# Keep the identity scalar native for ordinary symbolic products; explicit rational
# inputs enter through `_exact_coeff` and remain exact when they are combined later.
# A bare atom (symbol / array index / `conj(atom)`) as a single-monomial Coeff.
@inline _atom_coeff(x::SymbolicUtils.BasicSymbolic) =
    _poly_coeff(Poly(Monomial[Monomial(_ONE_C, SymbolicUtils.BasicSymbolic[x], Rational{Int}[1])]))
# An unrecognized symbolic value, kept as one raw symbolic expression tree.
@inline _sym_leaf(x::SymbolicUtils.BasicSymbolic) = _from_raw(x; normalize = false)

# A fractional power `base^r`. Native only for a numeric base or a single-atom
# unit-scalar monomial (giving that atom a rational exponent); any other base would
# need to distribute the radical (unsound), so it becomes a symbolic leaf.
function _rational_power(basearg, r::Rational{Int}, x)
    base = _recognize(basearg)
    _is_native(base) && return _native(base.z^r)
    if base.tail isa Poly && length(base.tail.terms) == 1
        m = base.tail.terms[1]
        if length(m.syms) == 1 && (m.scalar == _ONE_C || m.scalar == _ONE_R)
            _is_phase(only(m.syms)) &&
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            return _poly_coeff(
                Poly(Monomial[Monomial(_ONE_C, m.syms, Rational{Int}[m.exps[1] * r])])
            )
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
    _is_imaginary_unit(x) && return _CNUM_IM
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
            return _pow_cnum_integer(_recognize(args[1]), Int(pv))
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
        numerator = _recognize(args[1])
        denominator = _const_value(args[2])
        if denominator isa Integer
            return _mul_cnum(numerator, _to_cnum(1 // denominator))
        end
        return numerator / _recognize(args[2])
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

@noinline function _canonical_phase_monomial(
        scalar::CoeffScalar,
        syms::Vector{SymbolicUtils.BasicSymbolic},
        exps::Vector{Rational{Int}},
    )
    ordinary_syms = SymbolicUtils.BasicSymbolic[]
    ordinary_exps = Rational{Int}[]
    sizehint!(ordinary_syms, length(syms) - 1)
    sizehint!(ordinary_exps, length(exps) - 1)
    angle = _NUM_ZERO
    @inbounds for i in eachindex(syms)
        symbol = syms[i]
        exponent = exps[i]
        if _is_phase(symbol)
            denominator(exponent) == 1 ||
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            argument = Num(only(SymbolicUtils.arguments(symbol)))
            angle += _scale_expanded(argument, numerator(exponent))
        else
            push!(ordinary_syms, symbol)
            push!(ordinary_exps, exponent)
        end
    end
    ordinary_syms, ordinary_exps = _merge_factor_list(ordinary_syms, ordinary_exps)
    phase = _phase_coeff_expanded(angle)
    if phase.tail isa Native
        return Monomial(
            _scalar_mul(scalar, phase.z), ordinary_syms, ordinary_exps,
        )
    end
    phase_poly = phase.tail::Poly
    phase_term = only(phase_poly.terms)
    append!(ordinary_syms, phase_term.syms)
    append!(ordinary_exps, phase_term.exps)
    _sort_factors!(ordinary_syms, ordinary_exps)
    return Monomial(
        _scalar_mul(scalar, phase_term.scalar), ordinary_syms, ordinary_exps,
    )
end

@noinline function _scaled_phase_monomial(
        scalar::CoeffScalar,
        symbol::Num,
        exponent::Rational{Int},
    )
    denominator(exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(symbol))))
    phase = _phase_coeff_expanded(_scale_expanded(argument, numerator(exponent)))
    if phase.tail isa Native
        return Monomial(
            _scalar_mul(scalar, phase.z), _EMPTY_SYMS, _EMPTY_EXPS,
        )
    end
    phase_term = only((phase.tail::Poly).terms)
    return Monomial(
        _scalar_mul(scalar, phase_term.scalar),
        phase_term.syms,
        phase_term.exps,
    )
end

@noinline function _merged_phase_monomial(
        scalar::CoeffScalar,
        left::Num,
        left_exponent::Rational{Int},
        right::Num,
        right_exponent::Rational{Int},
    )
    denominator(left_exponent) == denominator(right_exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    left_argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(left))))
    right_argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(right))))
    left_power = numerator(left_exponent)
    right_power = numerator(right_exponent)
    common_negative = left_power < 0 && right_power < 0
    if common_negative
        left_power = -left_power
        right_power = -right_power
    end
    angle = _scale_expanded(left_argument, left_power) +
        _scale_expanded(right_argument, right_power)
    phase = _phase_coeff_expanded(angle)
    if phase.tail isa Native
        return Monomial(
            _scalar_mul(scalar, phase.z), _EMPTY_SYMS, _EMPTY_EXPS,
        )
    end
    phase_term = only((phase.tail::Poly).terms)
    exponents = common_negative ? -phase_term.exps : phase_term.exps
    return Monomial(
        _scalar_mul(scalar, phase_term.scalar),
        phase_term.syms,
        exponents,
    )
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
            scalar = _scalar_mul(scalar, ca.z)
        elseif t isa Poly && length(t.terms) == 1
            m = t.terms[1]
            scalar = _scalar_mul(scalar, m.scalar)
            append!(syms, m.syms)
            append!(exps, m.exps)
        else
            other = _mul_cnum_slow(other, ca)
            have_other = true
        end
    end
    iszero(scalar) && return _CNUM_ZERO
    phase_count = count(_is_phase, syms)
    monomial = if phase_count <= 1
        ms, me = _merge_factor_list(syms, exps)
        Monomial(_normalize_scalar(scalar), ms, me)
    else
        _canonical_phase_monomial(scalar, syms, exps)
    end
    mono = _from_poly(Monomial[monomial])
    return have_other ? _mul_cnum_slow(mono, other) : mono
end

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = _to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = _to_cnum(x)

function _pure_phase_data(c::Coeff)::Union{Nothing, Tuple{CoeffScalar, Num, Float64}}
    tail = c.tail
    tail isa Poly || return nothing
    length(tail.terms) == 1 || return nothing
    monomial = only(tail.terms)
    length(monomial.syms) == 1 || return nothing
    symbol = only(monomial.syms)
    _is_phase(symbol) || return nothing
    exponent = only(monomial.exps)
    denominator(exponent) == 1 || return nothing
    abs2(monomial.scalar) == 1.0 || return nothing
    argument = Num(only(SymbolicUtils.arguments(symbol)))
    angle = _scale_expanded(argument, numerator(exponent))
    sine_sign = 1.0
    if _leading_sign(angle) < 0
        angle = _negate_expanded(angle)
        sine_sign = -1.0
    end
    return (monomial.scalar, angle, sine_sign)
end

@inline function _scaled_trig(a::Real, trig, angle::Num)::Num
    iszero(a) && return _NUM_ZERO
    isone(a) && return trig(angle)
    a == -1 && return -trig(angle)
    return _num_from_scalar(a) * trig(angle)
end

function Base.real(c::Coeff)::Num
    _is_native(c) && return _num_from_float(real(c.z))
    phase = _pure_phase_data(c)
    if phase === nothing
        c.tail isa Poly && return real(_poly_to_num(c.tail))
        expr = _normalize_phase(_raw_expression(c))
        re, _ = _cnum_is_real(c) ? (expr, 0) : _raw_realimag(expr)
        return Num(_normalize_phase(re))
    end
    scalar, angle, sine_sign = phase
    return _scaled_trig(real(scalar), cos, angle) -
        _scaled_trig(imag(scalar) * sine_sign, sin, angle)
end

function Base.imag(c::Coeff)::Num
    _is_native(c) && return _num_from_float(imag(c.z))
    phase = _pure_phase_data(c)
    if phase === nothing
        c.tail isa Poly && return imag(_poly_to_num(c.tail))
        expr = _normalize_phase(_raw_expression(c))
        _, im = _cnum_is_real(c) ? (0, 0) : _raw_realimag(expr)
        return Num(_normalize_phase(im))
    end
    scalar, angle, sine_sign = phase
    return _scaled_trig(real(scalar) * sine_sign, sin, angle) +
        _scaled_trig(imag(scalar), cos, angle)
end

function Base.abs(c::Coeff)::Num
    _is_native(c) && return _num_from_float(abs(c.z))
    _pure_phase_data(c) === nothing && throw(MethodError(abs, (c,)))
    return _NUM_ONE
end

function Base.abs2(c::Coeff)::Num
    _is_native(c) && return _num_from_float(abs2(c.z))
    _pure_phase_data(c) === nothing && throw(MethodError(abs2, (c,)))
    return _NUM_ONE
end

@inline function _realimag(c::Coeff)
    _is_native(c) && return (_num_from_float(real(c.z)), _num_from_float(imag(c.z)))
    cn = to_num(c)
    return (real(cn), imag(cn))
end

@inline function _raw_parts(c::Coeff)
    tail = c.tail::RawSymbolicCoeff
    _cnum_is_real(c) && return (Num(tail.expr), _NUM_ZERO)
    return _realimag(c)
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
    expr = _normalize_phase(t.expr)
    re, im = _cnum_is_real(c) ? (expr, 0) : _raw_realimag(expr)
    return Complex(Num(_normalize_phase(re)), Num(_normalize_phase(im)))
end

Base.show(io::IO, c::Coeff) = show(io, to_num(c))

# Branch on the tail type so each `isequal` / `hash` call sees a concrete operand.
function Base.isequal(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Native
        return tb isa Native && isequal(a.z, b.z)
    elseif ta isa Poly
        return tb isa Poly && isequal(ta, tb)
    else
        return tb isa RawSymbolicCoeff && ta.real_slot == tb.real_slot &&
            isequal(ta.expr, tb.expr)
    end
end
Base.:(==)(a::Coeff, b::Coeff) = isequal(a, b)
function Base.hash(c::Coeff, h::UInt)
    t = c.tail
    t isa Native && return hash(c.z, hash(:CoeffNative, h))
    t isa Poly && return hash(t, hash(:CoeffSym, h))
    return hash(t.real_slot, hash(t.expr, hash(:CoeffRaw, h)))
end

# Coefficients are routinely compared against plain numbers / `Complex{Num}`
# (e.g. `prefactor(x) == 2`, `q[key] == _CNUM_ONE`); promote the number side.
Base.isequal(a::Coeff, b::Number) = isequal(a, _to_cnum(b))
Base.isequal(a::Number, b::Coeff) = isequal(_to_cnum(a), b)
Base.:(==)(a::Coeff, b::Number) = isequal(a, _to_cnum(b))
Base.:(==)(a::Number, b::Coeff) = isequal(_to_cnum(a), b)

Base.iszero(c::Coeff) = _iszero_cnum(c)
Base.conj(c::Coeff) = _conj_cnum(c)

function Base.inv(c::Coeff)::Coeff
    tail = c.tail
    tail isa Native && return _native(inv(c.z))
    if tail isa Poly && length(tail.terms) == 1
        monomial = only(tail.terms)
        return _poly_coeff(
            Poly(
                Monomial[
                    Monomial(
                        _scalar_inv(monomial.scalar),
                        monomial.syms,
                        -monomial.exps,
                    ),
                ]
            )
        )
    end
    return _from_raw(inv(tail.expr); real_slot = tail.real_slot)
end

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
function Base.:/(a::Coeff, b::Coeff)::Coeff
    (_is_native(a) && _is_native(b)) && return _native(a.z / b.z)
    if b.tail isa Native && a.tail isa Poly
        return _from_poly(_poly_scale(a.tail.terms, inv(b.z)))
    end
    if b.tail isa Poly && length(b.tail.terms) == 1
        inverse = inv(b)
        inverse_tail = inverse.tail::Poly
        if a.tail isa Native
            return _from_poly(_poly_scale(inverse_tail.terms, a.z))
        elseif a.tail isa Poly
            return _from_poly(_poly_mul(a.tail.terms, inverse_tail.terms))
        end
    end
    return _from_raw(
        _raw_expression(a) / _raw_expression(b);
        real_slot = _cnum_is_real(a) && _cnum_is_real(b),
    )
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

function _raw_conj(x)
    x isa Number && return conj(x)
    SymbolicUtils.symtype(x) <: Real && return x
    _is_phase(x) && return _expim_expanded(-only(SymbolicUtils.arguments(x)))
    _is_conj_call(x) && return only(SymbolicUtils.arguments(x))
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === complex
            return _raw_complex(_raw_conj(args[1]), -_raw_conj(args[2]))
        end
        op === (+) && return foldl(+, map(_raw_conj, args))
        op === (*) && return foldl(*, map(_raw_conj, args))
        if op === (/)
            return _raw_conj(args[1]) / _raw_conj(args[2])
        end
    end
    return conj(x)
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
    return _from_raw(_raw_conj(t.expr); real_slot = t.real_slot)
end

@inline function _iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, _NUM_ZERO)
end

@inline function _iszero_cnum(c::Coeff)
    _is_native(c) && return iszero(c.z)
    c.tail isa Poly && return false   # a canonical Poly never sums to zero
    value = _const_value(c.tail.expr)
    return value isa Number && iszero(value)
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
    return true
end

# Structural `a == -b`, used to recognize exact cancellation without a CAS round-trip.
@inline function _raw_is_negative_of(
        positive::SymbolicUtils.BasicSymbolic, negative::SymbolicUtils.BasicSymbolic,
    )
    positive_value = _const_value(positive)
    negative_value = _const_value(negative)
    positive_value isa Number && negative_value isa Number &&
        positive_value == -negative_value && return true
    SymbolicUtils.iscall(negative) || return false
    op = SymbolicUtils.operation(negative)
    args = SymbolicUtils.arguments(negative)
    if op === (*) && length(args) >= 2 && _const_value(args[1]) == -1
        if SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (*)
            positive_args = SymbolicUtils.arguments(positive)
            length(args) - 1 == length(positive_args) || return false
            @inbounds for i in eachindex(positive_args)
                isequal(args[i + 1], positive_args[i]) || return false
            end
            return true
        end
        return length(args) == 2 && isequal(args[2], positive)
    end
    if op === (+) && SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (+)
        positive_args = SymbolicUtils.arguments(positive)
        length(args) == length(positive_args) || return false
        @inbounds for i in eachindex(positive_args)
            _raw_is_negative_of(positive_args[i], args[i]) || return false
        end
        return true
    end
    if op === (/) && length(args) == 2 &&
            SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (/)
        positive_args = SymbolicUtils.arguments(positive)
        length(positive_args) == 2 || return false
        return isequal(args[2], positive_args[2]) &&
            (
            _raw_is_negative_of(positive_args[1], args[1]) ||
                _raw_is_negative_of(args[1], positive_args[1])
        )
    end
    return false
end

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
    if a.tail isa RawSymbolicCoeff && b.tail isa RawSymbolicCoeff
        return a.tail.real_slot == b.tail.real_slot &&
            (
            _raw_is_negative_of(a.tail.expr, b.tail.expr) ||
                _raw_is_negative_of(b.tail.expr, a.tail.expr)
        )
    end
    ar, ai = _realimag(a)
    br, bi = _realimag(b)
    return isequal(ar, -br) && isequal(ai, -bi)
end

@inline function _mul_cnum(a::Coeff, b::Coeff)
    (_is_native(a) && _is_native(b)) && return _native(a.z * b.z)
    return _mul_cnum_slow(a, b)
end

# Native and polynomial fast paths first, then multiply the intact raw expressions.
@noinline function _mul_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return _from_poly(_poly_mul(ta.terms, tb.terms))
    elseif ta isa Poly && tb isa Native
        return _from_poly(_poly_scale(ta.terms, b.z))
    elseif tb isa Poly && ta isa Native
        return _from_poly(_poly_scale(tb.terms, a.z))
    end
    return _from_raw_arithmetic(
        _raw_expression(a) * _raw_expression(b), _cnum_is_real(a) && _cnum_is_real(b),
    )
end

@inline function _neg_cnum(a::Coeff)
    t = a.tail
    t isa Native && return _native(-a.z)
    t isa Poly && return _from_poly(_poly_scale(t.terms, -_ONE_C))
    return _from_raw_arithmetic(-t.expr, t.real_slot)
end

@inline function _add_cnum(a::Coeff, b::Coeff)
    (_is_native(a) && _is_native(b)) && return _native(a.z + b.z)
    # Skip add-by-zero: `_recognize` folds every sum from `_CNUM_ZERO`, so without this each
    # fold would splice a throwaway zero `Monomial` into the Poly and merge it away.
    _is_native(a) && iszero(a.z) && return b
    _is_native(b) && iszero(b.z) && return a
    # Raw symbolic addition deliberately avoids a CAS round-trip. Recover the common exact
    # cancellation case here, which is needed by scalar identities such as the off-diagonal
    # entries of a symbolic orthogonal matrix product.
    (a.tail isa RawSymbolicCoeff && b.tail isa RawSymbolicCoeff && _isneg_cnum(a, b)) &&
        return _CNUM_ZERO
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
    return _from_raw_arithmetic(
        _raw_expression(a) + _raw_expression(b), _cnum_is_real(a) && _cnum_is_real(b),
    )
end

@inline function _pow_cnum_nonnegative(base::Coeff, n::Int)
    result = _CNUM_ONE
    while n > 0
        if isodd(n)
            result = if _is_native(result) && _is_native(base)
                _native(result.z * base.z)
            else
                _mul_cnum_slow(result, base)
            end
        end
        n >>= 1
        if n > 0
            base = _is_native(base) ? _native(base.z * base.z) : _mul_cnum_slow(base, base)
        end
    end
    return result
end

function _pow_cnum_integer(base::Coeff, n::Int)::Coeff
    n >= 0 && return _pow_cnum_nonnegative(base, n)
    n == typemin(Int) && return _to_cnum(to_num(base)^n)
    return _CNUM_ONE / _pow_cnum_nonnegative(base, -n)
end

Base.:^(base::Coeff, n::Integer)::Coeff = _pow_cnum_integer(base, Int(n))

@inline function _rewrite_complex_slots(x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    SymbolicUtils.operation(x) === complex || return nothing
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return nothing
    return args[1] + im * args[2]
end

const _COMPLEX_SLOT_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(_rewrite_complex_slots),
)

@inline _lower_complex_slots(x) = _COMPLEX_SLOT_REWRITER(x)
@inline _derivative_expression(c::Coeff) = _lower_complex_slots(_raw_expression(c))

(D::Symbolics.Differential)(c::Coeff) =
    _is_native(c) ? D(to_num(c)) : D(_derivative_expression(c))

function Symbolics.derivative(c::Coeff, var; simplify = false, kwargs...)::Coeff
    _is_native(c) && return _CNUM_ZERO
    input = _derivative_expression(c)
    differentiated = Symbolics.expand_derivatives(
        Symbolics.Differential(var)(input), simplify; kwargs...,
    )
    return _to_cnum(differentiated)
end

@inline _factor_cnum(s::SymbolicUtils.BasicSymbolic, e::Rational{Int}) =
    _poly_coeff(Poly(Monomial[Monomial(_ONE_C, SymbolicUtils.BasicSymbolic[s], Rational{Int}[e])]))

@inline function _euler_cos(argument)
    phase = _phase_coeff(argument)
    return _mul_cnum(_CNUM_HALF_EXACT, _add_cnum(phase, _conj_cnum(phase)))
end

@inline function _euler_sin(argument)
    phase = _phase_coeff(argument)
    difference = _add_cnum(phase, _neg_cnum(_conj_cnum(phase)))
    return _mul_cnum(_mul_cnum(_CNUM_NEG_IM, _CNUM_HALF_EXACT), difference)
end

function _phase_from_exponential_argument(argument)::Union{Coeff, Nothing}
    exponent = _recognize(argument)
    real_part, imaginary_part = _realimag(exponent)
    _iszero_num(real_part) || return nothing
    imaginary = SymbolicUtils.unwrap(imaginary_part)
    SymbolicUtils.symtype(imaginary) <: Real || return nothing
    return _phase_coeff(imaginary_part)
end

function _exponential_tree(x)::Coeff
    u = SymbolicUtils.unwrap(x)
    u isa Number && return _to_cnum(u)
    u isa SymbolicUtils.BasicSymbolic || return _to_cnum(u)
    SymbolicUtils.iscall(u) || return _recognize(u)
    op = SymbolicUtils.operation(u)
    args = SymbolicUtils.arguments(u)
    if op === exp && length(args) == 1
        phase = _phase_from_exponential_argument(only(args))
        phase === nothing || return phase
        return _recognize(u)
    elseif op === cis && length(args) == 1
        return _phase_coeff(only(args))
    elseif op === cos && length(args) == 1
        return _euler_cos(only(args))
    elseif op === sin && length(args) == 1
        return _euler_sin(only(args))
    elseif op === (+)
        result = _CNUM_ZERO
        for arg in args
            result = _add_cnum(result, _exponential_tree(arg))
        end
        return result
    elseif op === (*)
        result = _CNUM_ONE
        for arg in args
            result = _mul_cnum(result, _exponential_tree(arg))
        end
        return result
    elseif op === (/) && length(args) == 2
        return _exponential_tree(args[1]) / _exponential_tree(args[2])
    elseif op === (^) && length(args) == 2
        exponent = _const_value(args[2])
        exponent isa Integer || return _recognize(u)
        return _pow_cnum_integer(_exponential_tree(args[1]), Int(exponent))
    end
    return _recognize(u)
end

function _exponential_monomial(m::Monomial)
    result = _to_cnum(m.scalar)
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        exponent = m.exps[i]
        factor = _factor_cnum(symbol, exponent)
        if denominator(exponent) == 1 && SymbolicUtils.iscall(symbol)
            op = SymbolicUtils.operation(symbol)
            if op === cos || op === sin
                n = numerator(exponent)
                base = op === cos ?
                    _euler_cos(only(SymbolicUtils.arguments(symbol))) :
                    _euler_sin(only(SymbolicUtils.arguments(symbol)))
                factor = _pow_cnum_integer(base, n)
            end
        end
        result = _mul_cnum(result, factor)
    end
    return result
end

function _exponential_cnum(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    if tail isa Poly
        result = _native(c.z)
        for monomial in tail.terms
            result = _add_cnum(result, _exponential_monomial(monomial))
        end
        return result
    end
    return _exponential_tree(tail.expr)
end

@inline function _phase_trigonometric(symbol::SymbolicUtils.BasicSymbolic, n::Int)
    argument = only(SymbolicUtils.arguments(symbol))
    angle = expand(Num(n) * Num(argument))
    sine_sign = _CNUM_ONE
    if _leading_sign(angle) < 0
        angle = expand(-angle)
        sine_sign = _CNUM_NEG1
    end
    sine = _mul_cnum(sine_sign, _to_cnum(sin(angle)))
    return _add_cnum(_to_cnum(cos(angle)), _mul_cnum(_CNUM_IM, sine))
end

function _trigonometric_monomial(m::Monomial)
    result = _to_cnum(m.scalar)
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        exponent = m.exps[i]
        factor = if _is_phase(symbol) && denominator(exponent) == 1
            _phase_trigonometric(symbol, numerator(exponent))
        else
            _factor_cnum(symbol, exponent)
        end
        result = _mul_cnum(result, factor)
    end
    return result
end

@inline function _rewrite_phase_to_trig(x)
    _is_phase(x) || return nothing
    argument = only(SymbolicUtils.arguments(x))
    return cos(argument) + im * sin(argument)
end

const _TRIGONOMETRIC_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(_rewrite_phase_to_trig),
)

function _trigonometric_cnum(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    if tail isa RawSymbolicCoeff
        expanded = SymbolicUtils.expand(_TRIGONOMETRIC_REWRITER(tail.expr))
        return _from_raw(
            _simplify_raw(expanded); normalize = false,
        )
    end
    result = _native(c.z)
    for monomial in tail.terms
        result = _add_cnum(result, _trigonometric_monomial(monomial))
    end
    return result.tail isa RawSymbolicCoeff ? _from_raw(
            _simplify_raw(SymbolicUtils.expand(result.tail.expr)); normalize = false,
        ) : result
end

"""
    exponential_form(x)

Rewrite algebraic occurrences of `cos(θ)` and `sin(θ)` in a coefficient or quantum
expression using the exact phase atom [`expim`](@ref). The conversion is explicit and does
not affect ordinary display or [`simplify`](@ref).

See also [`trigonometric_form`](@ref).
"""
exponential_form(x::Coeff) = _exponential_cnum(x)
exponential_form(x::Num) = _exponential_tree(x)
exponential_form(x::SymbolicUtils.BasicSymbolic) = _exponential_tree(x)
exponential_form(x::Number) = x
exponential_form(x::Complex{Num}) = _exponential_cnum(_to_cnum(x))

"""
    PhaseTerm

One term in the finite phase decomposition of a [`Coeff`](@ref). It represents
`amplitude * expim(phase)`.
"""
struct PhaseTerm
    amplitude::Coeff
    phase::Num
end

@inline _constant_phase_term(c::Coeff) = PhaseTerm(c, _NUM_ZERO)

function _monomial_phase_term(m::Monomial)
    phase_index = _phase_factor_index(m.syms)
    for i in eachindex(m.syms)
        i == phase_index && continue
        _contains_phase(m.syms[i]) && throw(
            ArgumentError("a phase occurs inside an operation with no finite phase decomposition"),
        )
    end
    phase_index == 0 && return _constant_phase_term(_from_poly(Monomial[m]))

    exponent = m.exps[phase_index]
    denominator(exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    argument = only(SymbolicUtils.arguments(m.syms[phase_index]))
    phase = expand(Num(numerator(exponent)) * Num(argument))

    syms = copy(m.syms)
    exps = copy(m.exps)
    deleteat!(syms, phase_index)
    deleteat!(exps, phase_index)
    amplitude = _from_poly(Monomial[Monomial(m.scalar, syms, exps)])
    return PhaseTerm(amplitude, phase)
end

function _contains_phase(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    _is_phase(x) && return true
    SymbolicUtils.iscall(x) || return false
    return any(_contains_phase, SymbolicUtils.arguments(x))
end

function _multiply_phase_terms(left::Vector{PhaseTerm}, right::Vector{PhaseTerm})
    result = PhaseTerm[]
    sizehint!(result, length(left) * length(right))
    for a in left, b in right
        amplitude = _mul_cnum(a.amplitude, b.amplitude)
        iszero(amplitude) && continue
        push!(result, PhaseTerm(amplitude, expand(a.phase + b.phase)))
    end
    return result
end

function _phase_power_terms(base, exponent::Int)
    exponent < 0 && _contains_phase(base) && throw(
        ArgumentError("a phase-bearing denominator does not have a finite phase decomposition"),
    )
    exponent < 0 && return PhaseTerm[_constant_phase_term(_pow_cnum_integer(_recognize(base), exponent))]

    result = PhaseTerm[_constant_phase_term(_CNUM_ONE)]
    factor = _raw_phase_terms(base)
    for _ in 1:exponent
        result = _multiply_phase_terms(result, factor)
    end
    return result
end

function _raw_phase_terms(x)::Vector{PhaseTerm}
    value = _const_value(x)
    value isa Number && return PhaseTerm[_constant_phase_term(_to_cnum(value))]
    x isa SymbolicUtils.BasicSymbolic ||
        return PhaseTerm[_constant_phase_term(_to_cnum(x))]

    phase = _phase_power(x)
    if phase !== nothing
        argument, exponent = phase
        return PhaseTerm[PhaseTerm(_CNUM_ONE, expand(Num(exponent) * Num(argument)))]
    end

    SymbolicUtils.iscall(x) || return PhaseTerm[_constant_phase_term(_recognize(x))]
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        result = PhaseTerm[]
        for argument in args
            append!(result, _raw_phase_terms(argument))
        end
        return result
    elseif op === (*)
        result = PhaseTerm[_constant_phase_term(_CNUM_ONE)]
        for factor in args
            result = _multiply_phase_terms(result, _raw_phase_terms(factor))
        end
        return result
    elseif op === (/) && length(args) == 2
        numerator, denominator_ = args
        _contains_phase(denominator_) && throw(
            ArgumentError(
                "a phase-bearing denominator does not have a finite phase decomposition",
            ),
        )
        denominator_coeff = _recognize(denominator_)
        return PhaseTerm[
            PhaseTerm(term.amplitude / denominator_coeff, term.phase) for
                term in _raw_phase_terms(numerator)
        ]
    elseif op === (^) && length(args) == 2
        exponent = _const_value(args[2])
        exponent isa Integer || begin
            _contains_phase(args[1]) && throw(
                ArgumentError("a unit phase cannot have a non-integer power"),
            )
            return PhaseTerm[_constant_phase_term(_recognize(x))]
        end
        return _phase_power_terms(args[1], Int(exponent))
    elseif op === conj && length(args) == 1
        return PhaseTerm[
            PhaseTerm(conj(term.amplitude), -term.phase) for
                term in _raw_phase_terms(only(args))
        ]
    end

    _contains_phase(x) && throw(
        ArgumentError("a phase occurs inside an operation with no finite phase decomposition"),
    )
    return PhaseTerm[_constant_phase_term(_recognize(x))]
end

"""
    phase_terms(c::Coeff) -> Vector{PhaseTerm}

Decompose a coefficient into a finite sum of `amplitude * expim(phase)` terms. Trigonometric
factors are converted to exponential form first. A phase-bearing denominator or a phase
inside an unsupported nonlinear operation throws an `ArgumentError` because it does not
represent a finite phase polynomial.

The decomposition is exact:

```julia
c == sum(term.amplitude * expim(term.phase) for term in phase_terms(c))
```
"""
function phase_terms(c::Coeff)::Vector{PhaseTerm}
    normalized = _exponential_cnum(c)
    tail = normalized.tail
    if tail isa Native
        return iszero(normalized) ? PhaseTerm[] : PhaseTerm[_constant_phase_term(normalized)]
    elseif tail isa Poly
        return PhaseTerm[_monomial_phase_term(monomial) for monomial in tail.terms]
    end
    return _raw_phase_terms(tail.expr)
end

"""
    trigonometric_form(x)

Rewrite integer powers of [`expim(θ)`](@ref expim) in a coefficient or quantum expression
as `cos(nθ) + im*sin(nθ)`. Opposite phases then combine through ordinary coefficient
arithmetic. The conversion is explicit and leaves the stored input unchanged.

See also [`exponential_form`](@ref).
"""
trigonometric_form(x::Coeff) = _trigonometric_cnum(x)
trigonometric_form(x::Num) = _trigonometric_cnum(_to_cnum(x))
trigonometric_form(x::SymbolicUtils.BasicSymbolic) = _trigonometric_cnum(_to_cnum(x))
trigonometric_form(x::Number) = x
trigonometric_form(x::Complex{Num}) = _trigonometric_cnum(_to_cnum(x))

_sort_key(op::QSym) = (op.space_index, _name_rank(op.index.name_id))
