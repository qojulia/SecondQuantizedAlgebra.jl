# Zero-size tag marking the native fast path (the value lives inline in `slot`).
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
    Coeff

Coefficient representation for operator prefactors. A `Coeff` has three forms: a
native number (a small exact `GaussianRational{Int}` or an inexact `ComplexF64`), a
`Poly` parameter polynomial, and a raw SymbolicUtils fallback. The latter preserves a
single complex expression tree and is lowered to `Complex{Num}` only at public boundaries
(`to_num`).
"""
struct Coeff
    slot::NativeSlot
    tail::Union{Native, Poly, RawSymbolicCoeff}
end
const CNum = Coeff
const RawExpression = SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}

const NUM_ZERO = Num(0)
const NUM_ONE = Num(1)
const EXACT_ONE = one(ExactComplex)
const EXACT_IM = ExactComplex(0, 1)
const EXACT_NEG1 = ExactComplex(-1, 0)
const EMPTY_SYMS = SymbolicUtils.BasicSymbolic[]
const EMPTY_EXPS = Rational{Int}[]

# Adding 0.0+0.0im normalizes any signed zero (`-0.0 -> 0.0`) so that structurally
# equal coefficients (e.g. `conj(2.0)` vs `2.0`) stay `isequal` and hash identically.
@inline native(z::ComplexF64) = Coeff(NativeSlot(normalize_scalar(z)), NATIVE)
@inline native(z::ExactComplex) = Coeff(NativeSlot(z), NATIVE)
@inline native_scalar(c::Coeff) = native_scalar(c.slot)
const ZERO_SLOT = NativeSlot(zero(ExactComplex))
@inline symbolic(
    x::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}; real_slot::Bool = false,
) = Coeff(ZERO_SLOT, RawSymbolicCoeff(x, real_slot))
@inline poly_coeff(p::Poly) = Coeff(ZERO_SLOT, p)
@inline is_native(c::Coeff) = c.tail isa Native
@inline is_poly(c::Coeff) = c.tail isa Poly

const CNUM_ZERO = native(zero(ExactComplex))
const CNUM_ONE = native(EXACT_ONE)
const CNUM_NEG1 = native(EXACT_NEG1)
const CNUM_IM = native(EXACT_IM)
const CNUM_NEG_IM = native(-EXACT_IM)
const CNUM_HALF = native(ExactComplex(1 // 2, 0 // 1))

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
expim(x::Num) = phase_coeff(x)
expim(x::SymbolicUtils.BasicSymbolic) = phase_coeff(x)

@noinline function nonreal_phase_argument(x)
    throw(ArgumentError("`expim` requires a provably real argument; got `$x`"))
end

# A complex number is never accepted merely because its current imaginary part happens to
# be zero. The atom promises unit modulus structurally, so its domain has to be real by type,
# not by a value-dependent test.
expim(x::Number) = nonreal_phase_argument(x)

@inline is_phase(b) =
    b isa SymbolicUtils.BasicSymbolic &&
    SymbolicUtils.iscall(b) &&
    SymbolicUtils.operation(b) === expim

@inline is_imaginary_unit(x) = isequal(x, Symbolics.IM)

@inline function phase_factor_index(syms)::Int
    @inbounds for i in eachindex(syms)
        is_phase(syms[i]) && return i
    end
    return 0
end

# The explicit scalar `shape` matters: without it the term is shaped `Unknown` and every
# later `Complex{Num}` addition it takes part in fails on a shape mismatch.
const SCALAR_SHAPE = UnitRange{Int}[]
raw_complex(re, im) = SymbolicUtils.term(
    complex, SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im);
    type = Complex{Real}, shape = SCALAR_SHAPE,
)
expim_expanded(x) = SymbolicUtils.term(
    expim, SymbolicUtils.unwrap(x); type = Complex{Real}, shape = SCALAR_SHAPE,
)
# `expand` first: `(ω + 2J)*t` and `ω*t + 2J*t` are one phase and must intern to one atom.
expim_symbolic(x) = expim_expanded(expand(x))

# `type` and `shape` take part in hash-consing, so without these every `maketerm` rebuild
# (each `substitute`, each `Postwalk`) recomputes them from the generic fallbacks and mints
# a *different* atom, one whose symtype is `Real` and whose `conj` is therefore the identity.
SymbolicUtils.promote_symtype(::typeof(expim), ::SymbolicUtils.TypeT) = Complex{Real}
SymbolicUtils.promote_shape(::typeof(expim), ::SymbolicUtils.ShapeT) =
    SymbolicUtils.ShapeVecT()

# Without a rule `expand_derivatives` hits the global `nothing` fallback and leaves an inert
# `Differential` node. The body must yield a bare `BasicSymbolic`, hence `expim`.
Symbolics.@register_derivative expim(x) 1 im * expim_symbolic(x)

function leading_sign(v)
    u = SymbolicUtils.unwrap(v)
    if u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(u)
        op = SymbolicUtils.operation(u)
        args = SymbolicUtils.arguments(u)
        (op === (+) && !isempty(args)) && return leading_sign(args[1])
        if op === (*)
            for f in args
                n = const_value(SymbolicUtils.unwrap(f))
                n isa Real && !iszero(n) && return sign(n)
            end
            return 0.0
        end
    end
    n = const_value(u)
    return n isa Real ? sign(n) : 0.0
end

# Never route this through `Symbolics.wrap`: `Num` cannot hold a complex symtype, so wrapping
# splits the phase into `real`/`imag` halves and loses the single factor everything below
# depends on. The argument is oriented so `expim(-x)` interns as `expim(x)` at exponent `-1`;
# without that, two modes rotating at opposite rates carry unrelated atoms and never cancel.
function negate_expanded(a::Num)::Num
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result -= Num(term)
        end
        return result
    end
    return -a
end

function scale_expanded(a::Num, n::Int)::Num
    isone(n) && return a
    n == -1 && return negate_expanded(a)
    raw = SymbolicUtils.unwrap(a)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        result = NUM_ZERO
        for term in SymbolicUtils.arguments(raw)
            result += n * Num(term)
        end
        return result
    end
    return n * a
end

function phase_coeff_expanded(a::Num)
    u = SymbolicUtils.unwrap(a)
    v = const_value(u)
    # A phase over a literal is its value. Interning it instead would keep every coefficient
    # it touches off the native tier, so numerically cancelling terms would never fold.
    v isa Real && return iszero(v) ? CNUM_ONE : native(ComplexF64(exp(im * v)))
    v isa Number && return nonreal_phase_argument(a)
    (u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.symtype(u) <: Real) ||
        return nonreal_phase_argument(a)
    neg = leading_sign(a) < 0
    c = atom_coeff(expim_expanded(neg ? negate_expanded(a) : a))
    return neg ? conj_cnum(c) : c
end

phase_coeff(x) = phase_coeff_expanded(Num(expand(x)))

# Elementary functions with an exact value at argument `0`. Folds the `exp(0)` Symbolics
# leaves after Euler-expanding `exp(im*ω*t)`. Spelled as `===` (not `in` a tuple) to stay
# statically resolved; user-registered functions (`pulse(t)`) have no value here.
@inline is_one_at_zero(op) = op === exp || op === cos || op === cosh
@inline is_zero_at_zero(op) = op === sin || op === tan || op === sinh || op === tanh

# Concrete numeric content of an (unwrapped) symbolic value, else `nothing`
# (`BasicSymbolic` constants from substitute/simplify count, keeping them native).
@inline function const_value(v)
    v isa Number && return v
    (v isa SymbolicUtils.BasicSymbolic && SymbolicUtils.isconst(v)) && return v.val
    return nothing
end
@inline function phase_power(x)
    is_phase(x) && return (only(SymbolicUtils.arguments(x)), 1)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    SymbolicUtils.operation(x) === (^) || return nothing
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return nothing
    is_phase(args[1]) || return nothing
    exponent = const_value(args[2])
    exponent isa Integer || return nothing
    return (only(SymbolicUtils.arguments(args[1])), Int(exponent))
end

# One bottom-up rewrite step for the identities SymbolicUtils cannot infer for the
# package-local `expim` operation. Multiplication is handled as one AC node so phase
# collection is independent of factor order and tree grouping.
function rewrite_phase(x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if length(args) == 1 && is_phase(args[1])
        argument = only(SymbolicUtils.arguments(args[1]))
        op === conj && return expim_expanded(-argument)
        op === real && return cos(argument)
        op === imag && return sin(argument)
        (op === abs || op === abs2) && return 1
    elseif op === (/) && length(args) == 2
        angle = 0
        found_phase = false
        denominator = args[2]
        denominator_phase = phase_power(denominator)
        if denominator_phase !== nothing
            argument, exponent = denominator_phase
            return args[1] * expim_expanded(SymbolicUtils.expand(-exponent * argument))
        elseif SymbolicUtils.iscall(denominator) && SymbolicUtils.operation(denominator) === (*)
            factors = SymbolicUtils.arguments(denominator)
            ordinary_count = 0
            for factor in factors
                phase = phase_power(factor)
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
                return args[1] * expim_expanded(SymbolicUtils.expand(-angle))
            end
            ordinary_denominator = Any[]
            sizehint!(ordinary_denominator, ordinary_count)
            for factor in factors
                phase_power(factor) === nothing && push!(ordinary_denominator, factor)
            end
            denominator = foldl(*, ordinary_denominator; init = 1)
            return (args[1] / denominator) * expim_expanded(SymbolicUtils.expand(-angle))
        end
    elseif op === (*)
        ordinary = Any[]
        angle = 0
        phase_factors = 0
        changed = false
        for factor in args
            phase = phase_power(factor)
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
        value = const_value(expanded_angle)
        phase = value isa Number && iszero(value) ? 1 : expim_expanded(expanded_angle)
        # Build the real amplitude before attaching the complex phase.  Multiplying the
        # phase by an exact real rational first makes SymbolicUtils promote `1//2` to a
        # floating-point complex constant.
        isempty(ordinary) && return phase
        amplitude = foldl(*, ordinary; init = 1)
        phase === 1 && return amplitude
        # Keep the real amplitude in an explicit complex slot while rebuilding sums: a
        # real term and a complex phase term otherwise get promoted together as Float64.
        return raw_complex(
            amplitude::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}, 0 // 1,
        ) * phase
    elseif op === (^) && length(args) == 2 && is_phase(args[1])
        exponent = const_value(args[2])
        exponent isa Integer || return nothing
        argument = only(SymbolicUtils.arguments(args[1]))
        iszero(exponent) && return 1
        return expim_expanded(SymbolicUtils.expand(exponent * argument))
    end
    return nothing
end

const PHASE_NORMALIZER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Fixpoint(SymbolicUtils.Rewriters.Postwalk(rewrite_phase)),
)

normalize_phase(x) = PHASE_NORMALIZER(x)
function simplify_raw(x; kwargs...)
    normalized = normalize_phase(x)
    simplified = SymbolicUtils.simplify(normalized; kwargs...)
    return normalize_phase(simplified)
end

function from_raw(x; normalize::Bool = true, real_slot::Bool = false)::Coeff
    value = const_value(x)
    value isa Number && return to_cnum(value)
    x isa SymbolicUtils.BasicSymbolic || return to_cnum(x)
    expr = normalize ? normalize_phase(x) : x
    value = const_value(expr)
    value isa Number && return to_cnum(value)
    # Complex-valued intermediate expressions can acquire an explicit zero imaginary
    # slot while expanding exact rational coefficients.  Canonicalize that case back to
    # the real symbolic tree so equivalent expressions do not differ only by `complex(..., 0)`.
    real_part, imaginary_part = raw_realimag(expr)
    imaginary_value = const_value(imaginary_part)
    if imaginary_value isa Number && iszero(imaginary_value)
        real_expr = real_part isa SymbolicUtils.BasicSymbolic ? real_part :
            SymbolicUtils.unwrap(Num(real_part))
        real_value = const_value(real_expr)
        real_value isa Number && return to_cnum(real_value)
        return symbolic(real_expr; real_slot = true)
    end
    is_phase(expr) && return phase_coeff(only(SymbolicUtils.arguments(expr)))
    return symbolic(expr; real_slot)
end

@inline from_raw_arithmetic(
    x::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}, real_slot::Bool,
)::Coeff = symbolic(x; real_slot)

@inline num_from_scalar(x::Int) = Num(x)
@inline num_from_scalar(x::Float64) = Num(x)
@inline num_from_scalar(x::Rational{Int}) =
    denominator(x) == 1 ? Num(Int(numerator(x))) : Num(x)
function num_from_scalar(x::Rational{BigInt})
    fits_small(x) && return num_from_scalar(Rational{Int}(x))
    return denominator(x) == 1 ? Num(numerator(x)) : Num(x)
end

@inline function constant_term(::Type{E}, z)::Monomial{E} where {E <: ExactScalar}
    return Monomial{E}(as_tier(E, z), EMPTY_SYMS, EMPTY_EXPS)
end

@inline scalar_coeff(x::NativeScalar)::Coeff = native(x)
function scalar_coeff(x::BigExactComplex)::Coeff
    fits_small(x) && return native(as_tier(ExactComplex, x))
    return poly_coeff(Poly(Monomial{BigExactComplex}[constant_term(BigExactComplex, x)]))
end

@inline exact_value(::Type{E}, re::T, im::T) where {E <: ExactScalar, T <: Union{Int, Rational{Int}}} =
    E(re, im)
@inline exact_rational(x::Integer) = Rational{BigInt}(x)
@inline exact_rational(x::Rational) = Rational{BigInt}(x)
@inline exact_complex(re::Union{Integer, Rational}, im::Union{Integer, Rational}) =
    scalar_coeff(BigExactComplex(exact_rational(re), exact_rational(im)))

to_cnum(x::Coeff) = x
to_cnum(x::Num) = recognize(SymbolicUtils.unwrap(x))
to_cnum(x::Rational{Int}) = isinf(x) ? native(ComplexF64(x)) : tiered(exact_value, x, 0 // 1)
to_cnum(x::Complex{Int}) = tiered(exact_value, real(x), imag(x))
to_cnum(x::ExactComplex) = scalar_coeff(x)
@inline to_cnum(x::Int) = is_edge(x) ? exact_complex(x, 0) : native(unsafe_gaussian(x, 0, 1))
to_cnum(x::Union{Bool, Int8, Int16, Int32, UInt8, UInt16, UInt32}) = to_cnum(Int(x))
to_cnum(x::Complex{Bool}) = native(ExactComplex(Int(real(x)), Int(imag(x))))
to_cnum(x::Integer) = exact_complex(x, 0)
to_cnum(x::Rational) = isinf(x) ? native(ComplexF64(x)) : exact_complex(x, 0)
function to_cnum(x::Complex{<:Union{Integer, Rational}})
    (isinf(real(x)) || isinf(imag(x))) && return native(ComplexF64(x))
    return exact_complex(real(x), imag(x))
end
# Native only when the value round-trips through ComplexF64 with no loss; non-rational
# values that cannot be represented faithfully (for example, large bignums) stay symbolic.
function to_cnum(x::Real)
    z = ComplexF64(x)
    return z == x ? native(z) : symbolic(SymbolicUtils.unwrap(Num(x)))
end
function to_cnum(x::Complex)
    z = ComplexF64(x)
    z == x && return native(z)
    raw_re = SymbolicUtils.unwrap(Num(real(x)))
    raw_im = SymbolicUtils.unwrap(Num(imag(x)))
    return symbolic(raw_re + im * raw_im)
end
function to_cnum(x::Complex{Num})
    re, im = real(x), imag(x)
    raw_re, raw_im = SymbolicUtils.unwrap(re), SymbolicUtils.unwrap(im)
    im_value = const_value(raw_im)
    if is_phase(raw_re) && im_value isa Number && iszero(im_value)
        return phase_coeff(only(SymbolicUtils.arguments(raw_re)))
    end
    # The fields of `Complex{Num}` are explicit real and imaginary coefficient slots,
    # even when a contained symbolic atom has the conservative `Number` symtype.
    # Canonicalize both through the coefficient algebra so materializing and then
    # rebuilding a coefficient does not switch between Raw and Poly tiers.
    return cnum(re, im)
end
to_cnum(x::SymbolicUtils.BasicSymbolic) = recognize(x)


# Canonicalizing constructor from real/imag `Num` parts (`re + im*i`), used by the
# symbolic boundaries (substitute / conj / change_index) that may yield a polynomial.
function mul_by_im(c::Coeff)::Coeff
    tail = c.tail
    tail isa Native && return mul_native(native_scalar(c), EXACT_IM)
    tail isa Poly && return tiered(poly_scale, tail, EXACT_IM)
    raw = (im * tail.expr)::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    return from_raw_arithmetic(raw, false)
end

cnum(re::Num, im::Num) =
    add_cnum(recognize(SymbolicUtils.unwrap(re)), mul_by_im(recognize(SymbolicUtils.unwrap(im))))

# === Parameter-polynomial tier: folding, materialization, recognition ===

# Polynomial monomial scalars retain exact Gaussian rationals; only an explicitly
# floating-point input or an operation involving one enters the ComplexF64 path.

# Fold a canonical term list into a Coeff: empty -> zero, one constant term ->
# native, else a Poly. Inputs are already canonical, so no re-sort here.
function from_poly(terms::Vector{Monomial{ExactComplex}})
    isempty(terms) && return CNUM_ZERO
    length(terms) == 1 && isempty(terms[1].syms) && return scalar_coeff(terms[1].scalar)
    return poly_coeff(Poly(terms))
end
function from_poly(terms::Vector{Monomial{BigExactComplex}})
    all(fits_small, terms) && return from_poly(as_tier(ExactComplex, terms))
    return poly_coeff(Poly(terms))
end
@inline function fits_small(m::Monomial{BigExactComplex})
    scalar = m.scalar
    return scalar isa ComplexF64 || fits_small(scalar)
end

@inline tier_result(terms::Vector{<:Monomial}) = from_poly(terms)
@inline tier_result(z::CoeffScalar) = scalar_coeff(z)

@inline is_big(::Vector{Monomial{BigExactComplex}}) = true
@inline is_big(::BigExactComplex) = true
@inline is_big(c::Coeff) = c.tail isa Poly && !is_small(c.tail)
@inline is_big(parts::Vector{Coeff}) = any(is_big, parts)
@inline is_big(_) = false

@inline function tiered(f::F, a::A)::Coeff where {F, A}
    is_big(a) && return big_tier(f, a)
    try
        return tier_result(f(ExactComplex, a))
    catch err
        err isa OverflowError || rethrow()
        return big_tier(f, a)
    end
end
@inline function tiered(f::F, a::A, b::B)::Coeff where {F, A, B}
    (is_big(a) || is_big(b)) && return big_tier(f, a, b)
    try
        return tier_result(f(ExactComplex, a, b))
    catch err
        err isa OverflowError || rethrow()
        return big_tier(f, a, b)
    end
end
@inline function tiered(f::F, a::A, b::B, c::C)::Coeff where {F, A, B, C}
    (is_big(a) || is_big(b) || is_big(c)) && return big_tier(f, a, b, c)
    try
        return tier_result(f(ExactComplex, a, b, c))
    catch err
        err isa OverflowError || rethrow()
        return big_tier(f, a, b, c)
    end
end
@noinline function big_tier(f::F, a::A)::Coeff where {F, A}
    return tier_result(f(BigExactComplex, as_tier(BigExactComplex, a)))
end
@noinline function big_tier(f::F, a::A, b::B)::Coeff where {F, A, B}
    E = BigExactComplex
    return tier_result(f(E, as_tier(E, a), as_tier(E, b)))
end
@noinline function big_tier(f::F, a::A, b::B, c::C)::Coeff where {F, A, B, C}
    E = BigExactComplex
    return tier_result(f(E, as_tier(E, a), as_tier(E, b), as_tier(E, c)))
end

@inline is_small(p::Poly) = p.terms isa Vector{Monomial{ExactComplex}}
@inline small_terms(p::Poly) = p.terms::Vector{Monomial{ExactComplex}}
@inline function big_terms(p::Poly)::Vector{Monomial{BigExactComplex}}
    terms = p.terms
    terms isa Vector{Monomial{ExactComplex}} && return as_tier(BigExactComplex, terms)
    return terms
end
@inline function poly_terms(::Type{ExactComplex}, p::Poly)::Vector{Monomial{ExactComplex}}
    is_small(p) || throw(OverflowError("exact scalar exceeds the small tier"))
    return small_terms(p)
end
@inline poly_terms(::Type{BigExactComplex}, p::Poly) = big_terms(p)

@inline function tiered(f::F, a::Poly)::Coeff where {F}
    is_small(a) || return tier_result(f(BigExactComplex, big_terms(a)))
    return tiered(f, small_terms(a))
end
@inline function tiered(f::F, a::Poly, b::B)::Coeff where {F, B}
    is_small(a) || return big_tier(f, big_terms(a), b)
    return tiered(f, small_terms(a), b)
end
@inline function tiered(f::F, a::Poly, b::Poly)::Coeff where {F}
    (is_small(a) && is_small(b)) || return big_tier(f, big_terms(a), big_terms(b))
    return tiered(f, small_terms(a), small_terms(b))
end
@inline function tiered(f::F, a::Poly, b::B, c::C)::Coeff where {F, B, C}
    is_small(a) || return big_tier(f, big_terms(a), b, c)
    return tiered(f, small_terms(a), b, c)
end

function radical_factors(m::Monomial)
    factors = Tuple{Rational{Int}, Int}[]
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        is_radical_atom(s) || continue
        e = m.exps[i]
        p = s.val::Int
        k = 0
        for j in eachindex(factors)
            if factors[j][1] == e
                k = j
                break
            end
        end
        if k == 0
            push!(factors, (e, p))
        else
            product, overflow = Base.mul_with_overflow(factors[k][2], p)
            overflow ? push!(factors, (e, p)) : (factors[k] = (e, product))
        end
    end
    return insertion_sort!(factors, isless)
end

function radical_expression(n::Int, f::Rational{Int})::RawExpression
    base = SymbolicUtils.Const{SymbolicUtils.SymReal}(n)
    f == 1 // 2 && return sqrt(base)
    f == 1 // 3 && return cbrt(base)
    return SymbolicUtils.term(^, base, f; type = Real)
end

# One monomial term -> Complex{Num}. The integer part of each exponent is built by
# repeated multiply/divide (both infer `Num`, unlike `Num ^ Int` which infers `Any`).
function term_to_num(m::Monomial)
    prod = NUM_ONE
    @inbounds for i in eachindex(m.syms)
        s = m.syms[i]
        e = m.exps[i]
        is_radical_atom(s) && continue
        # `exp(-im*x)` rather than `1 / exp(im*x)`: same value, and it keeps an inverse
        # phase readable. Re-recognising it lands back on this atom, since `recognize`
        # re-orients every `expim`.
        if e < 0 && is_phase(s)
            s = expim_symbolic(-only(SymbolicUtils.arguments(s)))
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
    for (f, n) in radical_factors(m)
        prod = prod * Num(radical_expression(n, f))
    end
    # Guard the zero halves: `0 * x` does not always fold for a complex-symtype factor,
    # and an unfolded `0*expim(...)` would survive all the way to display.
    scalar = m.scalar
    real_scalar, imag_scalar = real(scalar), imag(scalar)
    re = iszero(real_scalar) ? NUM_ZERO : num_from_scalar(real_scalar) * prod
    imag_ = iszero(imag_scalar) ? NUM_ZERO : num_from_scalar(imag_scalar) * prod
    return Complex(re, imag_)
end

# Sum the terms; the only place a polynomial lowers to SymbolicUtils.
poly_to_num(p::Poly) = terms_to_num(p.terms)
function terms_to_num(terms::Vector{<:Monomial})
    isempty(terms) && return Complex(NUM_ZERO, NUM_ZERO)
    acc = term_to_num(terms[1])
    @inbounds for i in 2:length(terms)
        acc = acc + term_to_num(terms[i])
    end
    return acc
end

@inline function with_raw_number(g::G, z::ExactComplex)::RawExpression where {G}
    if iszero(z.im)
        isone(z.den) && return g(z.re)::RawExpression
        return g(z.re // z.den)::RawExpression
    end
    isone(z.den) && return g(Complex(z.re, z.im))::RawExpression
    return g(Complex(z.re // z.den, z.im // z.den))::RawExpression
end
@inline function with_raw_number(g::G, z::BigExactComplex)::RawExpression where {G}
    fits_small(z) && return with_raw_number(g, as_tier(ExactComplex, z))
    iszero(z.im) && return g(real(z))::RawExpression
    return g(Complex(real(z), imag(z)))::RawExpression
end
@inline with_raw_number(g::G, z::ComplexF64) where {G} = g(z)::RawExpression

@inline multiply_factors(factors::Vector{RawExpression}) =
    reduce((x, y) -> (x * y)::RawExpression, factors)::RawExpression

function raw_power(symbol::RawExpression, exponent::Rational{Int})::RawExpression
    denominator(exponent) == 1 && return symbol^numerator(exponent)
    return symbol^exponent
end

function factors_to_raw(m::Monomial)::RawExpression
    factors = RawExpression[]
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]::RawExpression
        is_radical_atom(symbol) && continue
        # The polynomial convention for a non-real-symtype atom is carried by the
        # `real_slot` bit on the raw coefficient. Do not wrap the factor here: SymbolicUtils
        # simplifies `complex(z, 0)` back to `z` before the bit can be observed.
        push!(factors, raw_power(symbol, m.exps[i]))
    end
    for (f, n) in radical_factors(m)
        push!(factors, radical_expression(n, f))
    end
    return multiply_factors(factors)
end

function term_to_raw(m::Monomial)::RawExpression
    factors = factors_to_raw(m)
    isone(m.scalar) && return factors
    return with_raw_number(v -> v * factors, m.scalar)
end

is_constant_term(m::Monomial) = isempty(m.syms)

raw_sum(terms::Vector{Monomial{E}}, range::AbstractUnitRange{Int}) where {E} =
    reduce((x, y) -> (x + y)::RawExpression, (term_to_raw(terms[i]) for i in range))::RawExpression

poly_to_raw(p::Poly) = terms_to_raw(p.terms)
function terms_to_raw(terms::Vector{Monomial{E}})::RawExpression where {E}
    is_constant_term(first(terms)) || return raw_sum(terms, eachindex(terms))
    constant = first(terms).scalar
    length(terms) == 1 &&
        return with_raw_number(SymbolicUtils.Const{SymbolicUtils.SymReal}, constant)
    variable = raw_sum(terms, 2:lastindex(terms))
    return with_raw_number(v -> variable + v, constant)
end

@inline function raw_tail(c::Coeff)::RawExpression
    tail = c.tail
    tail isa Poly && return poly_to_raw(tail)
    tail isa RawSymbolicCoeff && return tail.expr
    return with_raw_number(SymbolicUtils.Const{SymbolicUtils.SymReal}, native_scalar(c))
end

@inline function raw_binary(f::F, a::Coeff, b::Coeff)::RawExpression where {F}
    if is_native(a)
        y = raw_tail(b)
        return with_raw_number(v -> f(v, y), native_scalar(a))
    elseif is_native(b)
        x = raw_tail(a)
        return with_raw_number(v -> f(x, v), native_scalar(b))
    end
    return f(raw_tail(a), raw_tail(b))::RawExpression
end


function raw_realimag(x)
    x isa Number && return (real(x), imag(x))
    value = const_value(x)
    value isa Number && return (real(value), imag(value))
    is_imaginary_unit(x) && return (0, 1)
    SymbolicUtils.symtype(x) <: Real && return (x, 0)
    is_phase(x) && begin
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
                ar, ai = raw_realimag(argument)
                re += ar
                im += ai
            end
            return (re, im)
        elseif op === (*)
            re, im = 1, 0
            for argument in args
                ar, ai = raw_realimag(argument)
                re, im = re * ar - im * ai, re * ai + im * ar
            end
            return (re, im)
        elseif op === (/) && length(args) == 2 &&
                SymbolicUtils.symtype(args[2]) <: Real
            re, im = raw_realimag(args[1])
            return (re / args[2], im / args[2])
        elseif op === conj
            re, im = raw_realimag(only(args))
            return (re, -im)
        end
    end
    # An opaque Number-symtype expression may itself be complex-valued. Preserve the
    # public `real`/`imag` semantics here; the polynomial tier makes its real-slot
    # convention explicit in `term_to_raw` before this fallback is reached.
    return (real(x), imag(x))
end

poly_is_real(p::Poly) = terms_are_real(p.terms)
function terms_are_real(terms::Vector{<:Monomial})
    for monomial in terms
        scalar_is_real(monomial) || return false
        for symbol in monomial.syms
            (is_phase(symbol) || is_imaginary_unit(symbol)) && return false
        end
    end
    return true
end

@inline raw_is_real(t::RawSymbolicCoeff) =
    t.real_slot || SymbolicUtils.symtype(t.expr) <: Real

@inline function cnum_is_real(c::Coeff)
    t = c.tail
    t isa Native && return iszero(imag(native_scalar(c)))
    t isa Poly && return poly_is_real(t)
    return raw_is_real(t)
end

# An "atom" is an irreducible scalar the polynomial tier treats as one opaque
# variable: a symbol, an array index (`ω[i]`), or a non-algebraic one-arg call on an
# atom (`real(g)`, `imag(g)`, `sqrt`, `exp`, `conj`, ...). Algebraic ops (`+ * ^ /`,
# `complex`) are decomposed by `recognize` instead, keeping their structure native.
@inline function is_atom(b)
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
    return length(args) == 1 && is_atom(only(args))
end

# Keep the identity scalar native for ordinary symbolic products; explicit rational
# inputs enter through `scalar_coeff` and remain exact when they are combined later.
# A bare atom (symbol / array index / `conj(atom)`) as a single-monomial Coeff.
@inline atom_coeff(x::SymbolicUtils.BasicSymbolic) = poly_coeff(
    Poly(
        Monomial{ExactComplex}[
            Monomial{ExactComplex}(EXACT_ONE, SymbolicUtils.BasicSymbolic[x], Rational{Int}[1]),
        ],
    ),
)
# An unrecognized symbolic value, kept as one raw symbolic expression tree.
@inline sym_leaf(x::SymbolicUtils.BasicSymbolic) = from_raw(x; normalize = false)

function integer_root(n::Integer, q::Int)::BigInt
    target = BigInt(n)
    guess = round(BigInt, BigFloat(target)^(one(BigFloat) / q))
    for m in (guess - 1, guess, guess + 1)
        m >= 0 && m^q == target && return m
    end
    return guess
end
@inline is_integer_power(n::Integer, q::Int) = integer_root(n, q)^q == n

function is_rational_power(value::Rational{BigInt}, r::Rational{Int})
    q = denominator(r)
    (is_integer_power(numerator(value), q) && is_integer_power(denominator(value), q)) ||
        return false
    return !(iszero(value) && numerator(r) < 0)
end
function rational_power_value(value::Rational{BigInt}, r::Rational{Int})
    q = denominator(r)
    return (integer_root(numerator(value), q) // integer_root(denominator(value), q))^numerator(r)
end

@inline function monomial_terms(
        ::Type{E}, scalar::CoeffScalar, syms::Vector{SymbolicUtils.BasicSymbolic},
        exps::Vector{Rational{Int}},
    )::Vector{Monomial{E}} where {E <: ExactScalar}
    return Monomial{E}[Monomial{E}(as_tier(E, scalar), syms, exps)]
end

function numeric_radical(value::Rational{BigInt}, r::Rational{Int}, x)::Coeff
    is_rational_power(value, r) && return to_cnum(rational_power_value(value, r))
    value > 0 || return radical_leaf(x)
    numerator_factors = Tuple{Int, Int}[]
    denominator_factors = Tuple{Int, Int}[]
    prime_factorization!(numerator_factors, numerator(value)) || return radical_leaf(x)
    prime_factorization!(denominator_factors, denominator(value)) || return radical_leaf(x)
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    for (p, k) in numerator_factors
        push!(syms, SymbolicUtils.Const{SymbolicUtils.SymReal}(p)); push!(exps, k * r)
    end
    for (p, k) in denominator_factors
        push!(syms, SymbolicUtils.Const{SymbolicUtils.SymReal}(p)); push!(exps, -k * r)
    end
    sort_factors!(syms, exps)
    return tiered(monomial_terms, EXACT_ONE, syms, exps)
end

@inline radical_leaf(x::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}) =
    symbolic(x; real_slot = true)

@inline function is_exact_number(x)
    value = const_value(x)
    return value isa Integer || value isa Rational
end
@inline exact_number(x) = Rational{BigInt}(const_value(x)::Union{Integer, Rational})

# A fractional power `base^r`. Native only for a floating-point base or a single-atom
# unit-scalar monomial (giving that atom a rational exponent); any other base would
# need to distribute the radical (unsound), so it becomes a symbolic leaf.
function negative_radical(value::Rational{BigInt}, r::Rational{Int}, x, real_root::Bool)::Coeff
    magnitude = numeric_radical(-value, r, x)
    if magnitude.tail isa RawSymbolicCoeff
        real_root && return radical_leaf(x)
        return native(ComplexF64(Float64(value))^r)
    end
    unit = real_root ? EXACT_NEG1 : EXACT_IM^numerator(r)
    tail = magnitude.tail
    tail isa Poly && return tiered(poly_scale, tail, unit)
    return mul_native(native_scalar(magnitude), unit)
end

function rational_power(basearg, r::Rational{Int}, x)
    if is_exact_number(basearg)
        value = exact_number(basearg)
        value >= 0 && return numeric_radical(value, r, x)
        denominator(r) == 2 && return negative_radical(value, r, x, false)
    end
    base = recognize(basearg)
    is_native(base) && return native(native_float(native_scalar(base))^r)
    if base.tail isa Poly && length(base.tail.terms) == 1
        m = only(base.tail.terms)
        if length(m.syms) == 1 && isone(m.scalar)
            is_phase(only(m.syms)) &&
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            return tiered(monomial_terms, EXACT_ONE, m.syms, Rational{Int}[m.exps[1] * r])
        end
    end
    return sym_leaf(x)
end

# Total recognizer: evaluate a symbolic expression tree in the Coeff algebra.
# Numbers/atoms map to native/Poly; `+ * ^(int) /` compose through the coefficient
# arithmetic; radicals fold to rational exponents; everything else is a symbolic
# leaf. Always returns a `Coeff` (no `nothing` sentinel).
recognize(x::Number)::Coeff = to_cnum(x)
recognize(x::Num)::Coeff = recognize(SymbolicUtils.unwrap(x))
recognize(x)::Coeff = to_cnum(x)
function recognize(x::SymbolicUtils.BasicSymbolic)::Coeff
    SymbolicUtils.isconst(x) && return recognize(x.val)
    is_imaginary_unit(x) && return CNUM_IM
    SymbolicUtils.issym(x) && return atom_coeff(x)
    SymbolicUtils.iscall(x) || return sym_leaf(x)
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        return recognize_sum(args)
    elseif op === (*)
        return recognize_prod(args)
    elseif op === (^)
        length(args) == 2 || return sym_leaf(x)
        # `isa` guards narrow the `Any` exponent before arithmetic (a helper call on
        # the `Any` value would force a runtime dispatch).
        pv = const_value(args[2])
        if pv isa Integer
            return pow_cnum_integer(recognize(args[1]), Int(pv))
        elseif pv isa Rational{Int}
            return rational_power(args[1], pv, x)
        end
        return sym_leaf(x)
    elseif op === getindex
        return atom_coeff(x)
    elseif op === expim
        # Through `phase_coeff`, not `atom_coeff`: a phase recovered from a lowered
        # expression has to orient the same way as one built directly, or the two spellings
        # become unrelated atoms and stop cancelling.
        length(args) == 1 && return phase_coeff(only(args))
        return sym_leaf(x)
    elseif op === conj
        # Fold conj of a constant (e.g. conj(0.0) left by substituting a complex
        # parameter to a real value) so the coefficient can collapse to zero.
        inner = recognize(args[1])
        is_native(inner) && return native(conj(native_scalar(inner)))
        # Without this, `conj(expim(x))` interns as a fresh atom that never cancels against
        # the phase it came from.
        is_phase(args[1]) && return conj_cnum(inner)
        return is_atom(x) ? atom_coeff(x) : sym_leaf(x)
    elseif op === (/)
        length(args) == 2 || return sym_leaf(x)
        numerator = recognize(args[1])
        denominator = const_value(args[2])
        if denominator isa Integer
            return mul_cnum(numerator, to_cnum(1 // denominator))
        end
        return numerator / recognize(args[2])
    elseif op === complex
        length(args) == 2 && return cnum(Num(args[1]), Num(args[2]))
        return sym_leaf(x)
    elseif op === sqrt
        length(args) == 1 && return rational_power(only(args), 1 // 2, x)
        return sym_leaf(x)
    elseif op === cbrt
        length(args) == 1 || return sym_leaf(x)
        radicand = only(args)
        if is_exact_number(radicand) && exact_number(radicand) < 0
            return negative_radical(exact_number(radicand), 1 // 3, x, true)
        end
        return rational_power(radicand, 1 // 3, x)
    end
    # Fold an elementary function of a literal zero (`exp(0) -> 1`, `sin(0) -> 0`, ...).
    if length(args) == 1 && (is_one_at_zero(op) || is_zero_at_zero(op))
        a1 = recognize(only(args))
        (is_native(a1) && iszero(native_scalar(a1))) && return is_one_at_zero(op) ? CNUM_ONE : CNUM_ZERO
    end
    # Irreducible one-arg call on an atom (`exp`, `sin`, `real`, `imag`, ...): keep it
    # native as an opaque integer-exponent atom. Radicals are handled above instead.
    return is_atom(x) ? atom_coeff(x) : sym_leaf(x)
end

function recognize_sum(args)::Coeff
    parts = Coeff[]
    znative = CNUM_ZERO
    sym = CNUM_ZERO
    have_sym = false
    for a in args
        ca = recognize(a)
        t = ca.tail
        if t isa Native
            znative = add_cnum(znative, ca)
        elseif t isa Poly
            push!(parts, ca)
        else
            sym = add_cnum(sym, ca)
            have_sym = true
        end
    end
    iszero(znative) || push!(parts, znative)
    poly = tiered(summed_terms, parts)
    return have_sym ? add_cnum(poly, sym) : poly
end

function summed_terms(::Type{E}, parts::Vector{Coeff}) where {E}
    terms = Monomial{E}[]
    for c in parts
        t = c.tail
        if t isa Native
            push!(terms, constant_term(E, native_scalar(c)))
        else
            append!(terms, poly_terms(E, t::Poly))
        end
    end
    return canonical_terms!(terms)
end

# Insertion-sort a factor list (syms + matching exps) in place by objectid key.
# Avoids Base's `sortperm` machinery for the abstract `BasicSymbolic` eltype, whose
# inference is a nontrivial slice of first-call latency; these lists are tiny.
function sort_factors!(syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}})
    @inbounds for i in 2:length(syms)
        s = syms[i]; e = exps[i]; k = fkey(s)
        j = i - 1
        while j >= 1 && fkey(syms[j]) > k
            syms[j + 1] = syms[j]; exps[j + 1] = exps[j]
            j -= 1
        end
        syms[j + 1] = s; exps[j + 1] = e
    end
    return nothing
end

function merge_factor_list(syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}})
    n = length(syms)
    n <= 1 && return (syms, exps)
    sort_factors!(syms, exps)
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

@noinline function canonical_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        syms::Vector{SymbolicUtils.BasicSymbolic},
        exps::Vector{Rational{Int}},
    )::Monomial{E} where {E <: ExactScalar}
    ordinary_syms = SymbolicUtils.BasicSymbolic[]
    ordinary_exps = Rational{Int}[]
    sizehint!(ordinary_syms, length(syms) - 1)
    sizehint!(ordinary_exps, length(exps) - 1)
    angle = NUM_ZERO
    @inbounds for i in eachindex(syms)
        symbol = syms[i]
        exponent = exps[i]
        if is_phase(symbol)
            denominator(exponent) == 1 ||
                throw(ArgumentError("a unit phase cannot have a fractional power"))
            argument = Num(only(SymbolicUtils.arguments(symbol)))
            angle += scale_expanded(argument, numerator(exponent))
        else
            push!(ordinary_syms, symbol)
            push!(ordinary_exps, exponent)
        end
    end
    ordinary_syms, ordinary_exps = merge_factor_list(ordinary_syms, ordinary_exps)
    phase = phase_coeff_expanded(angle)
    if phase.tail isa Native
        return Monomial{E}(scalar_mul(scalar, as_tier(E, native_scalar(phase))), ordinary_syms, ordinary_exps)
    end
    phase_term = only(poly_terms(E, phase.tail::Poly))
    append!(ordinary_syms, phase_term.syms)
    append!(ordinary_exps, phase_term.exps)
    sort_factors!(ordinary_syms, ordinary_exps)
    return Monomial{E}(
        scalar_mul(scalar, phase_term.scalar), ordinary_syms, ordinary_exps,
    )
end

@noinline function scaled_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        symbol::Num,
        exponent::Rational{Int},
    )::Monomial{E} where {E <: ExactScalar}
    denominator(exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    argument = Num(only(SymbolicUtils.arguments(SymbolicUtils.unwrap(symbol))))
    phase = phase_coeff_expanded(scale_expanded(argument, numerator(exponent)))
    phase.tail isa Native && return constant_term(E, scalar_mul(scalar, as_tier(E, native_scalar(phase))))
    phase_term = only(poly_terms(E, phase.tail::Poly))
    return Monomial{E}(
        scalar_mul(scalar, phase_term.scalar),
        phase_term.syms,
        phase_term.exps,
    )
end

@noinline function merged_phase_monomial(
        ::Type{E},
        scalar::Union{ComplexF64, E},
        left::Num,
        left_exponent::Rational{Int},
        right::Num,
        right_exponent::Rational{Int},
    )::Monomial{E} where {E <: ExactScalar}
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
    angle = scale_expanded(left_argument, left_power) +
        scale_expanded(right_argument, right_power)
    phase = phase_coeff_expanded(angle)
    phase.tail isa Native && return constant_term(E, scalar_mul(scalar, as_tier(E, native_scalar(phase))))
    phase_term = only(poly_terms(E, phase.tail::Poly))
    exponents = common_negative ? -phase_term.exps : phase_term.exps
    return Monomial{E}(
        scalar_mul(scalar, phase_term.scalar),
        phase_term.syms,
        exponents,
    )
end

function recognize_prod(args)::Coeff
    factors = Coeff[]
    other = CNUM_ONE
    have_other = false
    for a in args
        ca = recognize(a)
        t = ca.tail
        if t isa Native || (t isa Poly && length(t.terms) == 1)
            push!(factors, ca)
        else
            other = mul_cnum_slow(other, ca)
            have_other = true
        end
    end
    mono = tiered(product_terms, factors)
    return have_other ? mul_cnum_slow(mono, other) : mono
end

function product_terms(::Type{E}, factors::Vector{Coeff}) where {E}
    scalar::Union{ComplexF64, E} = one(E)
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    for c in factors
        t = c.tail
        if t isa Native
            scalar = scalar_mul(scalar, as_tier(E, native_scalar(c)))
        else
            m = only(poly_terms(E, t::Poly))
            scalar = scalar_mul(scalar, m.scalar)
            append!(syms, m.syms)
            append!(exps, m.exps)
        end
    end
    iszero(scalar) && return Monomial{E}[]
    count(is_phase, syms) <= 1 || return Monomial{E}[canonical_phase_monomial(E, scalar, syms, exps)]
    ms, me = merge_factor_list(syms, exps)
    return Monomial{E}[Monomial{E}(normalize_scalar(scalar), ms, me)]
end

Base.convert(::Type{Coeff}, x::Coeff) = x
Base.convert(::Type{Coeff}, x::Complex{Num}) = to_cnum(x)
Base.convert(::Type{Coeff}, x::Number) = to_cnum(x)

# `Matrix{Coeff}(I, n, n)` builds its diagonal through the type, not through `convert`.
Coeff(x::Number) = to_cnum(x)

# LinearAlgebra routes its scalar fallbacks through `::Number`, which `Coeff` is not, so each
# hook a coefficient matrix or vector reaches is supplied here.

LinearAlgebra.dot(a::Coeff, b::Coeff) = mul_cnum(conj_cnum(a), b)
LinearAlgebra.dot(a::Coeff, b::Number) = LinearAlgebra.dot(a, to_cnum(b))
LinearAlgebra.dot(a::Number, b::Coeff) = LinearAlgebra.dot(to_cnum(a), b)

# `abs` is exact for native and pure-phase coefficients and throws otherwise, so a norm over
# symbolic coefficients refuses rather than inventing a magnitude it cannot order.
function LinearAlgebra.norm(c::Coeff, p::Real = 2)
    p == 0 && return iszero(c) ? NUM_ZERO : NUM_ONE
    return abs(c)
end

LinearAlgebra.symmetric(c::Coeff, ::Symbol = :U) = c
LinearAlgebra.symmetric_type(::Type{Coeff}) = Coeff
# A Hermitian matrix has a real diagonal, which is what the scalar hook enforces.
LinearAlgebra.hermitian(c::Coeff, ::Symbol = :U) = to_cnum(real(c))
LinearAlgebra.hermitian_type(::Type{Coeff}) = Coeff

LinearAlgebra.rmul!(A::AbstractArray{Coeff}, b::Coeff) = (A .= A .* b)
LinearAlgebra.lmul!(a::Coeff, B::AbstractArray{Coeff}) = (B .= a .* B)

# A symbolic coefficient has no magnitude order, so the pivoted LU `det` would otherwise fall
# back to cannot run. Expand by minors instead, as Symbolics does for `Num`.
function LinearAlgebra.det(A::AbstractMatrix{Coeff})
    LinearAlgebra.checksquare(A)
    rows, cols = axes(A)
    isempty(rows) && return CNUM_ONE
    length(rows) == 1 && return A[first(rows), first(cols)]
    if istriu(A) || istril(A)
        return foldl((acc, ij) -> acc * A[ij[1], ij[2]], zip(rows, cols); init = CNUM_ONE)
    end
    top = first(rows)
    rest = rows[(begin + 1):end]
    acc = CNUM_ZERO
    for (k, j) in enumerate(cols)
        term = A[top, j] * LinearAlgebra.det(@view A[rest, filter(!=(j), cols)])
        acc = isodd(k) ? acc + term : acc - term
    end
    return acc
end

@inline unit_modulus(z::ComplexF64) = abs2(z) == 1.0
@inline unit_modulus(z::ExactComplex) =
    widemul(z.re, z.re) + widemul(z.im, z.im) == widemul(z.den, z.den)
@inline unit_modulus(z::BigExactComplex) = isone(abs2(z))

function is_pure_phase(c::Coeff)::Bool
    tail = c.tail
    tail isa Poly || return false
    length(tail.terms) == 1 || return false
    monomial = only(tail.terms)
    (length(monomial.syms) == 1 && is_phase(only(monomial.syms))) || return false
    denominator(only(monomial.exps)) == 1 || return false
    return unit_modulus(monomial.scalar)
end

function pure_phase_parts(c::Coeff)::Tuple{CoeffScalar, Num, Bool}
    monomial = only((c.tail::Poly).terms)
    symbol = only(monomial.syms)
    exponent = only(monomial.exps)
    scalar = monomial.scalar
    argument = Num(only(SymbolicUtils.arguments(symbol)))
    angle = scale_expanded(argument, numerator(exponent))
    flip_sine = false
    if leading_sign(angle) < 0
        angle = negate_expanded(angle)
        flip_sine = true
    end
    return (scalar, angle, flip_sine)
end

@inline function scaled_trig(a::Real, trig, angle::Num)::Num
    iszero(a) && return NUM_ZERO
    isone(a) && return trig(angle)
    a == -1 && return -trig(angle)
    return num_from_scalar(a) * trig(angle)
end

function Base.real(c::Coeff)::Num
    is_native(c) && return num_from_scalar(real(native_scalar(c)))
    if !is_pure_phase(c)
        c.tail isa Poly && return real(poly_to_num(c.tail))
        expr = normalize_phase(raw_tail(c))
        re, _ = cnum_is_real(c) ? (expr, 0) : raw_realimag(expr)
        return Num(normalize_phase(re))
    end
    scalar, angle, flip_sine = pure_phase_parts(c)
    sine = flip_sine ? -imag(scalar) : imag(scalar)
    return scaled_trig(real(scalar), cos, angle) - scaled_trig(sine, sin, angle)
end

function Base.imag(c::Coeff)::Num
    is_native(c) && return num_from_scalar(imag(native_scalar(c)))
    if !is_pure_phase(c)
        c.tail isa Poly && return imag(poly_to_num(c.tail))
        expr = normalize_phase(raw_tail(c))
        _, im = cnum_is_real(c) ? (0, 0) : raw_realimag(expr)
        return Num(normalize_phase(im))
    end
    scalar, angle, flip_sine = pure_phase_parts(c)
    sine = flip_sine ? -real(scalar) : real(scalar)
    return scaled_trig(sine, sin, angle) + scaled_trig(imag(scalar), cos, angle)
end

function Base.abs(c::Coeff)::Num
    is_native(c) && return num_from_scalar(abs(native_scalar(c)))
    is_pure_phase(c) || throw(MethodError(abs, (c,)))
    return NUM_ONE
end

function Base.abs2(c::Coeff)::Num
    is_native(c) && return num_from_scalar(abs2(native_scalar(c)))
    is_pure_phase(c) || throw(MethodError(abs2, (c,)))
    return NUM_ONE
end

@inline function realimag(c::Coeff)
    is_native(c) && return (num_from_scalar(real(native_scalar(c))), num_from_scalar(imag(native_scalar(c))))
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
    t isa Native && return Complex(num_from_scalar(real(native_scalar(c))), num_from_scalar(imag(native_scalar(c))))
    t isa Poly && return poly_to_num(t)
    expr = normalize_phase(t.expr)
    re, im = cnum_is_real(c) ? (expr, 0) : raw_realimag(expr)
    return Complex(Num(normalize_phase(re)), Num(normalize_phase(im)))
end

Base.show(io::IO, c::Coeff) = show(io, to_num(c))

# Branch on the tail type so each `isequal` / `hash` call sees a concrete operand.
function Base.isequal(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Native
        return tb isa Native && isequal(native_scalar(a), native_scalar(b))
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
    t isa Native && return hash(native_scalar(c), hash(:CoeffNative, h))
    t isa Poly && return hash(t, hash(:CoeffSym, h))
    return hash(t.real_slot, hash(t.expr, hash(:CoeffRaw, h)))
end

# Coefficients are routinely compared against plain numbers / `Complex{Num}`
# (e.g. `get_prefactor(x) == 2`, `q[key] == CNUM_ONE`); promote the number side.
Base.isequal(a::Coeff, b::Number) = isequal(a, to_cnum(b))
Base.isequal(a::Number, b::Coeff) = isequal(to_cnum(a), b)
Base.:(==)(a::Coeff, b::Number) = isequal(a, to_cnum(b))
Base.:(==)(a::Number, b::Coeff) = isequal(to_cnum(a), b)

# Without this a coefficient falls into the iterable branch of `broadcastable`, and every
# array-scalar broadcast (`A .* c`) fails on `length(::Coeff)`.
Base.broadcastable(c::Coeff) = Ref(c)
# A coefficient is immutable. `matmul2x2!` on Julia 1.10 copies its operands.
Base.copy(c::Coeff) = c

Base.iszero(c::Coeff) = iszero_cnum(c)
Base.isone(c::Coeff) = isequal(c, CNUM_ONE)
Base.conj(c::Coeff) = conj_cnum(c)
Base.adjoint(c::Coeff) = conj_cnum(c)
Base.transpose(c::Coeff) = c
Base.zero(::Type{Coeff}) = CNUM_ZERO
Base.one(::Type{Coeff}) = CNUM_ONE
Base.zero(::Coeff) = CNUM_ZERO
Base.one(::Coeff) = CNUM_ONE
Base.oneunit(::Type{Coeff}) = CNUM_ONE
Base.oneunit(::Coeff) = CNUM_ONE

function native_div(a::NativeScalar, b::NativeScalar)::Coeff
    if a isa ExactComplex && b isa ExactComplex && !iszero(b)
        return tiered(exact_quotient, a, b)
    end
    return native(native_float(a) / native_float(b))
end
@inline exact_quotient(::Type{E}, a::E, b::E) where {E} = a * inv(b)

function inverse_terms(::Type{E}, terms::Vector{Monomial{E}}) where {E}
    monomial = only(terms)
    return monomial_terms(E, inv(monomial.scalar), monomial.syms, -monomial.exps)
end
@inline divided_terms(::Type{E}, terms::Vector{Monomial{E}}, z::CoeffScalar) where {E} =
    poly_scale(terms, inv(as_tier(E, z)))

function Base.inv(c::Coeff)::Coeff
    tail = c.tail
    tail isa Native && return native_div(EXACT_ONE, native_scalar(c))
    if tail isa Poly
        length(tail.terms) == 1 && return tiered(inverse_terms, tail)
        return CNUM_ONE / c
    end
    return from_raw(inv(tail.expr); real_slot = tail.real_slot)
end

# Coefficients flow through downstream code (and tests) as numbers; support the
# usual scalar arithmetic, promoting any `Number` operand into a `Coeff` first.
Base.:-(c::Coeff) = neg_cnum(c)
Base.:+(a::Coeff, b::Coeff) = add_cnum(a, b)
Base.:+(a::Coeff, b::Number) = add_cnum(a, to_cnum(b))
Base.:+(a::Number, b::Coeff) = add_cnum(to_cnum(a), b)
Base.:-(a::Coeff, b::Coeff) = add_cnum(a, neg_cnum(b))
Base.:-(a::Coeff, b::Number) = add_cnum(a, neg_cnum(to_cnum(b)))
Base.:-(a::Number, b::Coeff) = add_cnum(to_cnum(a), neg_cnum(b))
Base.:*(a::Coeff, b::Coeff) = mul_cnum(a, b)
Base.:*(a::Coeff, b::Number) = mul_cnum(a, to_cnum(b))
Base.:*(a::Number, b::Coeff) = mul_cnum(to_cnum(a), b)
function Base.:/(a::Coeff, b::Coeff)::Coeff
    (is_native(a) && is_native(b)) && return native_div(native_scalar(a), native_scalar(b))
    ta, tb = a.tail, b.tail
    tb isa Native && ta isa Poly && return tiered(divided_terms, ta, native_scalar(b))
    if tb isa Poly && length(tb.terms) == 1
        inverse = inv(b)
        inverse_tail = inverse.tail
        inverse_tail isa Poly || return mul_cnum(a, inverse)
        if ta isa Native
            return tiered(poly_scale, inverse_tail, native_scalar(a))
        elseif ta isa Poly
            return tiered(poly_mul, ta, inverse_tail)
        end
        any(is_radical_atom, only(inverse_tail.terms).syms) && return mul_cnum(a, inverse)
    end
    return from_raw(
        raw_binary(/, a, b);
        real_slot = cnum_is_real(a) && cnum_is_real(b),
    )
end
Base.:/(a::Coeff, b::Number) = a / to_cnum(b)
Base.:/(a::Number, b::Coeff) = to_cnum(a) / b

# `conj(conj(x)) == x`, so unwrap an existing `conj(...)` rather than nesting a
# second one (which never folds and survives downstream).
is_conj_call(x) =
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === conj

function raw_conj(x)
    x isa Number && return conj(x)
    SymbolicUtils.symtype(x) <: Real && return x
    is_phase(x) && return expim_expanded(-only(SymbolicUtils.arguments(x)))
    is_conj_call(x) && return only(SymbolicUtils.arguments(x))
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === complex
            return raw_complex(raw_conj(args[1]), -raw_conj(args[2]))
        end
        op === (+) && return foldl(+, map(raw_conj, args))
        op === (*) && return foldl(*, map(raw_conj, args))
        if op === (/)
            return raw_conj(args[1]) / raw_conj(args[2])
        end
    end
    return conj(x)
end

# Conjugate an atom factor: real-symtype atoms are self-conjugate, an existing
# `conj(...)` unwraps (involution), else wrap in `conj(...)` (still an atom to
# `is_atom`).
@inline function conj_atom(s::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.symtype(s) <: Real && return s
    is_conj_call(s) && return SymbolicUtils.arguments(s)[1]
    return SymbolicUtils.unwrap(conj(s))
end

# Native conjugation of a Poly: conjugate each scalar and atom, re-sort the rekeyed
# factors by `objectid`, re-canonicalize. Avoids the `to_num` round-trip per term.
conj_poly(p::Poly) = tiered(conj_terms, p)

function conj_terms(::Type{E}, source::Vector{Monomial{E}}) where {E}
    terms = Vector{Monomial{E}}(undef, length(source))
    @inbounds for k in eachindex(source)
        m = source[k]
        n = length(m.syms)
        if n == 0
            terms[k] = conj_monomial(m, m.syms, m.exps)
            continue
        end
        nsyms = Vector{SymbolicUtils.BasicSymbolic}(undef, n)
        nexps = copy(m.exps)
        for i in 1:n
            s = m.syms[i]
            # Unimodular: conjugation flips the exponent and keeps the atom, which is what
            # lets `p * conj(p)` cancel instead of growing a second unrelated factor.
            if is_phase(s)
                nsyms[i] = s
                nexps[i] = -nexps[i]
            else
                nsyms[i] = conj_atom(s)
            end
        end
        sort_factors!(nsyms, nexps)
        terms[k] = conj_monomial(m, nsyms, nexps)
    end
    return canonical_terms!(terms)
end

@inline function conj_cnum(c::Coeff)
    t = c.tail
    t isa Native && return native(conj(native_scalar(c)))
    t isa Poly && return conj_poly(t)
    return from_raw(raw_conj(t.expr); real_slot = t.real_slot)
end

@inline function iszero_num(x::Num)
    v = SymbolicUtils.unwrap(x)
    v isa Number && return iszero(v)
    return isequal(x, NUM_ZERO)
end

@inline function iszero_cnum(c::Coeff)
    is_native(c) && return iszero(native_scalar(c))
    c.tail isa Poly && return false   # a canonical Poly never sums to zero
    value = const_value(c.tail.expr)
    return value isa Number && iszero(value)
end

# Structural `a == -b`, used to recognize exact cancellation without a CAS round-trip.
@inline function raw_is_negative_of(
        positive::SymbolicUtils.BasicSymbolic, negative::SymbolicUtils.BasicSymbolic,
    )
    positive_value = const_value(positive)
    negative_value = const_value(negative)
    positive_value isa Number && negative_value isa Number &&
        positive_value == -negative_value && return true
    SymbolicUtils.iscall(negative) || return false
    op = SymbolicUtils.operation(negative)
    args = SymbolicUtils.arguments(negative)
    if op === (*) && length(args) >= 2 && const_value(args[1]) == -1
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
            raw_is_negative_of(positive_args[i], args[i]) || return false
        end
        return true
    end
    if op === (/) && length(args) == 2 &&
            SymbolicUtils.iscall(positive) && SymbolicUtils.operation(positive) === (/)
        positive_args = SymbolicUtils.arguments(positive)
        length(positive_args) == 2 || return false
        return isequal(args[2], positive_args[2]) &&
            (
            raw_is_negative_of(positive_args[1], args[1]) ||
                raw_is_negative_of(args[1], positive_args[1])
        )
    end
    return false
end

@inline function isneg_cnum(a::RawSymbolicCoeff, b::RawSymbolicCoeff)
    return a.real_slot == b.real_slot &&
        (
        raw_is_negative_of(a.expr, b.expr) ||
            raw_is_negative_of(b.expr, a.expr)
    )
end

@inline native_float(z::NativeScalar) = z isa ComplexF64 ? z : ComplexF64(z)

@noinline wide_product(a::ExactComplex, b::ExactComplex)::Coeff =
    scalar_coeff(BigExactComplex(a) * BigExactComplex(b))
@noinline wide_sum(a::ExactComplex, b::ExactComplex)::Coeff =
    scalar_coeff(BigExactComplex(a) + BigExactComplex(b))

@inline function mul_native(a::NativeScalar, b::NativeScalar)::Coeff
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = mul_checked(a, b)
        return overflow ? wide_product(a, b) : native(z)
    end
    return native(native_float(a) * native_float(b))
end
@inline function add_native(a::NativeScalar, b::NativeScalar)::Coeff
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = add_checked(a, b)
        return overflow ? wide_sum(a, b) : native(z)
    end
    return native(native_float(a) + native_float(b))
end

@inline function mul_cnum(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return mul_native(native_scalar(a), native_scalar(b))
    return mul_cnum_slow(a, b)
end

# Native and polynomial fast paths first, then multiply the intact raw expressions.
@noinline function mul_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return tiered(poly_mul, ta, tb)
    elseif ta isa Poly && tb isa Native
        return tiered(poly_scale, ta, native_scalar(b))
    elseif tb isa Poly && ta isa Native
        return tiered(poly_scale, tb, native_scalar(a))
    end
    raw = if is_native(a) || is_native(b)
        raw_binary(*, a, b)
    else
        raw_product(raw_tail(a), raw_tail(b))
    end
    return from_raw_arithmetic(raw, cnum_is_real(a) && cnum_is_real(b))
end

function is_numeric_radical(x::RawExpression)::Bool
    SymbolicUtils.iscall(x) || return false
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if (op === sqrt || op === cbrt) && length(args) == 1
        arg = only(args)
        return arg isa RawExpression && is_exact_number(arg) && exact_number(arg) >= 0
    elseif op === (^) && length(args) == 2
        base, exponent = args[1], args[2]
        (base isa RawExpression && exponent isa RawExpression) || return false
        return is_exact_number(base) && const_value(exponent) isa Rational{Int} &&
            exact_number(base) >= 0
    end
    return false
end

function radical_factor(x::RawExpression)::Coeff
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    op === sqrt && return numeric_radical(exact_number(only(args)), 1 // 2, x)
    op === cbrt && return numeric_radical(exact_number(only(args)), 1 // 3, x)
    return numeric_radical(exact_number(args[1]), const_value(args[2])::Rational{Int}, x)
end

@inline is_polynomial_coeff(c::Coeff) = !(c.tail isa RawSymbolicCoeff)

function mul_polynomial(a::Coeff, b::Coeff)::Coeff
    ta, tb = a.tail, b.tail
    (ta isa Native && tb isa Native) && return mul_native(native_scalar(a), native_scalar(b))
    ta isa Native && return tiered(poly_scale, tb::Poly, native_scalar(a))
    tb isa Native && return tiered(poly_scale, ta::Poly, native_scalar(b))
    return tiered(poly_mul, ta::Poly, tb::Poly)
end

function split_radicals(x::RawExpression)::Tuple{Coeff, Vector{RawExpression}}
    product = SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === (*)
    factors = product ? SymbolicUtils.arguments(x) : (x,)
    radicals = CNUM_ONE
    rest = RawExpression[]
    for f in factors
        f isa RawExpression || return (CNUM_ONE, RawExpression[x])
        if is_numeric_radical(f)
            c = radical_factor(f)
            if is_polynomial_coeff(c)
                radicals = mul_polynomial(radicals, c)
                continue
            end
        end
        push!(rest, f)
    end
    return (radicals, rest)
end

function raw_product(x::RawExpression, y::RawExpression)::RawExpression
    rx, restx = split_radicals(x)
    ry, resty = split_radicals(y)
    (is_native(rx) && is_native(ry)) && return (x * y)::RawExpression
    result = raw_tail(mul_polynomial(rx, ry))
    for f in restx
        result = (result * f)::RawExpression
    end
    for f in resty
        result = (result * f)::RawExpression
    end
    return result
end

@inline negated_terms(terms::Vector{Monomial{E}}) where {E} =
    Monomial{E}[Monomial{E}(normalize_scalar(-m.scalar), m.syms, m.exps) for m in terms]

@inline function neg_cnum(a::Coeff)
    t = a.tail
    t isa Native && return native(-native_scalar(a))
    t isa Poly && return poly_coeff(Poly(negated_terms(t.terms)))
    raw = (-t.expr)::SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}
    return from_raw_arithmetic(raw, t.real_slot)
end

@inline function add_cnum(a::Coeff, b::Coeff)
    (is_native(a) && is_native(b)) && return add_native(native_scalar(a), native_scalar(b))
    # Skip add-by-zero: `recognize` folds every sum from `CNUM_ZERO`, so without this each
    # fold would splice a throwaway zero `Monomial` into the Poly and merge it away.
    is_native(a) && iszero(native_scalar(a)) && return b
    is_native(b) && iszero(native_scalar(b)) && return a
    # Raw symbolic addition deliberately avoids a CAS round-trip. Recover the common exact
    # cancellation case here, which is needed by scalar identities such as the off-diagonal
    # entries of a symbolic orthogonal matrix product.
    if a.tail isa RawSymbolicCoeff && b.tail isa RawSymbolicCoeff
        isneg_cnum(a.tail, b.tail) && return CNUM_ZERO
    end
    return add_cnum_slow(a, b)
end

# Polynomial addition is a native merge (no escalation): this is what makes the
# tier pay off on sum-heavy workloads, where the single-monomial design regressed.
@noinline function add_cnum_slow(a::Coeff, b::Coeff)
    ta, tb = a.tail, b.tail
    if ta isa Poly && tb isa Poly
        return tiered(poly_add, ta, tb)
    elseif ta isa Poly && tb isa Native
        return tiered(constant_added_terms, ta, native_scalar(b))
    elseif tb isa Poly && ta isa Native
        return tiered(constant_added_terms, tb, native_scalar(a))
    end
    raw = raw_binary(+, a, b)
    return from_raw_arithmetic(raw, cnum_is_real(a) && cnum_is_real(b))
end

@inline constant_added_terms(::Type{E}, terms::Vector{Monomial{E}}, z::CoeffScalar) where {E} =
    poly_add(terms, Monomial{E}[constant_term(E, z)])

@inline function pow_cnum_nonnegative(base::Coeff, n::Int)
    result = CNUM_ONE
    while n > 0
        if isodd(n)
            result = if is_native(result) && is_native(base)
                mul_native(native_scalar(result), native_scalar(base))
            else
                mul_cnum_slow(result, base)
            end
        end
        n >>= 1
        if n > 0
            base = is_native(base) ? mul_native(native_scalar(base), native_scalar(base)) : mul_cnum_slow(base, base)
        end
    end
    return result
end

function pow_cnum_integer(base::Coeff, n::Int)::Coeff
    n >= 0 && return pow_cnum_nonnegative(base, n)
    n == typemin(Int) && return to_cnum(to_num(base)^n)
    return CNUM_ONE / pow_cnum_nonnegative(base, -n)
end

Base.:^(base::Coeff, n::Integer)::Coeff = pow_cnum_integer(base, Int(n))

@inline function rewrite_complex_slots(x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    SymbolicUtils.iscall(x) || return nothing
    SymbolicUtils.operation(x) === complex || return nothing
    args = SymbolicUtils.arguments(x)
    length(args) == 2 || return nothing
    return args[1] + im * args[2]
end

const COMPLEX_SLOT_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(rewrite_complex_slots),
)

@inline lower_complex_slots(x) = COMPLEX_SLOT_REWRITER(x)
@inline derivative_expression(c::Coeff) = lower_complex_slots(raw_tail(c))

(D::Symbolics.Differential)(c::Coeff) =
    is_native(c) ? D(to_num(c)) : D(derivative_expression(c))

function Symbolics.derivative(c::Coeff, var; simplify = false, kwargs...)::Coeff
    is_native(c) && return CNUM_ZERO
    input = derivative_expression(c)
    differentiated = Symbolics.expand_derivatives(
        Symbolics.Differential(var)(input), simplify; kwargs...,
    )
    return to_cnum(differentiated)
end

@inline factor_cnum(s::SymbolicUtils.BasicSymbolic, e::Rational{Int}) =
    tiered(monomial_terms, EXACT_ONE, SymbolicUtils.BasicSymbolic[s], Rational{Int}[e])

@inline function euler_cos(argument)
    phase = phase_coeff(argument)
    return mul_cnum(CNUM_HALF, add_cnum(phase, conj_cnum(phase)))
end

@inline function euler_sin(argument)
    phase = phase_coeff(argument)
    difference = add_cnum(phase, neg_cnum(conj_cnum(phase)))
    return mul_cnum(mul_cnum(CNUM_NEG_IM, CNUM_HALF), difference)
end

function phase_from_exponential_argument(argument)::Union{Coeff, Nothing}
    exponent = recognize(argument)
    real_part, imaginary_part = realimag(exponent)
    iszero_num(real_part) || return nothing
    imaginary = SymbolicUtils.unwrap(imaginary_part)
    SymbolicUtils.symtype(imaginary) <: Real || return nothing
    return phase_coeff(imaginary_part)
end

function exponential_tree(x)::Coeff
    u = SymbolicUtils.unwrap(x)
    u isa Number && return to_cnum(u)
    u isa SymbolicUtils.BasicSymbolic || return to_cnum(u)
    SymbolicUtils.iscall(u) || return recognize(u)
    op = SymbolicUtils.operation(u)
    args = SymbolicUtils.arguments(u)
    if op === exp && length(args) == 1
        phase = phase_from_exponential_argument(only(args))
        phase === nothing || return phase
        return recognize(u)
    elseif op === cis && length(args) == 1
        return phase_coeff(only(args))
    elseif op === cos && length(args) == 1
        return euler_cos(only(args))
    elseif op === sin && length(args) == 1
        return euler_sin(only(args))
    elseif op === (+)
        result = CNUM_ZERO
        for arg in args
            result = add_cnum(result, exponential_tree(arg))
        end
        return result
    elseif op === (*)
        result = CNUM_ONE
        for arg in args
            result = mul_cnum(result, exponential_tree(arg))
        end
        return result
    elseif op === (/) && length(args) == 2
        return exponential_tree(args[1]) / exponential_tree(args[2])
    elseif op === (^) && length(args) == 2
        exponent = const_value(args[2])
        exponent isa Integer || return recognize(u)
        return pow_cnum_integer(exponential_tree(args[1]), Int(exponent))
    end
    return recognize(u)
end

function exponential_monomial(m::Monomial)
    result = to_cnum(m.scalar)
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        exponent = m.exps[i]
        factor = factor_cnum(symbol, exponent)
        if denominator(exponent) == 1 && SymbolicUtils.iscall(symbol)
            op = SymbolicUtils.operation(symbol)
            if op === cos || op === sin
                n = numerator(exponent)
                base = op === cos ?
                    euler_cos(only(SymbolicUtils.arguments(symbol))) :
                    euler_sin(only(SymbolicUtils.arguments(symbol)))
                factor = pow_cnum_integer(base, n)
            end
        end
        result = mul_cnum(result, factor)
    end
    return result
end

function exponential_cnum(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    if tail isa Poly
        result = CNUM_ZERO
        for monomial in tail.terms
            result = add_cnum(result, exponential_monomial(monomial))
        end
        return result
    end
    return exponential_tree(tail.expr)
end

@inline function phase_trigonometric(symbol::SymbolicUtils.BasicSymbolic, n::Int)
    argument = only(SymbolicUtils.arguments(symbol))
    angle = expand(Num(n) * Num(argument))
    sine_sign = CNUM_ONE
    if leading_sign(angle) < 0
        angle = expand(-angle)
        sine_sign = CNUM_NEG1
    end
    sine = mul_cnum(sine_sign, to_cnum(sin(angle)))
    return add_cnum(to_cnum(cos(angle)), mul_cnum(CNUM_IM, sine))
end

function trigonometric_monomial(m::Monomial)
    result = to_cnum(m.scalar)
    @inbounds for i in eachindex(m.syms)
        symbol = m.syms[i]
        exponent = m.exps[i]
        factor = if is_phase(symbol) && denominator(exponent) == 1
            phase_trigonometric(symbol, numerator(exponent))
        else
            factor_cnum(symbol, exponent)
        end
        result = mul_cnum(result, factor)
    end
    return result
end

@inline function rewrite_phase_to_trig(x)
    is_phase(x) || return nothing
    argument = only(SymbolicUtils.arguments(x))
    return cos(argument) + im * sin(argument)
end

const TRIGONOMETRIC_REWRITER = SymbolicUtils.Rewriters.PassThrough(
    SymbolicUtils.Rewriters.Postwalk(rewrite_phase_to_trig),
)

function trigonometric_cnum(c::Coeff)
    tail = c.tail
    tail isa Native && return c
    if tail isa RawSymbolicCoeff
        expanded = SymbolicUtils.expand(TRIGONOMETRIC_REWRITER(tail.expr))
        return from_raw(
            simplify_raw(expanded); normalize = false,
        )
    end
    result = CNUM_ZERO
    for monomial in tail.terms
        result = add_cnum(result, trigonometric_monomial(monomial))
    end
    return result.tail isa RawSymbolicCoeff ? from_raw(result.tail.expr) : result
end

"""
    exponential_form(x)

Rewrite algebraic occurrences of `cos(θ)` and `sin(θ)` in a coefficient or quantum
expression using the exact phase atom [`expim`](@ref). The conversion is explicit and does
not affect ordinary display or [`simplify`](@ref).

See also [`trigonometric_form`](@ref).
"""
exponential_form(x::Coeff) = exponential_cnum(x)
exponential_form(x::Num) = exponential_tree(x)
exponential_form(x::SymbolicUtils.BasicSymbolic) = exponential_tree(x)
exponential_form(x::Number) = x
exponential_form(x::Complex{Num}) = exponential_cnum(to_cnum(x))

"""
    PhaseTerm

One term in the finite phase decomposition of a [`Coeff`](@ref). It represents
`amplitude * expim(phase)`.
"""
struct PhaseTerm
    amplitude::Coeff
    phase::Num
end

@inline constant_phase_term(c::Coeff) = PhaseTerm(c, NUM_ZERO)

function monomial_phase_term(m::Monomial)
    phase_index = phase_factor_index(m.syms)
    for i in eachindex(m.syms)
        i == phase_index && continue
        contains_phase(m.syms[i]) && throw(
            ArgumentError("a phase occurs inside an operation with no finite phase decomposition"),
        )
    end
    phase_index == 0 && return constant_phase_term(from_poly([m]))

    exponent = m.exps[phase_index]
    denominator(exponent) == 1 ||
        throw(ArgumentError("a unit phase cannot have a fractional power"))
    argument = only(SymbolicUtils.arguments(m.syms[phase_index]))
    phase = expand(Num(numerator(exponent)) * Num(argument))

    syms = copy(m.syms)
    exps = copy(m.exps)
    deleteat!(syms, phase_index)
    deleteat!(exps, phase_index)
    amplitude = from_poly([with_factors(m, syms, exps)])
    return PhaseTerm(amplitude, phase)
end

function contains_phase(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    is_phase(x) && return true
    SymbolicUtils.iscall(x) || return false
    return any(contains_phase, SymbolicUtils.arguments(x))
end

function multiply_phase_terms(left::Vector{PhaseTerm}, right::Vector{PhaseTerm})
    result = PhaseTerm[]
    sizehint!(result, length(left) * length(right))
    for a in left, b in right
        amplitude = mul_cnum(a.amplitude, b.amplitude)
        iszero(amplitude) && continue
        push!(result, PhaseTerm(amplitude, expand(a.phase + b.phase)))
    end
    return result
end

function phase_power_terms(base, exponent::Int)
    exponent < 0 && contains_phase(base) && throw(
        ArgumentError("a phase-bearing denominator does not have a finite phase decomposition"),
    )
    exponent < 0 && return PhaseTerm[constant_phase_term(pow_cnum_integer(recognize(base), exponent))]

    result = PhaseTerm[constant_phase_term(CNUM_ONE)]
    factor = raw_phase_terms(base)
    for _ in 1:exponent
        result = multiply_phase_terms(result, factor)
    end
    return result
end

function raw_phase_terms(x)::Vector{PhaseTerm}
    value = const_value(x)
    value isa Number && return PhaseTerm[constant_phase_term(to_cnum(value))]
    x isa SymbolicUtils.BasicSymbolic ||
        return PhaseTerm[constant_phase_term(to_cnum(x))]

    phase = phase_power(x)
    if phase !== nothing
        argument, exponent = phase
        return PhaseTerm[PhaseTerm(CNUM_ONE, expand(Num(exponent) * Num(argument)))]
    end

    SymbolicUtils.iscall(x) || return PhaseTerm[constant_phase_term(recognize(x))]
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op === (+)
        result = PhaseTerm[]
        for argument in args
            append!(result, raw_phase_terms(argument))
        end
        return result
    elseif op === (*)
        result = PhaseTerm[constant_phase_term(CNUM_ONE)]
        for factor in args
            result = multiply_phase_terms(result, raw_phase_terms(factor))
        end
        return result
    elseif op === (/) && length(args) == 2
        numerator, denominator_ = args
        contains_phase(denominator_) && throw(
            ArgumentError(
                "a phase-bearing denominator does not have a finite phase decomposition",
            ),
        )
        denominator_coeff = recognize(denominator_)
        return PhaseTerm[
            PhaseTerm(term.amplitude / denominator_coeff, term.phase) for
                term in raw_phase_terms(numerator)
        ]
    elseif op === (^) && length(args) == 2
        exponent = const_value(args[2])
        exponent isa Integer || begin
            contains_phase(args[1]) && throw(
                ArgumentError("a unit phase cannot have a non-integer power"),
            )
            return PhaseTerm[constant_phase_term(recognize(x))]
        end
        return phase_power_terms(args[1], Int(exponent))
    elseif op === conj && length(args) == 1
        return PhaseTerm[
            PhaseTerm(conj(term.amplitude), -term.phase) for
                term in raw_phase_terms(only(args))
        ]
    end

    contains_phase(x) && throw(
        ArgumentError("a phase occurs inside an operation with no finite phase decomposition"),
    )
    return PhaseTerm[constant_phase_term(recognize(x))]
end

"""
    phase_terms(c::Coeff) -> Vector{PhaseTerm}

Decompose a coefficient into a finite sum of `amplitude * expim(phase)` terms. Trigonometric
factors are converted to exponential form first. A phase-bearing denominator or a phase
inside an unsupported nonlinear operation throws an `ArgumentError` because it does not
represent a finite phase polynomial.

The decomposition is exact after conversion to exponential form:

```julia
exponential_form(c) == sum(term.amplitude * expim(term.phase) for term in phase_terms(c))
```
"""
function phase_terms(c::Coeff)::Vector{PhaseTerm}
    normalized = exponential_cnum(c)
    tail = normalized.tail
    if tail isa Native
        return iszero(normalized) ? PhaseTerm[] : PhaseTerm[constant_phase_term(normalized)]
    elseif tail isa Poly
        return PhaseTerm[monomial_phase_term(monomial) for monomial in tail.terms]
    end
    return raw_phase_terms(tail.expr)
end

"""
    trigonometric_form(x)

Rewrite integer powers of [`expim(θ)`](@ref expim) in a coefficient or quantum expression
as `cos(nθ) + im*sin(nθ)`. Opposite phases then combine through ordinary coefficient
arithmetic. The conversion is explicit and leaves the stored input unchanged.

See also [`exponential_form`](@ref).
"""
trigonometric_form(x::Coeff) = trigonometric_cnum(x)
trigonometric_form(x::Num) = trigonometric_cnum(to_cnum(x))
trigonometric_form(x::SymbolicUtils.BasicSymbolic) = trigonometric_cnum(to_cnum(x))
trigonometric_form(x::Number) = x
trigonometric_form(x::Complex{Num}) = trigonometric_cnum(to_cnum(x))

sort_key(op::QSym) = (op.space_index, name_rank(op.index.name_id))
