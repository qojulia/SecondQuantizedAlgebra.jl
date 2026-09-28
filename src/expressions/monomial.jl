"""
    Monomial

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`). A radical atom (a hashconsed `Const{SymReal}`
of a prime `Int`) always keeps its exponent in `(0, 1)`: any integer part is
folded into `scalar` (`canonical_monomial`, `radical_power`), so `√2·√2`
normalizes to the scalar `2`, and `√2·√3`/`√6`/`√12/2` all reduce to the same
two prime atoms `Const(2)^(1/2)`, `Const(3)^(1/2)` and compare `isequal`.
"""
# Exact Gaussian rationals live alongside the floating-point fast path. The small tier
# covers the literals produced by Julia's `//`; intermediate growth that overflows it
# promotes to the big tier, and every big result that fits is demoted again, so each
# exact value has exactly one stored representation (small XOR big, never both).
const ExactComplex = Complex{Rational{Int}}
const BigExactComplex = Complex{Rational{BigInt}}
const CoeffScalar = Union{ComplexF64, ExactComplex, BigExactComplex}
const ExactScalar = Union{ExactComplex, BigExactComplex}
const SmallScalar = Union{ComplexF64, ExactComplex}

# Largest integer below which every integer is a Float64: native integer values up to it
# are exact, larger integer-valued floats are not.
const MAX_EXACT_FLOAT = maxintfloat(Float64)

# The big scalar lives in its own field: a `BigInt` member would turn the small isbits union
# into a boxed pointer and allocate on every monomial of the fast path.
struct Monomial
    small::Union{ComplexF64, ExactComplex}
    wide::Union{Nothing, BigExactComplex}
    syms::Vector{SymbolicUtils.BasicSymbolic}   # sorted by objectid, distinct
    exps::Vector{Rational{Int}}                 # matching nonzero exponents
end

@inline Monomial(z::SmallScalar, syms, exps) = Monomial(z, nothing, syms, exps)
function Monomial(z::BigExactComplex, syms, exps)
    small = canonical_exact(z)
    small isa ExactComplex && return Monomial(small, nothing, syms, exps)
    return Monomial(zero(ComplexF64), small, syms, exps)
end

@inline function term_scalar(m::Monomial)::CoeffScalar
    wide = m.wide
    return wide === nothing ? m.small : wide
end

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactComplex) = z
@inline normalize_scalar(z::BigExactComplex) = canonical_exact(z)

@inline fits_int(x::Rational{BigInt}) =
    -typemax(Int) <= numerator(x) <= typemax(Int) && denominator(x) <= typemax(Int)

# Demote a big exact value to the small tier whenever it fits (I3).
function canonical_exact(z::BigExactComplex)::ExactScalar
    re, im = real(z), imag(z)
    (fits_int(re) && fits_int(im)) || return z
    return ExactComplex(
        Rational{Int}(Int(numerator(re)), Int(denominator(re))),
        Rational{Int}(Int(numerator(im)), Int(denominator(im))),
    )
end

@inline widen_exact(z::ExactComplex) = BigExactComplex(z)
@inline widen_exact(z::BigExactComplex) = z

# Small-tier arithmetic reports overflow through a flag instead of an `OverflowError`, so
# the common case pays for no exception handler; an overflowing operation is redone in the
# big tier. A small value never holds a `typemin(Int)` numerator (the big tier keeps it),
# so negation is always safe and a flagged overflow is the only way out of the tier.
@inline function coprime_parts(a::Int, b::Int)
    g = gcd(a, b)
    return (div(a, g), div(b, g))
end

@inline function small_rational(n::Int, d::Int, overflow::Bool)
    return (Base.unsafe_rational(n, d), overflow | (n == typemin(Int)))
end

@inline function checked_rational_mul(x::Rational{Int}, y::Rational{Int})
    xn, yd = coprime_parts(numerator(x), denominator(y))
    xd, yn = coprime_parts(denominator(x), numerator(y))
    n, o1 = Base.Checked.mul_with_overflow(xn, yn)
    d, o2 = Base.Checked.mul_with_overflow(xd, yd)
    return small_rational(n, d, o1 | o2)
end

@inline function checked_rational_add(x::Rational{Int}, y::Rational{Int})
    xd, yd = coprime_parts(denominator(x), denominator(y))
    a, o1 = Base.Checked.mul_with_overflow(numerator(x), yd)
    b, o2 = Base.Checked.mul_with_overflow(numerator(y), xd)
    n, o3 = Base.Checked.add_with_overflow(a, b)
    d, o4 = Base.Checked.mul_with_overflow(denominator(x), yd)
    (o1 | o2 | o3 | o4) && return (zero(Rational{Int}), true)
    g = gcd(n, d)
    return small_rational(div(n, g), div(d, g), false)
end

@inline function checked_exact_mul(a::ExactComplex, b::ExactComplex)
    ar, ai, br, bi = real(a), imag(a), real(b), imag(b)
    rr, o1 = checked_rational_mul(ar, br)
    ii, o2 = checked_rational_mul(ai, bi)
    ri, o3 = checked_rational_mul(ar, bi)
    ir, o4 = checked_rational_mul(ai, br)
    (o1 | o2 | o3 | o4) && return (a, true)
    re, o5 = checked_rational_add(rr, -ii)
    im, o6 = checked_rational_add(ri, ir)
    return (ExactComplex(re, im), o5 | o6)
end

@inline function checked_exact_add(a::ExactComplex, b::ExactComplex)
    re, o1 = checked_rational_add(real(a), real(b))
    im, o2 = checked_rational_add(imag(a), imag(b))
    return (ExactComplex(re, im), o1 | o2)
end

@inline function exact_mul(a::ExactScalar, b::ExactScalar)::ExactScalar
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = checked_exact_mul(a, b)
        overflow || return z
    end
    return canonical_exact(widen_exact(a) * widen_exact(b))
end

@inline function exact_add(a::ExactScalar, b::ExactScalar)::ExactScalar
    if a isa ExactComplex && b isa ExactComplex
        z, overflow = checked_exact_add(a, b)
        overflow || return z
    end
    return canonical_exact(widen_exact(a) + widen_exact(b))
end

function exact_inv(a::ExactComplex)::ExactScalar
    z = try
        inv(a)
    catch err
        err isa OverflowError || rethrow()
        return canonical_exact(inv(widen_exact(a)))
    end
    (numerator(real(z)) == typemin(Int) || numerator(imag(z)) == typemin(Int)) &&
        return widen_exact(z)
    return z
end
exact_inv(a::BigExactComplex)::ExactScalar = canonical_exact(inv(a))

# `z` raised to an integer power by repeated squaring in the exact tier: an overflow of the
# small tier promotes the running product rather than throwing (I5) — the same mechanism
# `radical_power` below uses to fold a radical's integer exponent into the scalar exactly.
function exact_pow(z::ExactScalar, n::Int)::ExactScalar
    e = abs(n)
    result::ExactScalar = ExactComplex(1 // 1, 0 // 1)
    base::ExactScalar = z
    while e > 0
        isodd(e) && (result = exact_mul(result, base))
        e >>= 1
        e > 0 && (base = exact_mul(base, base))
    end
    return n >= 0 ? result : exact_inv(result)
end

scalar_conj(a::ExactComplex) = conj(a)
scalar_conj(a::BigExactComplex)::ExactScalar = canonical_exact(conj(a))
scalar_conj(a::ComplexF64) = conj(a)

@inline to_float_scalar(z::ComplexF64) = z
@inline to_float_scalar(z::ExactScalar) = ComplexF64(Float64(real(z)), Float64(imag(z)))

# Native integer/Gaussian-integer factors (1, -1, and `im`) do not introduce
# inexactness when multiplied by an exact scalar. Genuine non-integral floats do, and
# so do integer-valued floats beyond `2^53`, which no longer denote a unique integer (I4/I7).
@inline function integer_scalar(z::ComplexF64)
    re, im = real(z), imag(z)
    (abs(re) <= MAX_EXACT_FLOAT && abs(im) <= MAX_EXACT_FLOAT) || return nothing
    (isinteger(re) && isinteger(im)) || return nothing
    return ExactComplex(Int(re) // 1, Int(im) // 1)
end

# Whether the Float64 product of two native values is exact whenever both are integers:
# every partial product and partial sum then stays within `2^53`.
@inline native_product_exact(a::ComplexF64, b::ComplexF64) =
    (abs(real(a)) + abs(imag(a))) * (abs(real(b)) + abs(imag(b))) <= MAX_EXACT_FLOAT
# A rounded integer sum lands at or beyond `2^53`, so a smaller result is exact.
@inline native_sum_exact(s::ComplexF64) =
    abs(real(s)) < MAX_EXACT_FLOAT && abs(imag(s)) < MAX_EXACT_FLOAT

@noinline function native_mul_wide(a::ComplexF64, b::ComplexF64)::CoeffScalar
    ea, eb = integer_scalar(a), integer_scalar(b)
    (ea === nothing || eb === nothing) && return normalize_scalar(a * b)
    return exact_mul(ea, eb)
end

@noinline function native_add_wide(a::ComplexF64, b::ComplexF64)::CoeffScalar
    ea, eb = integer_scalar(a), integer_scalar(b)
    (ea === nothing || eb === nothing) && return normalize_scalar(a + b)
    return exact_add(ea, eb)
end

@inline function scalar_inv(z::ComplexF64)
    exact = integer_scalar(z)
    return exact === nothing ? inv(z) : exact_inv(exact)
end
@inline scalar_inv(z::ExactScalar) = exact_inv(z)

# General scalar arithmetic for the paths off the float fast path. One unspecialized entry
# per operation branches on the operand forms by hand: a call split over the nine
# combinations of two `CoeffScalar`s exceeds the union-split limit and would dispatch at
# run time, while this signature is invoked statically from a union-typed call site.
function scalar_mul(@nospecialize(a::CoeffScalar), @nospecialize(b::CoeffScalar))::CoeffScalar
    if a isa ComplexF64
        if b isa ComplexF64
            native_product_exact(a, b) && return normalize_scalar(a * b)
            return native_mul_wide(a, b)
        end
        return exact_float_mul(b, a)
    end
    b isa ComplexF64 && return exact_float_mul(a, b)
    return exact_mul(a, b)
end

@inline function exact_float_mul(a::ExactScalar, b::ComplexF64)::CoeffScalar
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(to_float_scalar(a) * b) : exact_mul(a, ib)
end

function scalar_add(@nospecialize(a::CoeffScalar), @nospecialize(b::CoeffScalar))::CoeffScalar
    if a isa ComplexF64
        if b isa ComplexF64
            s = a + b
            native_sum_exact(s) && return normalize_scalar(s)
            return native_add_wide(a, b)
        end
        return exact_float_add(b, a)
    end
    b isa ComplexF64 && return exact_float_add(a, b)
    return exact_add(a, b)
end

@inline function exact_float_add(a::ExactScalar, b::ComplexF64)::CoeffScalar
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(to_float_scalar(a) + b) : exact_add(a, ib)
end

# Hot-path scalar arithmetic builds its result straight into a monomial. A local holding
# the whole `CoeffScalar` union is boxed, since the union mixes isbits and `BigInt`
# members, so these read the isbits `small` field and leave the big tier to a fallback.
#
# Precondition shared by every helper below except `canonical_monomial`: `syms`/`exps` are
# already radical-canonical (I2 holds for every factor). Only `canonical_monomial` may see
# an arbitrary, not-yet-canonical triple.
@inline with_factors(m::Monomial, syms, exps) = Monomial(m.small, m.wide, syms, exps)
@inline scalar_iszero(m::Monomial) = m.wide === nothing && iszero(m.small)
@inline function scalar_is_real(m::Monomial)
    wide = m.wide
    return wide === nothing ? iszero(imag(m.small)) : iszero(imag(wide))
end
@inline function scalar_isequal(a::Monomial, b::Monomial)
    (a.wide === nothing && b.wide === nothing) && return isequal(a.small, b.small)
    return isequal(term_scalar(a), term_scalar(b))
end
@inline function hash_scalar(m::Monomial, h::UInt)
    wide = m.wide
    return wide === nothing ? hash(m.small, h) : hash(wide, h)
end
@inline function normalize_monomial(m::Monomial)
    small = m.small
    (m.wide === nothing && small isa ComplexF64) || return m
    return Monomial(normalize_scalar(small), nothing, m.syms, m.exps)
end
@inline function conj_monomial(m::Monomial, syms, exps)
    small = m.small
    (m.wide === nothing && small isa ComplexF64) &&
        return Monomial(conj(small), nothing, syms, exps)
    return Monomial(scalar_conj(term_scalar(m)), syms, exps)
end

@inline exact_operand(z::ComplexF64) = integer_scalar(z)
@inline exact_operand(z::ExactComplex) = z

@inline function exact_mul_monomial(x::ExactComplex, y::ExactComplex, syms, exps)::Monomial
    z, overflow = checked_exact_mul(x, y)
    overflow || return Monomial(z, nothing, syms, exps)
    return Monomial(widen_exact(x) * widen_exact(y), syms, exps)
end

@inline function exact_add_monomial(x::ExactComplex, y::ExactComplex, syms, exps)::Monomial
    z, overflow = checked_exact_add(x, y)
    overflow || return Monomial(z, nothing, syms, exps)
    return Monomial(widen_exact(x) + widen_exact(y), syms, exps)
end

@inline function mul_small(x::SmallScalar, y::SmallScalar, syms, exps)::Monomial
    (x isa ComplexF64 && y isa ComplexF64) && return Monomial(native_mul_wide(x, y), syms, exps)
    ex, ey = exact_operand(x), exact_operand(y)
    if ex === nothing || ey === nothing
        z = normalize_scalar(to_float_scalar(x) * to_float_scalar(y))
        return Monomial(z, nothing, syms, exps)
    end
    return exact_mul_monomial(ex, ey, syms, exps)
end

@inline function add_small(x::SmallScalar, y::SmallScalar, syms, exps)::Monomial
    (x isa ComplexF64 && y isa ComplexF64) && return Monomial(native_add_wide(x, y), syms, exps)
    ex, ey = exact_operand(x), exact_operand(y)
    if ex === nothing || ey === nothing
        z = normalize_scalar(to_float_scalar(x) + to_float_scalar(y))
        return Monomial(z, nothing, syms, exps)
    end
    return exact_add_monomial(ex, ey, syms, exps)
end

# Product of the two scalars carrying the given factors. Precondition: `syms`/`exps`
# already satisfy I2 (this does not fold radicals — use `canonical_monomial` when they
# might not, e.g. after merging two operands' factor lists).
@inline function mul_scalars(a::Monomial, b::Monomial, syms, exps)::Monomial
    x, y = a.small, b.small
    if a.wide === nothing && b.wide === nothing &&
            x isa ComplexF64 && y isa ComplexF64 && native_product_exact(x, y)
        return Monomial(normalize_scalar(x * y), nothing, syms, exps)
    end
    return mul_scalars_slow(a, b, syms, exps)
end
@noinline function mul_scalars_slow(a::Monomial, b::Monomial, syms, exps)::Monomial
    (a.wide === nothing && b.wide === nothing) && return mul_small(a.small, b.small, syms, exps)
    return Monomial(scalar_mul(term_scalar(a), term_scalar(b)), syms, exps)
end

# Sum of the scalars of two like terms, keeping the factors of `a`. Precondition: same as
# `mul_scalars`.
@inline function add_scalars(a::Monomial, b::Monomial)::Monomial
    x, y = a.small, b.small
    if a.wide === nothing && b.wide === nothing && x isa ComplexF64 && y isa ComplexF64
        s = x + y
        native_sum_exact(s) && return Monomial(normalize_scalar(s), nothing, a.syms, a.exps)
    end
    return add_scalars_slow(a, b)
end
@noinline function add_scalars_slow(a::Monomial, b::Monomial)::Monomial
    (a.wide === nothing && b.wide === nothing) && return add_small(a.small, b.small, a.syms, a.exps)
    return Monomial(scalar_add(term_scalar(a), term_scalar(b)), a.syms, a.exps)
end

@inline function scale_monomial(t::Monomial, z::SmallScalar)::Monomial
    x = t.small
    if t.wide === nothing && x isa ComplexF64 && z isa ComplexF64 && native_product_exact(x, z)
        return Monomial(normalize_scalar(x * z), nothing, t.syms, t.exps)
    end
    return scale_monomial_slow(t, z)
end
@noinline function scale_monomial_slow(t::Monomial, z::CoeffScalar)::Monomial
    (t.wide === nothing && !(z isa BigExactComplex)) && return mul_small(t.small, z, t.syms, t.exps)
    return Monomial(scalar_mul(term_scalar(t), z), t.syms, t.exps)
end

# === Radical folding (I2, generalized to compose with the big tier per I5) ===

# A numeric radical is a factor whose atom is the hashconsed `Const` of a *prime*, with an
# exponent in `(0, 1)`. Integer parts live in the exact scalar, so `√2·√2 = 2` and
# `1/√2 = √2/2` each have one representation. Factoring always reduces to prime bases
# (`radical_coeff`/`prime_factorization!`), so `√2·√3` and `√6` reach the *same* two-atom
# representation `Const(2)^(1/2) * Const(3)^(1/2)` — neither is ever combined into a single
# `Const(6)` atom.
@inline is_radical_atom(s::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.isconst(s) && s.val isa Int

# `p^n` folded into the scalar exactly (I5): repeated squaring over the exact tier, so a
# large exponent widens to `BigExactComplex` rather than falling back to `Float64`.
radical_power(p::Int, n::Int)::ExactScalar = exact_pow(ExactComplex(p // 1, 0 // 1), n)

@inline function needs_radical_fold(
        syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}},
    )
    @inbounds for i in eachindex(syms)
        is_radical_atom(syms[i]) || continue
        e = exps[i]
        (e <= 0 || e >= 1) && return true
    end
    return false
end

"""
    canonical_monomial(scalar, syms, exps) -> Monomial

The one function allowed to receive an arbitrary, not-yet-canonical
`(scalar, syms, exps)` triple. Folds every radical exponent into `(0, 1)`
(I2), routing the extracted integer power through the exact scalar tier so it
widens rather than overflowing or falling back to a float (I5), and then
dispatches to the tier-aware inner `Monomial` constructor on the resulting
scalar (I3). Every other scalar-only helper above requires `syms`/`exps` to
already be radical-canonical.
"""
function canonical_monomial(
        scalar::CoeffScalar,
        syms::Vector{SymbolicUtils.BasicSymbolic},
        exps::Vector{Rational{Int}},
    )::Monomial
    needs_radical_fold(syms, exps) || return Monomial(scalar, syms, exps)
    osyms = SymbolicUtils.BasicSymbolic[]
    oexps = Rational{Int}[]
    sizehint!(osyms, length(syms)); sizehint!(oexps, length(exps))
    @inbounds for i in eachindex(syms)
        s = syms[i]
        e = exps[i]
        if is_radical_atom(s)
            q = fld(numerator(e), denominator(e))
            if q != 0
                scalar = scalar_mul(scalar, radical_power(s.val::Int, q))
                e -= q
            end
            iszero(e) && continue
        end
        push!(osyms, s); push!(oexps, e)
    end
    return Monomial(scalar, osyms, oexps)
end

"""
    Poly

A sparse multivariate polynomial over named parameters (a sum of distinct
[`Monomial`](@ref) terms in canonical order), kept off SymbolicUtils hashconsing
and lowered to `Complex{Num}` only at the symbolic boundaries (see `poly_to_num`).
"""
struct Poly
    terms::Vector{Monomial}
end

# Factor identity key: SymbolicUtils hashconses, so `objectid`/`===` are exact and
# type-stable factor identity (unlike `hash`/`isequal` on abstract `BasicSymbolic`).
@inline fkey(s::SymbolicUtils.BasicSymbolic) = objectid(s)

@inline function same_factors(a::Monomial, b::Monomial)
    length(a.syms) == length(b.syms) || return false
    a.exps == b.exps || return false
    @inbounds for i in eachindex(a.syms)
        a.syms[i] === b.syms[i] || return false
    end
    return true
end

# Total order on monomial factor sets (factors are pre-sorted by objectid within
# each monomial), giving a canonical term order for `Poly` equality / hashing.
function term_less(a::Monomial, b::Monomial)
    la, lb = length(a.syms), length(b.syms)
    la != lb && return la < lb
    @inbounds for i in 1:la
        ha, hb = fkey(a.syms[i]), fkey(b.syms[i])
        ha != hb && return ha < hb
        a.exps[i] != b.exps[i] && return a.exps[i] < b.exps[i]
    end
    return false
end

# Merge two sorted factor lists, summing exponents and dropping cancellations.
function merge_factors(syma, expa, symb, expb)
    ia, ib = 1, 1
    na, nb = length(syma), length(symb)
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    sizehint!(syms, na + nb); sizehint!(exps, na + nb)
    @inbounds while ia <= na || ib <= nb
        if ib > nb || (ia <= na && fkey(syma[ia]) < fkey(symb[ib]))
            push!(syms, syma[ia]); push!(exps, expa[ia]); ia += 1
        elseif ia > na || fkey(syma[ia]) > fkey(symb[ib])
            push!(syms, symb[ib]); push!(exps, expb[ib]); ib += 1
        else
            e = expa[ia] + expb[ib]
            e != 0 && (push!(syms, syma[ia]); push!(exps, e))
            ia += 1; ib += 1
        end
    end
    return (syms, exps)
end

# Two phase-bearing factors that do not simply cancel: the phase arguments combine
# symbolically. Kept out of `term_mul` so its common path holds no boxed scalar.
@noinline function phase_term_mul(a::Monomial, b::Monomial, phase_a::Int, phase_b::Int)
    scalar = scalar_mul(term_scalar(a), term_scalar(b))
    if length(a.syms) == 1 && length(b.syms) == 1 &&
            a.syms[phase_a] === b.syms[phase_b]
        exponent = a.exps[phase_a] + b.exps[phase_b]
        return scaled_phase_monomial(scalar, Num(a.syms[phase_a]), exponent)
    end
    if length(a.syms) == 1 && length(b.syms) == 1
        return merged_phase_monomial(
            scalar,
            Num(a.syms[phase_a]),
            a.exps[phase_a],
            Num(b.syms[phase_b]),
            b.exps[phase_b],
        )
    end
    syms = vcat(a.syms, b.syms)
    exps = vcat(a.exps, b.exps)
    return canonical_phase_monomial(scalar, syms, exps)
end

function term_mul(a::Monomial, b::Monomial)
    isempty(a.syms) && return mul_scalars(a, b, b.syms, b.exps)
    isempty(b.syms) && return mul_scalars(a, b, a.syms, a.exps)
    phase_a = phase_factor_index(a.syms)
    phase_b = phase_factor_index(b.syms)
    if phase_a != 0 && phase_b != 0
        # The common inverse pair stays on the ordinary identity merge: no symbolic
        # argument arithmetic and no additional allocation. The merge can still leave a
        # shared radical atom's exponent outside `(0, 1)` (e.g. two factors each
        # contributing a fractional power of the same prime), so this still folds.
        if a.syms[phase_a] === b.syms[phase_b] &&
                a.exps[phase_a] == -b.exps[phase_b]
            se = merge_factors(a.syms, a.exps, b.syms, b.exps)
            return canonical_monomial(scalar_mul(term_scalar(a), term_scalar(b)), se[1], se[2])
        end
        return phase_term_mul(a, b, phase_a, phase_b)
    end
    se = merge_factors(a.syms, a.exps, b.syms, b.exps)
    return canonical_monomial(scalar_mul(term_scalar(a), term_scalar(b)), se[1], se[2])
end

# Insertion sort by a strict-less predicate. The polynomial passes sort very short
# vectors; this keeps the whole `Base.Sort` machinery (ScratchQuickSort, partition!,
# issorted, ...) out of inference, which is a large chunk of first-call latency.
function insertion_sort!(v::AbstractVector, lt::F) where {F}
    @inbounds for i in 2:length(v)
        x = v[i]
        j = i - 1
        while j >= 1 && lt(x, v[j])
            v[j + 1] = v[j]
            j -= 1
        end
        v[j + 1] = x
    end
    return v
end

# Sort terms into canonical order, merge like-factor terms, drop zero scalars.
function canonical_terms!(terms::Vector{Monomial})
    isempty(terms) && return terms
    insertion_sort!(terms, term_less)
    w = 0
    @inbounds for r in eachindex(terms)
        t = terms[r]
        if w > 0 && same_factors(terms[w], t)
            terms[w] = add_scalars(terms[w], t)
        else
            w += 1
            terms[w] = normalize_monomial(t)
        end
    end
    resize!(terms, w)
    filter!(t -> !scalar_iszero(t), terms)
    return terms
end

# Sorted merge of two canonical term lists, dropping zero-scalar terms so the result
# stays canonical and zero-free (a stray zero would break Poly equality/hashing).
function poly_add(p::Vector{Monomial}, q::Vector{Monomial})
    out = Monomial[]
    sizehint!(out, length(p) + length(q))
    ip, iq = 1, 1
    np, nq = length(p), length(q)
    @inbounds while ip <= np || iq <= nq
        if iq > nq || (ip <= np && term_less(p[ip], q[iq]))
            t = p[ip]; ip += 1
            scalar_iszero(t) || push!(out, t)
        elseif ip > np || term_less(q[iq], p[ip])
            t = q[iq]; iq += 1
            scalar_iszero(t) || push!(out, t)
        else   # same factor set: sum scalars, drop exact cancellations
            s = add_scalars(p[ip], q[iq])
            scalar_iszero(s) || push!(out, s)
            ip += 1; iq += 1
        end
    end
    return out
end

function poly_mul(p::Vector{Monomial}, q::Vector{Monomial})
    if length(p) == 1 && length(q) == 1
        return Monomial[normalize_monomial(term_mul(p[1], q[1]))]
    end
    out = Monomial[]
    sizehint!(out, length(p) * length(q))
    for a in p, b in q
        push!(out, term_mul(a, b))
    end
    return canonical_terms!(out)
end

# Scale every term; preserves canonical order (factors unchanged, so tier-only).
function poly_scale(p::Vector{Monomial}, z::SmallScalar)
    iszero(z) && return Monomial[]
    return Monomial[scale_monomial(t, z) for t in p]
end
poly_scale(p::Vector{Monomial}, z::BigExactComplex) =
    Monomial[scale_monomial_slow(t, z) for t in p]

function Base.isequal(a::Poly, b::Poly)
    length(a.terms) == length(b.terms) || return false
    @inbounds for i in eachindex(a.terms)
        ta, tb = a.terms[i], b.terms[i]
        (scalar_isequal(ta, tb) && same_factors(ta, tb)) || return false
    end
    return true
end
Base.:(==)(a::Poly, b::Poly) = isequal(a, b)
function Base.hash(p::Poly, h::UInt)
    @inbounds for t in p.terms
        h = hash_scalar(t, h)
        for i in eachindex(t.syms)
            h = hash(t.exps[i], hash(fkey(t.syms[i]), h))
        end
    end
    return hash(:Poly, h)
end
