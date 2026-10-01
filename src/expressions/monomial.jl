"""
    Monomial

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`). A radical atom (a hashconsed `Const{SymReal}`
of a prime `Int`) always keeps its exponent in `(0, 1)`: any integer part is
folded into `scalar` (`canonical_monomial`, `radical_power`), so `√2·√2`
normalizes to the scalar `2`, and `√2·√3`/`√6`/`√24/2` all reduce to the same
two prime atoms `Const(2)^(1/2)`, `Const(3)^(1/2)` and compare `isequal`.
"""

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
        # argument arithmetic and no additional allocation.
        if a.syms[phase_a] === b.syms[phase_b] &&
                a.exps[phase_a] == -b.exps[phase_b]
            return merged_term_mul(a, b)
        end
        return phase_term_mul(a, b, phase_a, phase_b)
    end
    return merged_term_mul(a, b)
end

@inline function merged_term_mul(a::Monomial, b::Monomial)
    se = merge_factors(a.syms, a.exps, b.syms, b.exps)
    syms = se[1]::Vector{SymbolicUtils.BasicSymbolic}
    exps = se[2]::Vector{Rational{Int}}
    needs_radical_fold(syms, exps) || return mul_scalars(a, b, syms, exps)
    return canonical_monomial(scalar_mul(term_scalar(a), term_scalar(b)), syms, exps)
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

# Scale every term; preserves canonical order (factors unchanged).
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
