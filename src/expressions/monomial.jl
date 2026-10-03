"""
    Monomial{E}

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`). An exact `scalar` has the type `E` of its tier,
`ExactComplex` or `BigExactComplex`, and every term of one `Poly` shares that tier.
A radical atom (a hashconsed `Const{SymReal}` of a prime `Int`) always keeps its
exponent in `(0, 1)`: the constructor folds any integer part into `scalar`, so
`√2·√2` normalizes to the scalar `2`, and `√2·√3`/`√6`/`√24/2` all reduce to the
same two prime atoms `Const(2)^(1/2)`, `Const(3)^(1/2)` and compare `isequal`.
"""
struct Monomial{E <: ExactScalar}
    scalar::Union{ComplexF64, E}
    syms::Vector{SymbolicUtils.BasicSymbolic}   # sorted by objectid, distinct
    exps::Vector{Rational{Int}}                 # matching nonzero exponents

    function Monomial{E}(
            scalar::Union{ComplexF64, E},
            syms::Vector{SymbolicUtils.BasicSymbolic},
            exps::Vector{Rational{Int}},
        ) where {E <: ExactScalar}
        needs_radical_fold(syms, exps) || return new{E}(scalar, syms, exps)
        osyms = SymbolicUtils.BasicSymbolic[]
        oexps = Rational{Int}[]
        sizehint!(osyms, length(syms)); sizehint!(oexps, length(exps))
        @inbounds for i in eachindex(syms)
            s = syms[i]
            e = exps[i]
            if is_radical_atom(s)
                q = fld(numerator(e), denominator(e))
                if q != 0
                    scalar = scalar_mul(E, scalar, radical_power(E, s.val::Int, q))
                    e -= q
                end
                iszero(e) && continue
            end
            push!(osyms, s); push!(oexps, e)
        end
        return new{E}(scalar, osyms, oexps)
    end
end

@inline as_tier(::Type{E}, m::Monomial{E}) where {E <: ExactScalar} = m
@inline as_tier(::Type{E}, m::Monomial) where {E <: ExactScalar} =
    Monomial{E}(as_tier(E, m.scalar), m.syms, m.exps)
@inline as_tier(::Type{E}, terms::Vector{Monomial{E}}) where {E <: ExactScalar} = terms
as_tier(::Type{E}, terms::Vector{<:Monomial}) where {E <: ExactScalar} =
    Monomial{E}[as_tier(E, m) for m in terms]

@inline with_factors(
    m::Monomial{E}, syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}},
) where {E} = Monomial{E}(m.scalar, syms, exps)
@inline scalar_iszero(m::Monomial) = iszero(m.scalar)
@inline scalar_is_real(m::Monomial) = iszero(imag(m.scalar))
@inline scalar_isequal(a::Monomial, b::Monomial) = isequal(a.scalar, b.scalar)
@inline hash_scalar(m::Monomial, h::UInt) = hash(m.scalar, h)
@inline function normalize_monomial(m::Monomial{E}) where {E}
    scalar = m.scalar
    scalar isa ComplexF64 || return m
    return Monomial{E}(normalize_scalar(scalar), m.syms, m.exps)
end
@inline conj_monomial(
    m::Monomial{E}, syms::Vector{SymbolicUtils.BasicSymbolic}, exps::Vector{Rational{Int}},
) where {E} =
    Monomial{E}(conj(m.scalar), syms, exps)

@inline mul_scalars(a::Monomial{E}, b::Monomial{E}, syms, exps) where {E} =
    Monomial{E}(scalar_mul(E, a.scalar, b.scalar), syms, exps)
@inline add_scalars(a::Monomial{E}, b::Monomial{E}) where {E} =
    Monomial{E}(scalar_add(E, a.scalar, b.scalar), a.syms, a.exps)
@inline scale_monomial(t::Monomial{E}, z::Union{ComplexF64, E}) where {E} =
    Monomial{E}(scalar_mul(E, t.scalar, z), t.syms, t.exps)

"""
    Poly

A sparse multivariate polynomial over named parameters (a sum of distinct
[`Monomial`](@ref) terms in canonical order), kept off SymbolicUtils hashconsing
and lowered to `Complex{Num}` only at the symbolic boundaries (see `poly_to_num`).
"""
struct Poly
    terms::Union{Vector{Monomial{ExactComplex}}, Vector{Monomial{BigExactComplex}}}
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

@noinline function phase_term_mul(
        a::Monomial{E}, b::Monomial{E}, phase_a::Int, phase_b::Int,
    )::Monomial{E} where {E}
    scalar = scalar_mul(E, a.scalar, b.scalar)
    if length(a.syms) == 1 && length(b.syms) == 1 &&
            a.syms[phase_a] === b.syms[phase_b]
        exponent = a.exps[phase_a] + b.exps[phase_b]
        return scaled_phase_monomial(E, scalar, Num(a.syms[phase_a]), exponent)
    end
    if length(a.syms) == 1 && length(b.syms) == 1
        return merged_phase_monomial(
            E,
            scalar,
            Num(a.syms[phase_a]),
            a.exps[phase_a],
            Num(b.syms[phase_b]),
            b.exps[phase_b],
        )
    end
    syms = vcat(a.syms, b.syms)
    exps = vcat(a.exps, b.exps)
    return canonical_phase_monomial(E, scalar, syms, exps)
end

function term_mul(a::Monomial{E}, b::Monomial{E})::Monomial{E} where {E}
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

@inline function merged_term_mul(a::Monomial{E}, b::Monomial{E})::Monomial{E} where {E}
    se = merge_factors(a.syms, a.exps, b.syms, b.exps)
    syms = se[1]::Vector{SymbolicUtils.BasicSymbolic}
    exps = se[2]::Vector{Rational{Int}}
    return mul_scalars(a, b, syms, exps)
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
function canonical_terms!(terms::Vector{Monomial{E}}) where {E}
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
function poly_add(p::Vector{Monomial{E}}, q::Vector{Monomial{E}}) where {E}
    out = Monomial{E}[]
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

function poly_mul(p::Vector{Monomial{E}}, q::Vector{Monomial{E}}) where {E}
    if length(p) == 1 && length(q) == 1
        product = normalize_monomial(term_mul(p[1], q[1]))
        return scalar_iszero(product) ? Monomial{E}[] : Monomial{E}[product]
    end
    out = Monomial{E}[]
    sizehint!(out, length(p) * length(q))
    for a in p, b in q
        push!(out, term_mul(a, b))
    end
    return canonical_terms!(out)
end

# Scale every term; preserves canonical order (factors unchanged).
function poly_scale(p::Vector{Monomial{E}}, z::Union{ComplexF64, E}) where {E}
    iszero(z) && return Monomial{E}[]
    return filter!(!scalar_iszero, Monomial{E}[scale_monomial(t, z) for t in p])
end

poly_add(::Type{E}, p::Vector{Monomial{E}}, q::Vector{Monomial{E}}) where {E} = poly_add(p, q)
poly_mul(::Type{E}, p::Vector{Monomial{E}}, q::Vector{Monomial{E}}) where {E} = poly_mul(p, q)
poly_scale(::Type{E}, p::Vector{Monomial{E}}, z) where {E} = poly_scale(p, as_tier(E, z))

Base.isequal(a::Poly, b::Poly) = terms_isequal(a.terms, b.terms)
Base.:(==)(a::Poly, b::Poly) = isequal(a, b)
function terms_isequal(a::Vector{Monomial{A}}, b::Vector{Monomial{B}}) where {A, B}
    A === B || return false
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        ta, tb = a[i], b[i]
        (scalar_isequal(ta, tb) && same_factors(ta, tb)) || return false
    end
    return true
end

Base.hash(p::Poly, h::UInt) = hash(:Poly, hash_terms(p.terms, h))
function hash_terms(terms::Vector{<:Monomial}, h::UInt)
    @inbounds for t in terms
        h = hash_scalar(t, h)
        for i in eachindex(t.syms)
            h = hash(t.exps[i], hash(fkey(t.syms[i]), h))
        end
    end
    return h
end
