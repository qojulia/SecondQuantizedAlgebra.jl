"""
    Monomial{E}

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`). An exact `scalar` has the type `E` of its tier,
`ExactComplex` or `BigExactComplex`, and every term of one `Poly` shares that tier.
The `slot` field stores the scalar as an `E`, with a float packed into it (`pack_scalar`),
so the layout is concrete; read it with `monomial_scalar`.
A radical atom (a hashconsed `Const{SymReal}` of a prime, see `radical_atom`) always keeps its
exponent in `(0, 1)`: the constructor folds any integer part into `scalar`, so
`√2·√2` normalizes to the scalar `2`, and `√2·√3`/`√6`/`√24/2` all reduce to the
same two prime atoms `Const(2)^(1/2)`, `Const(3)^(1/2)` and compare `isequal`.
"""
struct Monomial{E <: ExactScalar}
    slot::E
    syms::Vector{RawExpression}
    exps::Vector{Rational{Int}}

    function Monomial{E}(
            scalar::Union{ComplexF64, E},
            syms::Vector{RawExpression},
            exps::Vector{Rational{Int}},
        ) where {E <: ExactScalar}
        needs_radical_fold(syms, exps) || return new{E}(pack_scalar(E, scalar), syms, exps)
        scalar, syms, exps = fold_radicals(E, scalar, syms, exps)
        return new{E}(pack_scalar(E, scalar), syms, exps)
    end

    Monomial{E}(::Unchecked, scalar::Union{ComplexF64, E}, m::Monomial) where {E <: ExactScalar} =
        new{E}(pack_scalar(E, scalar), m.syms, m.exps)
end

const EMPTY_SYMS = RawExpression[]
const EMPTY_EXPS = Rational{Int}[]

@inline monomial_scalar(m::Monomial) = unpack_scalar(m.slot)
@inline constant_term(::Type{E}, z) where {E <: ExactScalar} =
    Monomial{E}(as_tier(E, z), EMPTY_SYMS, EMPTY_EXPS)
@inline function fits_small(m::Monomial{BigExactComplex})
    scalar = monomial_scalar(m)
    return scalar isa ComplexF64 || fits_small(scalar)
end
@inline with_scalar(m::Monomial{E}, z::Union{ComplexF64, E}) where {E} =
    Monomial{E}(Unchecked(), z, m)

@inline as_tier(::Type{E}, m::Monomial{E}) where {E <: ExactScalar} = m
@inline as_tier(::Type{E}, m::Monomial) where {E <: ExactScalar} =
    Monomial{E}(Unchecked(), as_tier(E, monomial_scalar(m)), m)
@inline as_tier(::Type{E}, terms::Vector{Monomial{E}}) where {E <: ExactScalar} = terms
as_tier(::Type{E}, terms::Vector{<:Monomial}) where {E <: ExactScalar} =
    Monomial{E}[as_tier(E, m) for m in terms]

@inline function scalar_iszero(m::Monomial)
    s = m.slot
    iszero(s.re) && iszero(s.im) && return true
    return holds_float(s) && iszero(monomial_scalar(m))
end
@inline scalar_is_real(m::Monomial) = iszero(imag(monomial_scalar(m)))
@inline scalar_isequal(a::Monomial, b::Monomial) = isequal(monomial_scalar(a), monomial_scalar(b))
@inline hash_scalar(m::Monomial, h::UInt) = hash(monomial_scalar(m), h)
@inline function normalize_monomial(m::Monomial{E}) where {E}
    scalar = monomial_scalar(m)
    scalar isa ComplexF64 || return m
    return with_scalar(m, normalize_scalar(scalar))
end

@inline function add_scalars(a::Monomial{E}, b::Monomial{E}) where {E}
    sa, sb = a.slot, b.slot
    holds_float(sa) || holds_float(sb) || return with_scalar(a, sa + sb)
    return with_scalar(a, scalar_add(monomial_scalar(a), monomial_scalar(b)))
end
@inline scale_monomial(t::Monomial{E}, z::Union{ComplexF64, E}) where {E} =
    with_scalar(t, scalar_mul(monomial_scalar(t), z))

@inline function monomial_terms(
        ::Type{E}, scalar::CoeffScalar, syms::Vector{RawExpression},
        exps::Vector{Rational{Int}},
    ) where {E <: ExactScalar}
    return Monomial{E}[monomial(E, as_tier(E, scalar), syms, exps)]
end

function sort_factors!(syms::Vector{RawExpression}, exps::Vector{Rational{Int}})
    @inbounds for i in 2:length(syms)
        s = syms[i]; e = exps[i]; k = fkey(s)
        j = i - 1
        while j >= 1 && fkey(syms[j]) > k
            syms[j + 1] = syms[j]; exps[j + 1] = exps[j]
            j -= 1
        end
        j == i - 1 && continue
        syms[j + 1] = s; exps[j + 1] = e
    end
    return syms
end

function merge_factor_list(syms::Vector{RawExpression}, exps::Vector{Rational{Int}})
    n = length(syms)
    n <= 1 && return (syms, exps)
    sort_factors!(syms, exps)
    is_merged(syms, exps) && return (syms, exps)
    osyms = RawExpression[]
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

@inline function is_merged(syms::Vector{RawExpression}, exps::Vector{Rational{Int}})
    @inbounds for i in eachindex(syms)
        iszero(exps[i]) && return false
        i > 1 && syms[i] === syms[i - 1] && return false
    end
    return true
end

"""
    monomial(E, scalar, syms, exps)

The canonical monomial `scalar * ∏ symsᵢ^expsᵢ` from factors in any order. Equal factors
merge, zero exponents drop, several unit phases fold into one and radical exponents fold into
`(0, 1)`. It takes ownership of `syms` and `exps`; `Monomial{E}` itself expects factors that
are already sorted and merged.
"""
function monomial(
        ::Type{E}, scalar::Union{ComplexF64, E}, syms::Vector{RawExpression},
        exps::Vector{Rational{Int}},
    ) where {E <: ExactScalar}
    several_phases = count(is_phase, syms) > 1
    syms, exps = merge_factor_list(syms, exps)
    several_phases && return canonical_phase_monomial!(E, scalar, syms, exps)
    return Monomial{E}(normalize_scalar(scalar), syms, exps)
end

"""
    Poly

A sparse multivariate polynomial over named parameters (a sum of distinct
[`Monomial`](@ref) terms in canonical order), kept off SymbolicUtils hashconsing
and lowered to `Complex{Num}` only at the symbolic boundaries (see `poly_to_num`).
"""
struct Poly
    terms::Union{Vector{Monomial{ExactComplex}}, Vector{Monomial{BigExactComplex}}}
end

@inline fkey(s::SymbolicUtils.BasicSymbolic) = objectid(s)

@inline function same_factors(a::Monomial, b::Monomial)
    length(a.syms) == length(b.syms) || return false
    a.exps == b.exps || return false
    @inbounds for i in eachindex(a.syms)
        a.syms[i] === b.syms[i] || return false
    end
    return true
end

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

function merge_factors(a::Monomial, b::Monomial)
    ia, ib = 1, 1
    na, nb = length(a.syms), length(b.syms)
    syms = RawExpression[]
    exps = Rational{Int}[]
    sizehint!(syms, na + nb); sizehint!(exps, na + nb)
    @inbounds while ia <= na || ib <= nb
        if ib > nb || (ia <= na && fkey(a.syms[ia]) < fkey(b.syms[ib]))
            push!(syms, a.syms[ia]); push!(exps, a.exps[ia]); ia += 1
        elseif ia > na || fkey(a.syms[ia]) > fkey(b.syms[ib])
            push!(syms, b.syms[ib]); push!(exps, b.exps[ib]); ib += 1
        else
            e = a.exps[ia] + b.exps[ib]
            e != 0 && (push!(syms, a.syms[ia]); push!(exps, e))
            ia += 1; ib += 1
        end
    end
    return (syms, exps)
end

function term_mul(a::Monomial{E}, b::Monomial{E}) where {E}
    scalar = scalar_mul(monomial_scalar(a), monomial_scalar(b))
    isempty(a.syms) && return with_scalar(b, scalar)
    isempty(b.syms) && return with_scalar(a, scalar)
    syms, exps = merge_factors(a, b)
    (any(is_phase, a.syms) && any(is_phase, b.syms) && any(is_phase, syms)) ||
        return Monomial{E}(scalar, syms, exps)
    return canonical_phase_monomial!(E, scalar, syms, exps)
end

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
        else
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

function poly_scale(p::Vector{Monomial{E}}, z::Union{ComplexF64, E}) where {E}
    iszero(z) && return Monomial{E}[]
    out = Monomial{E}[]
    sizehint!(out, length(p))
    for t in p
        s = scale_monomial(t, z)
        scalar_iszero(s) || push!(out, s)
    end
    return out
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
