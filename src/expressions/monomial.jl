"""
    Monomial

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`).
"""
struct Monomial
    scalar::ComplexF64
    syms::Vector{SymbolicUtils.BasicSymbolic}   # sorted by objectid, distinct
    exps::Vector{Rational{Int}}                 # matching nonzero exponents
end

"""
    Poly

A native sparse multivariate polynomial over named parameters (a sum of distinct
[`Monomial`](@ref) terms in canonical order), kept off SymbolicUtils hashconsing
and lowered to `Complex{Num}` only at the symbolic boundaries (see `_poly_to_num`).
"""
struct Poly
    terms::Vector{Monomial}
end

# Factor identity key: SymbolicUtils hashconses, so `objectid`/`===` are exact and
# type-stable factor identity (unlike `hash`/`isequal` on abstract `BasicSymbolic`).
@inline _fkey(s::SymbolicUtils.BasicSymbolic) = objectid(s)

@inline function _same_factors(a::Monomial, b::Monomial)
    length(a.syms) == length(b.syms) || return false
    a.exps == b.exps || return false
    @inbounds for i in eachindex(a.syms)
        a.syms[i] === b.syms[i] || return false
    end
    return true
end

# Total order on monomial factor sets (factors are pre-sorted by objectid within
# each monomial), giving a canonical term order for `Poly` equality / hashing.
function _term_less(a::Monomial, b::Monomial)
    la, lb = length(a.syms), length(b.syms)
    la != lb && return la < lb
    @inbounds for i in 1:la
        ha, hb = _fkey(a.syms[i]), _fkey(b.syms[i])
        ha != hb && return ha < hb
        a.exps[i] != b.exps[i] && return a.exps[i] < b.exps[i]
    end
    return false
end

# Merge two sorted factor lists, summing exponents and dropping cancellations.
function _merge_factors(syma, expa, symb, expb)
    ia, ib = 1, 1
    na, nb = length(syma), length(symb)
    syms = SymbolicUtils.BasicSymbolic[]
    exps = Rational{Int}[]
    sizehint!(syms, na + nb); sizehint!(exps, na + nb)
    @inbounds while ia <= na || ib <= nb
        if ib > nb || (ia <= na && _fkey(syma[ia]) < _fkey(symb[ib]))
            push!(syms, syma[ia]); push!(exps, expa[ia]); ia += 1
        elseif ia > na || _fkey(syma[ia]) > _fkey(symb[ib])
            push!(syms, symb[ib]); push!(exps, expb[ib]); ib += 1
        else
            e = expa[ia] + expb[ib]
            e != 0 && (push!(syms, syma[ia]); push!(exps, e))
            ia += 1; ib += 1
        end
    end
    return (syms, exps)
end

function _term_mul(a::Monomial, b::Monomial)
    isempty(a.syms) && return Monomial(a.scalar * b.scalar, b.syms, b.exps)
    isempty(b.syms) && return Monomial(a.scalar * b.scalar, a.syms, a.exps)
    se = _merge_factors(a.syms, a.exps, b.syms, b.exps)
    return Monomial(a.scalar * b.scalar, se[1], se[2])
end

# Sort terms into canonical order, merge like-factor terms, drop zero scalars.
function _canonical_terms!(terms::Vector{Monomial})
    isempty(terms) && return terms
    sort!(terms; lt = _term_less)
    w = 0
    @inbounds for r in eachindex(terms)
        t = terms[r]
        if w > 0 && _same_factors(terms[w], t)
            s = terms[w].scalar + t.scalar + complex(0.0, 0.0)
            terms[w] = Monomial(s, terms[w].syms, terms[w].exps)
        else
            w += 1
            terms[w] = Monomial(t.scalar + complex(0.0, 0.0), t.syms, t.exps)
        end
    end
    resize!(terms, w)
    filter!(t -> t.scalar != 0, terms)
    return terms
end

# Sorted merge of two canonical term lists, dropping zero-scalar terms so the result
# stays canonical and zero-free (a stray zero would break Poly equality/hashing).
function _poly_add(p::Vector{Monomial}, q::Vector{Monomial})
    out = Monomial[]
    sizehint!(out, length(p) + length(q))
    ip, iq = 1, 1
    np, nq = length(p), length(q)
    @inbounds while ip <= np || iq <= nq
        if iq > nq || (ip <= np && _term_less(p[ip], q[iq]))
            t = p[ip]; ip += 1
            t.scalar != 0 && push!(out, t)
        elseif ip > np || _term_less(q[iq], p[ip])
            t = q[iq]; iq += 1
            t.scalar != 0 && push!(out, t)
        else   # same factor set: sum scalars, drop exact cancellations
            s = p[ip].scalar + q[iq].scalar + complex(0.0, 0.0)
            s != 0 && push!(out, Monomial(s, p[ip].syms, p[ip].exps))
            ip += 1; iq += 1
        end
    end
    return out
end

function _poly_mul(p::Vector{Monomial}, q::Vector{Monomial})
    if length(p) == 1 && length(q) == 1
        t = _term_mul(p[1], q[1])
        return Monomial[Monomial(t.scalar + complex(0.0, 0.0), t.syms, t.exps)]
    end
    out = Monomial[]
    sizehint!(out, length(p) * length(q))
    for a in p, b in q
        push!(out, _term_mul(a, b))
    end
    return _canonical_terms!(out)
end

# Scale every term; preserves canonical order (factors unchanged).
function _poly_scale(p::Vector{Monomial}, z::ComplexF64)
    iszero(z) && return Monomial[]
    return Monomial[Monomial(t.scalar * z + complex(0.0, 0.0), t.syms, t.exps) for t in p]
end

function Base.isequal(a::Poly, b::Poly)
    length(a.terms) == length(b.terms) || return false
    @inbounds for i in eachindex(a.terms)
        ta, tb = a.terms[i], b.terms[i]
        (isequal(ta.scalar, tb.scalar) && _same_factors(ta, tb)) || return false
    end
    return true
end
Base.:(==)(a::Poly, b::Poly) = isequal(a, b)
function Base.hash(p::Poly, h::UInt)
    @inbounds for t in p.terms
        h = hash(t.scalar, h)
        for i in eachindex(t.syms)
            h = hash(t.exps[i], hash(_fkey(t.syms[i]), h))
        end
    end
    return hash(:Poly, h)
end
