"""
    Monomial

One term of a parameter polynomial: a native `ComplexF64` scalar times a product
of named symbolic parameters with integer exponents, `scalar * ∏ symᵢ^expᵢ`.
Factors are stored sorted by their `objectid` and deduplicated.
"""
struct Monomial
    scalar::ComplexF64
    syms::Vector{SymbolicUtils.BasicSymbolic}   # sorted by objectid, distinct
    exps::Vector{Int}                           # matching nonzero exponents
end

"""
    Poly

A native sparse multivariate polynomial over named parameters: a sum of
[`Monomial`](@ref) terms with distinct factor sets, sorted into a canonical
order. Tier-2 coefficient representation: products and sums of parameter
polynomials are native (scalar arithmetic plus integer-exponent merges), so the
common symbolic coefficient never routes through SymbolicUtils hashconsing. It
lowers to `Complex{Num}` only at the symbolic boundaries (see `_poly_to_num`).
"""
struct Poly
    terms::Vector{Monomial}
end

# Factor identity key. SymbolicUtils hashconses (the SU4 default), so two
# structurally-equal symbols are the *same object*; `objectid` is therefore a
# stable identity key and `===` is exact factor equality. Both are builtins that
# take `Any` without dispatching, unlike `hash(::BasicSymbolic)` /
# `isequal(::BasicSymbolic, ::BasicSymbolic)` whose abstract element type forces
# runtime dispatch. Using them keeps every polynomial operation fully type-stable.
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
    exps = Int[]
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

_term_mul(a::Monomial, b::Monomial) =
    (se = _merge_factors(a.syms, a.exps, b.syms, b.exps); Monomial(a.scalar * b.scalar, se[1], se[2]))

# Sort terms into canonical order, merge like-factor terms, drop zero scalars.
function _canonical_terms(terms::Vector{Monomial})
    isempty(terms) && return terms
    sort!(terms; lt = _term_less)
    out = Monomial[]
    @inbounds for t in terms
        if !isempty(out) && _same_factors(out[end], t)
            s = out[end].scalar + t.scalar
            out[end] = Monomial(s + complex(0.0, 0.0), t.syms, t.exps)
        else
            push!(out, Monomial(t.scalar + complex(0.0, 0.0), t.syms, t.exps))
        end
    end
    filter!(t -> t.scalar != 0, out)
    return out
end

_poly_add(p::Vector{Monomial}, q::Vector{Monomial}) = _canonical_terms(vcat(p, q))

function _poly_mul(p::Vector{Monomial}, q::Vector{Monomial})
    out = Monomial[]
    sizehint!(out, length(p) * length(q))
    for a in p, b in q
        push!(out, _term_mul(a, b))
    end
    return _canonical_terms(out)
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
