"""
    Monomial

One term of a parameter polynomial: `scalar * ∏ symᵢ^expᵢ`. Factors are sorted by
`objectid` and deduplicated; `Rational{Int}` exponents let radicals of a single
atom merge (`sqrt(p)*sqrt(p) = p`).
"""
# Keep exact Gaussian rationals alongside the floating-point fast path.  A rational
# scalar is deliberately narrower than an arbitrary `Complex{Rational}`: this keeps
# the polynomial representation concrete while covering the exact literals produced
# by Julia's `//` operator.
const ExactComplex = Complex{Rational{Int}}
const CoeffScalar = Union{ComplexF64, ExactComplex}

struct Monomial
    scalar::CoeffScalar
    syms::Vector{SymbolicUtils.BasicSymbolic}   # sorted by objectid, distinct
    exps::Vector{Rational{Int}}                 # matching nonzero exponents
end

@inline normalize_scalar(z::ComplexF64) = z + complex(0.0, 0.0)
@inline normalize_scalar(z::ExactComplex) = z

# Native integer/Gaussian-integer factors (1, -1, and `im`) do not introduce
# inexactness when multiplied by an exact scalar.  Genuine non-integral floats do.
@inline function integer_scalar(z::ComplexF64)
    re, im = real(z), imag(z)
    isfinite(re) && isfinite(im) && isinteger(re) && isinteger(im) || return nothing
    abs(re) <= typemax(Int) && abs(im) <= typemax(Int) || return nothing
    return ExactComplex(Int(re) // 1, Int(im) // 1)
end

@inline function scalar_inv(z::ComplexF64)
    exact = integer_scalar(z)
    return exact === nothing ? inv(z) : inv(exact)
end
@inline scalar_inv(z::ExactComplex) = inv(z)

@inline function scalar_mul(a::ExactComplex, b::ComplexF64)
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(a * b) : a * ib
end
@inline scalar_mul(a::ComplexF64, b::ExactComplex) = scalar_mul(b, a)
@inline scalar_mul(a::ExactComplex, b::ExactComplex) = a * b
@inline scalar_mul(a::ComplexF64, b::ComplexF64) = normalize_scalar(a * b)

@inline function scalar_add(a::ExactComplex, b::ComplexF64)
    ib = integer_scalar(b)
    return ib === nothing ? normalize_scalar(a + b) : a + ib
end
@inline scalar_add(a::ComplexF64, b::ExactComplex) = scalar_add(b, a)
@inline scalar_add(a::ExactComplex, b::ExactComplex) = a + b
@inline scalar_add(a::ComplexF64, b::ComplexF64) = normalize_scalar(a + b)

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

function term_mul(a::Monomial, b::Monomial)
    scalar = scalar_mul(a.scalar, b.scalar)
    isempty(a.syms) && return Monomial(scalar, b.syms, b.exps)
    isempty(b.syms) && return Monomial(scalar, a.syms, a.exps)
    phase_a = phase_factor_index(a.syms)
    phase_b = phase_factor_index(b.syms)
    if phase_a != 0 && phase_b != 0
        # The common inverse pair stays on the ordinary identity merge: no symbolic
        # argument arithmetic and no additional allocation.
        if a.syms[phase_a] === b.syms[phase_b] &&
                a.exps[phase_a] == -b.exps[phase_b]
            se = merge_factors(a.syms, a.exps, b.syms, b.exps)
            return Monomial(scalar, se[1], se[2])
        end
        if length(a.syms) == 1 && length(b.syms) == 1 &&
                a.syms[phase_a] === b.syms[phase_b]
            exponent = a.exps[phase_a] + b.exps[phase_b]
            return scaled_phase_monomial(
                scalar, Num(a.syms[phase_a]), exponent,
            )
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
    se = merge_factors(a.syms, a.exps, b.syms, b.exps)
    return Monomial(scalar, se[1], se[2])
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
            s = scalar_add(terms[w].scalar, t.scalar)
            terms[w] = Monomial(s, terms[w].syms, terms[w].exps)
        else
            w += 1
            terms[w] = Monomial(normalize_scalar(t.scalar), t.syms, t.exps)
        end
    end
    resize!(terms, w)
    filter!(t -> t.scalar != 0, terms)
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
            t.scalar != 0 && push!(out, t)
        elseif ip > np || term_less(q[iq], p[ip])
            t = q[iq]; iq += 1
            t.scalar != 0 && push!(out, t)
        else   # same factor set: sum scalars, drop exact cancellations
            s = scalar_add(p[ip].scalar, q[iq].scalar)
            s != 0 && push!(out, Monomial(s, p[ip].syms, p[ip].exps))
            ip += 1; iq += 1
        end
    end
    return out
end

function poly_mul(p::Vector{Monomial}, q::Vector{Monomial})
    if length(p) == 1 && length(q) == 1
        t = term_mul(p[1], q[1])
        return Monomial[Monomial(normalize_scalar(t.scalar), t.syms, t.exps)]
    end
    out = Monomial[]
    sizehint!(out, length(p) * length(q))
    for a in p, b in q
        push!(out, term_mul(a, b))
    end
    return canonical_terms!(out)
end

# Scale every term; preserves canonical order (factors unchanged).
function poly_scale(p::Vector{Monomial}, z::CoeffScalar)
    iszero(z) && return Monomial[]
    return Monomial[Monomial(scalar_mul(t.scalar, z), t.syms, t.exps) for t in p]
end

function Base.isequal(a::Poly, b::Poly)
    length(a.terms) == length(b.terms) || return false
    @inbounds for i in eachindex(a.terms)
        ta, tb = a.terms[i], b.terms[i]
        (isequal(ta.scalar, tb.scalar) && same_factors(ta, tb)) || return false
    end
    return true
end
Base.:(==)(a::Poly, b::Poly) = isequal(a, b)
function Base.hash(p::Poly, h::UInt)
    @inbounds for t in p.terms
        h = hash(t.scalar, h)
        for i in eachindex(t.syms)
            h = hash(t.exps[i], hash(fkey(t.syms[i]), h))
        end
    end
    return hash(:Poly, h)
end
