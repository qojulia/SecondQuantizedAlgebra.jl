struct OrderedTerm
    prefactor::CNum
    ops::Vector{Op}
end

struct HeisenbergPair
    raising::OpKind
    lowering::OpKind
    K::CNum
end

const HEISENBERG_PAIRS = (
    HeisenbergPair(OP_CREATE, OP_DESTROY, CNUM_ONE),
    HeisenbergPair(OP_POSITION, OP_MOMENTUM, CNUM_NEG_IM),
)

"""
    normal_to_symmetric(expr)

Rewrite a normal-ordered expression in symmetric (Weyl) order.

The result is written with the same monomials, but each now stands for its fully
symmetrized product: `a' * a` reads as ``(a^\\dagger a + a a^\\dagger)/2`` and
`x * p` as ``(xp + px)/2``. This is the form phase-space methods such as the
truncated Wigner approximation work with. Bosonic modes and quadratures are
converted; other operators pass through unchanged. The inverse is
[`symmetric_to_normal`](@ref).

# Examples

```jldoctest
julia> hf = FockSpace(:f);

julia> @qnumbers a::Destroy(hf);

julia> normal_to_symmetric(a' * a)
-1//2 + a' * a

julia> hp = PhaseSpace(:osc);

julia> @qnumbers x::Position(hp) p::Momentum(hp);

julia> normal_to_symmetric(x * p)
1//2im + x * p
```

See also [`symmetric_to_normal`](@ref), [`normal_order`](@ref).
"""
function normal_to_symmetric(s::QAdd)
    return convert_ordering(s, -1 // 2)
end
normal_to_symmetric(op::QSym) = normal_to_symmetric(to_qadd(op))

"""
    symmetric_to_normal(expr)

Rewrite a symmetric (Weyl) ordered expression in normal order.

Each monomial of `expr` is read as its symmetrized product, as returned by
[`normal_to_symmetric`](@ref), and expanded back into normal order. The two
functions are exact inverses.

# Examples

```jldoctest
julia> hf = FockSpace(:f);

julia> @qnumbers a::Destroy(hf);

julia> symmetric_to_normal(a' * a)
1//2 + a' * a

julia> normal_to_symmetric(symmetric_to_normal(a' * a))
a' * a
```

See also [`normal_to_symmetric`](@ref), [`normal_order`](@ref).
"""
function symmetric_to_normal(s::QAdd)
    return convert_ordering(s, 1 // 2)
end
symmetric_to_normal(op::QSym) = symmetric_to_normal(to_qadd(op))

function convert_ordering(s::QAdd, α::Rational)
    d = QTermDict()
    for (term, c) in s.arguments
        convert_term!(d, c, term.ops, term.ne, α)
    end
    return QAdd(d, copy(s.indices))
end

"""
    convert_term!(d, c, ops, ne, α)

Expand one term under Weyl/symmetric ordering and accumulate into `d`.

Every site holding both operators of a pair in `HEISENBERG_PAIRS` is expanded by
`expand_site!`; terms without such a site pass through unchanged. A new
Heisenberg pair kind is one more entry in `HEISENBERG_PAIRS`.
"""
function convert_term!(
        d::QTermDict, c::CNum, ops::Vector{Op},
        ne::Vector{NonEqualPair}, α::Rational,
    )
    current = OrderedTerm[OrderedTerm(c, ops)]
    for pair in HEISENBERG_PAIRS
        for (key, m, n) in conjugate_sites(ops, pair)
            current = expand_site!(current, pair, key, m, n, α)
        end
    end
    for t in current
        addto!(d, t.ops, t.prefactor, ne)
    end
    return
end

function conjugate_sites(ops::Vector{Op}, pair::HeisenbergPair)
    counts = Dict{SiteKey, Tuple{Int, Int}}()
    for op in ops
        raising = op.kind === pair.raising
        (raising || op.kind === pair.lowering) || continue
        key = site_key(op)
        m, n = get(counts, key, (0, 0))
        counts[key] = raising ? (m + 1, n) : (m, n + 1)
    end
    sites = Tuple{SiteKey, Int, Int}[]
    for (key, (m, n)) in counts
        (m > 0 && n > 0) && push!(sites, (key, m, n))
    end
    return sites
end

function expand_site!(
        current::Vector{OrderedTerm}, pair::HeisenbergPair,
        key::SiteKey, m::Int, n::Int, α::Rational,
    )
    Kα = pair.K * to_cnum(α)
    kmax = min(m, n)
    next = OrderedTerm[]
    sizehint!(next, length(current) * (kmax + 1))
    for t in current
        push!(next, t)
        Kα_pow = CNUM_ONE
        for k in 1:kmax
            Kα_pow = Kα_pow * Kα
            coeff = to_cnum(binomial(m, k) * binomial(n, k) * factorial(k)) * Kα_pow
            new_ops = remove_at!(copy(t.ops), pair.raising, key, k)
            remove_at!(new_ops, pair.lowering, key, k)
            push!(next, OrderedTerm(t.prefactor * coeff, new_ops))
        end
    end
    return next
end

function remove_at!(ops::Vector{Op}, kind::OpKind, key::SiteKey, k::Int)
    removed = 0
    i = 1
    while i <= length(ops) && removed < k
        op = ops[i]
        if op.kind === kind && site_key(op) == key
            deleteat!(ops, i)
            removed += 1
        else
            i += 1
        end
    end
    return ops
end
