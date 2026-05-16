struct OrderedTerm
    prefactor::CNum
    ops::Vector{QSym}
end

# A Heisenberg pair at one site, with `m` left-operators and `n` right-operators
# satisfying `[right, left] = K` (c-number). The conversion formula reads
#
#   A^m B^n = sum_{k=0}^{min(m,n)} C(m,k) C(n,k) k! (K α)^k Sym(A^{m-k} B^{n-k})
#
# with `α = -1/2` for normal → symmetric and `α = +1/2` for symmetric → normal.
const _WEYL_K_FOCK = _CNUM_ONE       # [Destroy, Create] = 1
const _WEYL_K_PHASE = _CNUM_NEG_IM   # [Momentum, Position] = -i

# Each entry describes one (m, n) group on a single site, ready to feed into
# the combinatorial expansion. Fock pairs distinguish multiple ladder modes on
# the same Hilbert subspace by `op.name`; phase pairs collapse by site, since
# every PhaseSpace subsystem carries one canonical (x, p) pair.
const _FockSite = Tuple{Int, Symbol, Symbol, Int, Int}
const _PhaseSite = Tuple{Int, Symbol, Int, Int}

"""
    normal_to_symmetric(expr) -> QAdd

Convert `expr` from normal to symmetric (Weyl) ordering. Exact inverse of
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
    return _convert_ordering(s, -1 // 2)
end
normal_to_symmetric(op::QSym) = normal_to_symmetric(_single_qadd(_CNUM_ONE, QSym[op]))

"""
    symmetric_to_normal(expr) -> QAdd

Convert `expr` from symmetric (Weyl) to normal ordering. Exact inverse of
[`normal_to_symmetric`](@ref).

# Examples

```jldoctest
julia> hf = FockSpace(:f);

julia> @qnumbers a::Destroy(hf);

julia> symmetric_to_normal(normal_to_symmetric(a' * a))
a' * a

julia> hp = PhaseSpace(:osc);

julia> @qnumbers x::Position(hp) p::Momentum(hp);

julia> symmetric_to_normal(normal_to_symmetric(x * p))
x * p
```

See also [`normal_to_symmetric`](@ref), [`normal_order`](@ref).
"""
function symmetric_to_normal(s::QAdd)
    return _convert_ordering(s, 1 // 2)
end
symmetric_to_normal(op::QSym) = symmetric_to_normal(_single_qadd(_CNUM_ONE, QSym[op]))

function _convert_ordering(s::QAdd, α::Rational)
    d = QTermDict()
    for (term, c) in s.arguments
        _iszero_cnum(c) && continue
        _convert_term!(d, c, term.ops, term.ne, α)
    end
    return QAdd(d, copy(s.indices))
end

function _convert_term!(
        d::QTermDict, c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair}, α::Rational,
    )
    fock_sites = _count_fock_sites(ops)
    phase_sites = _count_phase_sites(ops)
    if isempty(fock_sites) && isempty(phase_sites)
        _addto!(d, ops, c, ne)
        return
    end
    current = OrderedTerm[OrderedTerm(c, ops)]
    for site in fock_sites
        current = _expand_fock_site!(current, site, α)
    end
    for site in phase_sites
        current = _expand_phase_site!(current, site, α)
    end
    for t in current
        _addto!(d, t.ops, t.prefactor, ne)
    end
    return
end

function _count_fock_sites(ops::Vector{QSym})
    counts = Dict{Tuple{Int, Symbol, Symbol}, Tuple{Int, Int}}()
    for op in ops
        if op isa Create
            k = (op.space_index, op.index.name, op.name)
            m, n = get(counts, k, (0, 0))
            counts[k] = (m + 1, n)
        elseif op isa Destroy
            k = (op.space_index, op.index.name, op.name)
            m, n = get(counts, k, (0, 0))
            counts[k] = (m, n + 1)
        end
    end
    sites = _FockSite[]
    for ((si, idxname, name), (m, n)) in counts
        (m > 0 && n > 0) && push!(sites, (si, idxname, name, m, n))
    end
    return sites
end

function _count_phase_sites(ops::Vector{QSym})
    counts = Dict{Tuple{Int, Symbol}, Tuple{Int, Int}}()
    for op in ops
        if op isa Position
            k = (op.space_index, op.index.name)
            m, n = get(counts, k, (0, 0))
            counts[k] = (m + 1, n)
        elseif op isa Momentum
            k = (op.space_index, op.index.name)
            m, n = get(counts, k, (0, 0))
            counts[k] = (m, n + 1)
        end
    end
    sites = _PhaseSite[]
    for ((si, idxname), (m, n)) in counts
        (m > 0 && n > 0) && push!(sites, (si, idxname, m, n))
    end
    return sites
end

@inline function _expand_site!(
        current::Vector{OrderedTerm},
        m::Int, n::Int, K::CNum, α::Rational, remove::F,
    ) where {F}
    Kα = _mul_cnum(K, _to_cnum(α))
    kmax = min(m, n)
    next = OrderedTerm[]
    sizehint!(next, length(current) * (kmax + 1))
    for t in current
        push!(next, OrderedTerm(t.prefactor, copy(t.ops)))
        Kα_pow = _CNUM_ONE
        for k in 1:kmax
            Kα_pow = _mul_cnum(Kα_pow, Kα)
            r = binomial(m, k) * binomial(n, k) * factorial(k)
            coeff = _mul_cnum(_to_cnum(r), Kα_pow)
            new_ops = remove(copy(t.ops), k)
            push!(next, OrderedTerm(_mul_cnum(t.prefactor, coeff), new_ops))
        end
    end
    return next
end

function _expand_fock_site!(current::Vector{OrderedTerm}, site::_FockSite, α::Rational)
    si, idxname, name, m, n = site
    return _expand_site!(
        current, m, n, _WEYL_K_FOCK, α,
        (ops, k) -> _remove_fock_at!(ops, si, idxname, name, k),
    )
end

function _expand_phase_site!(current::Vector{OrderedTerm}, site::_PhaseSite, α::Rational)
    si, idxname, m, n = site
    return _expand_site!(
        current, m, n, _WEYL_K_PHASE, α,
        (ops, k) -> _remove_phase_at!(ops, si, idxname, k),
    )
end

function _remove_fock_at!(
        ops::Vector{QSym}, si::Int, idxname::Symbol, name::Symbol, k::Int,
    )
    creates = 0
    i = 1
    while i <= length(ops) && creates < k
        op = ops[i]
        if op isa Create && op.space_index == si &&
                op.index.name == idxname && op.name == name
            deleteat!(ops, i)
            creates += 1
        else
            i += 1
        end
    end
    destroys = 0
    i = 1
    while i <= length(ops) && destroys < k
        op = ops[i]
        if op isa Destroy && op.space_index == si &&
                op.index.name == idxname && op.name == name
            deleteat!(ops, i)
            destroys += 1
        else
            i += 1
        end
    end
    return ops
end

function _remove_phase_at!(
        ops::Vector{QSym}, si::Int, idxname::Symbol, k::Int,
    )
    positions = 0
    i = 1
    while i <= length(ops) && positions < k
        op = ops[i]
        if op isa Position && op.space_index == si && op.index.name == idxname
            deleteat!(ops, i)
            positions += 1
        else
            i += 1
        end
    end
    momenta = 0
    i = 1
    while i <= length(ops) && momenta < k
        op = ops[i]
        if op isa Momentum && op.space_index == si && op.index.name == idxname
            deleteat!(ops, i)
            momenta += 1
        else
            i += 1
        end
    end
    return ops
end
