# Vacuum expectation values. The reference state is the product over subspaces of
# each subspace's own lowest state, which the algebra already knows: Fock |0⟩,
# NLevelSpace |g⟩, PhaseSpace the oscillator ground state fixed by x = (a+a†)/√2,
# Pauli/Spin the top weight |m = S⟩ (σz = +1), the first basis state of the
# `SpinBasis` that `to_numeric` builds, and CollectiveNLevelSpace all N particles
# in level 1.

const _VacSize = Union{Nothing, Real, AbstractDict}
const _VacAmp = Complex{Rational{Int}}

# Each free-index pair the algebra cannot resolve doubles the branch count, so cap
# the recursion rather than letting a pathological term run away.
const _VAC_MAX_SPLITS = 8

"""
    Vacuum()

The reference state ``|0\\rangle`` of a symbolic Hilbert space, for use as the
state argument of [`expect`](@ref expect(::QAdd, ::Vacuum)).

It is the product over subspaces of each subspace's own lowest state, so it needs
no dimensions, no cutoff, and no numeric backend:

| Subspace | State |
|:--|:--|
| [`FockSpace`](@ref) | the vacuum ``\\|0\\rangle`` |
| [`NLevelSpace`](@ref) | the declared `ground_state` ``\\|g\\rangle`` |
| [`PhaseSpace`](@ref) | the oscillator ground state implied by ``x = (a+a^\\dagger)/\\sqrt2`` |
| [`PauliSpace`](@ref) | the ``\\sigma_z = +1`` eigenstate |
| [`SpinSpace`](@ref) | the top weight ``\\|m = S\\rangle`` |
| [`CollectiveNLevelSpace`](@ref) | all ``N`` particles in level 1 |

Apart from an [`NLevelSpace`](@ref) whose ground state is not level 1, this is
the state that [`to_numeric`](@ref) represents as basis state 1 of every
subsystem.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> expect(a * a', Vacuum())
1
```

See also [`expect`](@ref expect(::QAdd, ::Vacuum)).
"""
struct Vacuum end

expect(x::_ScalarLike, ::Vacuum; spin::_VacSize = nothing, particles::_VacSize = nothing) =
    to_num(_to_cnum(x))
expect(op::QSym, v::Vacuum; spin::_VacSize = nothing, particles::_VacSize = nothing) =
    expect(_single_qadd(_CNUM_ONE, Op[op]), v; spin = spin, particles = particles)

"""
    expect(op, ::Vacuum; spin = nothing, particles = nothing) -> Complex{Num}

Vacuum expectation value ``\\langle 0 | op | 0 \\rangle`` as an exact symbolic
scalar, with no basis truncation and no numeric backend.

`op` is normal-ordered first, then every term is read against [`Vacuum`](@ref):
each same-site block is contracted on its subspace's reference state and blocks
on distinct sites factorize.

Because every state reachable from the vacuum is ``|\\psi\\rangle = O|0\\rangle``,
the bilinear ``\\langle\\psi|A|\\varphi\\rangle`` is
`expect(Oψ' * A * Oφ, Vacuum())`. Matrix elements, state overlaps, norms, and
zero-point energies are all this one call.

Two subspaces carry a size the algebra does not: pass `spin` for a
[`SpinSpace`](@ref) and `particles` for a [`CollectiveNLevelSpace`](@ref). Each
takes a single value for every such subspace, or an `AbstractDict` mapping
subspace position ([`acts_on`](@ref)) to a value when they differ.

A [`Σ`](@ref) contributes its index range once the operators are contracted. When
the summand does not reduce to a range factor, because the coefficient still
depends on the bound index or the scope carries a `≠` constraint, the result
keeps the sum symbolically as an [`is_indexed_sum`](@ref) node.

Two free indices the algebra cannot separate are resolved by evaluating both the
same-site and the distinct-site reading. When they agree the common value is
returned; when they differ the result carries the exact case split weighted by
[`kronecker_delta`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> expect(a' * a, Vacuum())
0

julia> expect(a^2 * (a')^2, Vacuum())
2

julia> expect(a^2 * (a')^3, Vacuum())
0
```

```jldoctest
julia> h = NLevelSpace(:atom, 2);

julia> σ(i, j) = Transition(h, :σ, i, j);

julia> expect(σ(1, 2) * σ(2, 1), Vacuum())
1

julia> expect(σ(2, 1) * σ(1, 2), Vacuum())
0
```

```jldoctest
julia> h = SpinSpace(:S);

julia> Sx = Spin(h, :S, 1); Sy = Spin(h, :S, 2); Sz = Spin(h, :S, 3);

julia> expect(Sx^2 + Sy^2 + Sz^2, Vacuum(); spin = 3 // 2)
3.75
```

See also [`Vacuum`](@ref), [`normal_order`](@ref), [`numeric_average`](@ref).
"""
function expect(q::QAdd, ::Vacuum; spin::_VacSize = nothing, particles::_VacSize = nothing)
    return to_num(_vac_qadd(normal_order(q), spin, particles, 0))
end

function _vac_qadd(nq::QAdd, spin::_VacSize, particles::_VacSize, depth::Int)
    total = _CNUM_ZERO
    for (term, c) in nq.arguments
        total = _add_cnum(total, _vac_term(term, c, nq.indices, spin, particles, depth))
    end
    return total
end

function _vac_term(
        term::QTerm, c::CNum, scope::Vector{Index},
        spin::_VacSize, particles::_VacSize, depth::Int,
    )
    pair = _vac_undetermined_pair(term.ops, term.ne)
    pair === nothing || return _vac_split(term, c, scope, pair, spin, particles, depth)
    f = _vac_ops(term.ops, term.ne, spin, particles)
    _iszero_cnum(f) && return _CNUM_ZERO
    return _vac_scope(_mul_cnum(c, f), term, c, scope)
end

# The first operator pair whose site relation the algebra could not resolve. Both
# must carry a real index; a bare operator paired against an indexed one on the
# same subspace has no site label to identify, so there is nothing to split on.
function _vac_undetermined_pair(ops::Vector{Op}, ne::Vector{NonEqualPair})
    n = length(ops)
    for i in 1:(n - 1), j in (i + 1):n
        _site_compare(ops[i], ops[j], ne) === Undetermined || continue
        α, β = ops[i].index, ops[j].index
        (has_index(α) && has_index(β)) || throw(
            ArgumentError(
                "vacuum expectation cannot factorize `$(ops[i])` against `$(ops[j])`: one of " *
                    "them carries no site index, so the algebra cannot tell whether they act on " *
                    "the same site. Give both an index, or keep them on separate subspaces."
            )
        )
        return (α, β)
    end
    return nothing
end

# Split on `α == β` versus `α ≠ β`. Both branches re-canonicalize (the diagonal
# one also renames the coefficient), so each resolves the pair and recursion
# terminates. Agreeing branches collapse, which is why an expression like
# ⟨σᵢᵍᵍ σⱼᵍᵍ⟩ comes back delta-free.
function _vac_split(
        term::QTerm, c::CNum, scope::Vector{Index}, pair::NonEqualPair,
        spin::_VacSize, particles::_VacSize, depth::Int,
    )
    depth < _VAC_MAX_SPLITS || throw(
        ArgumentError(
            "vacuum expectation gave up after $(_VAC_MAX_SPLITS) free-index case splits; " *
                "declare the site relations up front with `assume_distinct_index`"
        )
    )
    α, β = pair
    (α in scope || β in scope) && throw(
        ArgumentError(
            "vacuum expectation cannot split the bound index `$(index_name(α in scope ? α : β))` " *
                "against `$(index_name(α in scope ? β : α))`: a summation index must be resolved " *
                "by `Σ`, not by a case split"
        )
    )
    q = QAdd(QTermDict(_copy_key(term) => c), copy(scope))
    diag = _vac_qadd(normal_order(change_index(q, β, α)), spin, particles, depth + 1)
    off = _vac_qadd(normal_order(assume_distinct_index(q, [pair])), spin, particles, depth + 1)
    isequal(diag, off) && return diag
    δ = _to_cnum(kronecker_delta(α, β))
    return _add_cnum(
        _mul_cnum(δ, diag), _mul_cnum(_add_cnum(_CNUM_ONE, _neg_cnum(δ)), off),
    )
end

# Contracting the operators erases a bound index, so one that reached the term
# through an operator now contributes its range. A term that never depended on the
# index absorbed that range when the sum was built, so scaling it again would
# double-count. What is left over stays a symbolic sum.
function _vac_scope(value::CNum, term::QTerm, c::CNum, scope::Vector{Index})
    isempty(scope) && return value
    out = value
    kept = Index[]
    for idx in scope
        dep_c = _depends_on_index_term(c, _EMPTY_OPS, idx)
        dep_ops = any(op -> op.index == idx, term.ops)
        (dep_c || dep_ops) || continue
        if !dep_c && !any(p -> p[1] == idx || p[2] == idx, term.ne)
            out = _mul_cnum(out, _to_cnum(index_range(idx)))
        else
            push!(kept, idx)
        end
    end
    isempty(kept) && return out
    ne = NonEqualPair[p for p in term.ne if any(i -> p[1] == i || p[2] == i, kept)]
    return _vac_sum_node(out, kept, ne)
end

function _vac_sum_node(value::CNum, indices::Vector{Index}, ne::Vector{NonEqualPair})
    z = to_num(value)
    re, im_ = real(z), imag(z)
    wrap(part) = _iszero_num(part) ? _NUM_ZERO : Num(indexed_sum(part, indices; non_equal = ne))
    return _cnum(wrap(re), wrap(im_))
end

# Operators on distinct sites factorize; same-site runs are contiguous in canonical
# form, so one pass over the run boundaries splits the product into site blocks.
function _vac_ops(
        ops::Vector{Op}, ne::Vector{NonEqualPair}, spin::_VacSize, particles::_VacSize,
    )
    n = length(ops)
    n == 0 && return _CNUM_ONE
    acc = _CNUM_ONE
    i = 1
    while i <= n
        j = i
        while j < n && _site_compare(ops[j], ops[j + 1], ne) === Equal
            j += 1
        end
        acc = _mul_cnum(acc, _vac_site(ops, i, j, spin, particles))
        _iszero_cnum(acc) && return _CNUM_ZERO
        i = j + 1
    end
    return acc
end

function _vac_site(
        ops::Vector{Op}, i::Int, j::Int, spin::_VacSize, particles::_VacSize,
    )
    k = ops[i].kind
    if k === OP_DESTROY || k === OP_CREATE
        return _CNUM_ZERO                    # canonical (a†)ᵐaⁿ with m+n ≥ 1 kills |0⟩
    elseif k === OP_TRANSITION
        return _vac_transition(ops, i, j)
    elseif k === OP_PAULI
        return _vac_pauli(ops, i, j)
    elseif k === OP_POSITION || k === OP_MOMENTUM
        return _vac_quadrature(ops, i, j)
    elseif k === OP_SPIN
        return _vac_spin(ops, i, j, _vac_half_integer(ops[i], spin, "spin", "`spin = 1//2`"))
    else
        return _vac_collective(
            ops, i, j, _vac_count(ops[i], particles, "particles", "`particles = 4`"),
        )
    end
end

# ⟨g| σ^{a₁b₁} … σ^{aₙbₙ} |g⟩ = δ_{a₁g} (∏ δ_{bₖaₖ₊₁}) δ_{bₙg}.
function _vac_transition(ops::Vector{Op}, i::Int, j::Int)
    ops[i].l1 == ops[i].g || return _CNUM_ZERO
    for m in i:(j - 1)
        ops[m].l2 == ops[m + 1].l1 || return _CNUM_ZERO
    end
    return ops[j].l2 == ops[j].g ? _CNUM_ONE : _CNUM_ZERO
end

# Fold σⱼσₖ = δⱼₖI + iϵⱼₖₗσₗ over the run; `axis == 0` is the identity.
function _vac_pauli(ops::Vector{Op}, i::Int, j::Int)
    c = _CNUM_ONE
    axis = 0
    for m in i:j
        b = Int(ops[m].l1)
        if axis == 0
            axis = b
        elseif axis == b
            axis = 0
        else
            c = _mul_cnum(c, _to_cnum(im * _levi_civita[axis][b]))
            axis = 6 - axis - b
        end
    end
    return (axis == 0 || axis == 3) ? c : _CNUM_ZERO
end

@inline function _vac_add!(d::Dict{K, _VacAmp}, k::K, v::_VacAmp) where {K}
    iszero(v) && return d
    d[k] = get(d, k, zero(_VacAmp)) + v
    return d
end

# Rescaled Fock basis |n⟩′ = |n⟩/√(n!) makes the ladder action exactly rational:
# ã = a+a† sends |n⟩′ to |n-1⟩′ + (n+1)|n+1⟩′, and p̃ = i(a†-a) to
# i((n+1)|n+1⟩′ - |n-1⟩′). The 1/√2 per quadrature is pulled out as 2^(-N/2),
# rational because an odd-length product has zero amplitude on |0⟩ anyway.
function _vac_quadrature(ops::Vector{Op}, i::Int, j::Int)
    n = j - i + 1
    isodd(n) && return _CNUM_ZERO
    ket = Dict{Int, _VacAmp}(0 => one(_VacAmp))
    for k in j:-1:i
        p = ops[k].kind === OP_MOMENTUM
        nxt = Dict{Int, _VacAmp}()
        for (m, amp) in ket
            up = p ? im * amp : amp
            _vac_add!(nxt, m + 1, up * (m + 1))
            m > 0 && _vac_add!(nxt, m - 1, p ? -up : amp)
        end
        ket = nxt
    end
    amp = get(ket, 0, zero(_VacAmp))
    iszero(amp) && return _CNUM_ZERO
    return _to_cnum(amp * (1 // (1 << (n ÷ 2))))
end

# Rescaled weight basis |m⟩′ with S⁻|m⟩′ = |m-1⟩′ and
# S⁺|m⟩′ = (S(S+1) - m(m+1))|m+1⟩′, so every amplitude stays rational. The
# reference |S,S⟩ is |S⟩′, and S⁺'s coefficient vanishes there on its own; only
# S⁻ at the bottom weight needs the explicit cut.
function _vac_spin(ops::Vector{Op}, i::Int, j::Int, S::Rational{Int})
    casimir = S * (S + 1)
    ket = Dict{Rational{Int}, _VacAmp}(S => one(_VacAmp))
    for k in j:-1:i
        axis = ops[k].l1
        nxt = Dict{Rational{Int}, _VacAmp}()
        for (m, amp) in ket
            if axis == 3
                _vac_add!(nxt, m, amp * m)
                continue
            end
            # Sx = (S⁺+S⁻)/2, Sy = (S⁺-S⁻)/2i
            cp, cm = axis == 1 ? (_VacAmp(1 // 2, 0), _VacAmp(1 // 2, 0)) :
                (_VacAmp(0, -1 // 2), _VacAmp(0, 1 // 2))
            _vac_add!(nxt, m + 1, cp * amp * (casimir - m * (m + 1)))
            m > -S && _vac_add!(nxt, m - 1, cm * amp)
        end
        ket = nxt
    end
    amp = get(ket, S, zero(_VacAmp))
    return iszero(amp) ? _CNUM_ZERO : _to_cnum(amp)
end

# Sⁱʲ = bᵢ†bⱼ on the occupation basis, rescaled by ∏√(nₖ!) so that bₖ|n⟩′ = |n-1ₖ⟩′
# and bₖ†|n⟩′ = (nₖ+1)|n+1ₖ⟩′. The reference |N,0,…⟩ picks up the same rescaling on
# the bra and the ket, so it cancels and the amplitude is read off directly.
function _vac_collective(ops::Vector{Op}, i::Int, j::Int, N::Int)
    nlev = 1
    for k in i:j
        nlev = max(nlev, Int(ops[k].l1), Int(ops[k].l2))
    end
    start = zeros(Int, nlev)
    start[1] = N
    ket = Dict{Vector{Int}, _VacAmp}(start => one(_VacAmp))
    for k in j:-1:i
        a, b = Int(ops[k].l1), Int(ops[k].l2)
        nxt = Dict{Vector{Int}, _VacAmp}()
        for (occ, amp) in ket
            occ[b] > 0 || continue
            new_occ = copy(occ)
            new_occ[b] -= 1
            coeff = new_occ[a] + 1
            new_occ[a] += 1
            _vac_add!(nxt, new_occ, amp * coeff)
        end
        ket = nxt
    end
    amp = get(ket, start, zero(_VacAmp))
    return iszero(amp) ? _CNUM_ZERO : _to_cnum(amp)
end

# `spin` / `particles` accept one value for every matching subspace, or a Dict
# keyed by subspace position when they differ.
function _vac_size(op::Op, spec::_VacSize, name::String, hint::String)
    spec === nothing && throw(
        ArgumentError(
            "a vacuum expectation of this operator needs `$name`, which its Hilbert space does " *
                "not carry: pass $hint, or a `Dict` from subspace position to value when they differ"
        )
    )
    spec isa AbstractDict || return spec
    si = Int(op.space_index)
    haskey(spec, si) || throw(ArgumentError("`$name` has no entry for subspace $si"))
    return spec[si]
end

function _vac_half_integer(op::Op, spec::_VacSize, name::String, hint::String)
    s = _vac_size(op, spec, name, hint)
    (s isa Real && s >= 0 && isinteger(2s)) || throw(
        ArgumentError("`$name` must be a non-negative half-integer, got `$s`")
    )
    return Rational{Int}(Int(2s), 2)
end

function _vac_count(op::Op, spec::_VacSize, name::String, hint::String)
    n = _vac_size(op, spec, name, hint)
    (n isa Real && n >= 0 && isinteger(n)) || throw(
        ArgumentError("`$name` must be a non-negative integer, got `$n`")
    )
    return Int(n)
end
