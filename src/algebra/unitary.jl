const _SiteKey = Tuple{Int32, Index, Int32}

# The generators a transform covers on one site, taken from the rule keys rather than
# re-derived from a single `Op`. That is what lets a `PhaseSpace` pair and a
# `CollectiveTransition` work at all: their generating set is not recoverable from one
# member, but whoever built the rules was handed it.
struct SiteInfo
    key::_SiteKey
    generators::Vector{Op}
end

# An operator-valued displacement amplitude, with the sites it occupies precomputed. `U'XU`
# is `X` itself for every `X` commuting with the exponent, and something outside the class
# for every `X` that does not, so those sites are neither covered nor free: they are gated.
struct BlockedAmplitude
    amp::QAdd
    keys::Vector{_SiteKey}
end

"""
    UnitaryTransform

An exactly solvable unitary change of frame, stored as its action on operators.

The class covered is the one whose adjoint action closes on the generators of a site:
displacements, rotations and squeezes on a Fock or phase space, `SO(3)` rotations of a Pauli
or spin triple, and `U(n)` rotations of an N-level or collective transition set. Build one
with [`Displace`](@ref), [`Rotation`](@ref), [`Squeeze`](@ref), [`Bogoliubov`](@ref),
[`RotatingFrame`](@ref) or [`DressedFrame`](@ref), or from a generator with
`UnitaryTransform(G, θ)`; compose with `*`, invert with `inv`, and apply with
[`conjugate`](@ref) (observables) or [`transform`](@ref) (Hamiltonians).

Both directions are stored, so `inv` needs no symbolic reciprocal. The rule set must cover
every generator of every site it touches; [`conjugate`](@ref) throws otherwise rather than
return a half-transformed expression. A rule set keyed on an operator carrying a free index
covers that whole family, as the same local transform at every site of the index.

See also [`is_canonical`](@ref), [`gauge_term`](@ref), [`generators`](@ref),
[`constraints`](@ref).
"""
struct UnitaryTransform
    rules::Dict{Op, QAdd}
    inv_rules::Dict{Op, QAdd}
    sites::Vector{SiteInfo}
    gauge::QAdd
    time::Num
    reductions::Vector{ParamRelation}
    constraints::Vector{Num}
    blockers::Vector{BlockedAmplitude}
    wildcards::Vector{Index}

    function UnitaryTransform(
            rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, gauge::QAdd, time::Num,
            reductions::Vector{ParamRelation}, constraints::Vector{Num},
            blockers::Vector{BlockedAmplitude} = BlockedAmplitude[],
        )
        isempty(rules) && throw(ArgumentError("a `UnitaryTransform` needs at least one rule"))
        sites = _site_infos(rules)
        _check_complete(sites)
        _check_complete(_site_infos(inv_rules))
        wildcards = _wildcard_indices(rules)
        isempty(wildcards) || _check_families(rules, inv_rules, sites)
        for g in keys(rules)
            haskey(inv_rules, g) ||
                throw(ArgumentError("`inv_rules` is missing the generator `$g`"))
        end
        length(rules) == length(inv_rules) || throw(
            ArgumentError("`rules` and `inv_rules` must cover the same generators")
        )
        # A numeric parameter folds its relation's members to literals, which the reduction
        # must never see. Filtering here covers construction, composition and inversion.
        rels = all(_is_usable_rel, reductions) ? reductions :
            filter(_is_usable_rel, reductions)
        return new(
            rules, inv_rules, sites, gauge, time, rels, constraints, blockers, wildcards,
        )
    end

    # Trusted path, for rule sets whose completeness follows from an already-validated one:
    # `inv` swaps the two maps, `_merge_rules` unions two complete key sets, `_with_gauge`
    # reuses its input's verbatim. Re-deriving completeness there is pure cost, and every
    # time-dependent constructor rebuilds on top of a static one.
    global _unchecked_transform(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, gauge::QAdd, time::Num,
        reductions::Vector{ParamRelation}, constraints::Vector{Num},
        blockers::Vector{BlockedAmplitude} = BlockedAmplitude[],
    ) = new(
        rules, inv_rules, _site_infos(rules), gauge, time, reductions, constraints, blockers,
        _wildcard_indices(rules),
    )
end

const _site_key = site_key

function _site_infos(rules::Dict{Op, QAdd})
    out = SiteInfo[]
    for g in sort!(collect(keys(rules)))
        k = _site_key(g)
        i = findfirst(s -> s.key == k, out)
        i === nothing ? push!(out, SiteInfo(k, Op[g])) : push!(out[i].generators, g)
    end
    sort!(out; by = s -> (s.key[1], _index_key(s.key[2]), _name_rank(s.key[3])))
    return out
end

# === Indexed families ===
#
# A rule set keyed on an operator carrying a free index means the same local transform at
# every site of that index's range: `U = ⊗ᵢ Uᵢ`. Everything below follows from that sentence.
# An abstract index is the wildcard and a per-slot one (`i(3)`) stays an ordinary site, so
# both readings are available with no second spelling and no index minted on the user's
# behalf.

const _FamilyKey = Tuple{Int32, Int32}

_family_key(o::Op) = (o.space_index, _site_key(o)[3])
_is_wildcard(i::Index) = has_index(i) && index_slot(i) === nothing

function _wildcard_indices(rules::Dict{Op, QAdd})
    out = Index[]
    for g in keys(rules)
        _is_wildcard(g.index) && !any(i -> i == g.index, out) && push!(out, g.index)
    end
    return out
end

# The wildcard covering `o`'s family, or `nothing`.
function _wildcard_for(U::UnitaryTransform, o::Op)
    fk = _family_key(o)
    for s in U.sites
        _is_wildcard(s.key[2]) && (s.key[1], s.key[3]) == fk && return s.key[2]
    end
    return nothing
end

# A family rule image must stay on the site it acts on. Anything else is unsound rather than
# merely unsupported: `a_i -> cosh(r)*a_i + sinh(r)*b'` passes every same-site residual while
# `[ã_i, ã_j'] = -sinh(r)^2`, because the instances then share the off-family site. With the
# restriction the instances have disjoint support, so they commute for free and no
# cross-index residual is needed. It also keeps substitution clear of diagonal splitting:
# the multiset of `(space, index)` pairs per term is unchanged, so no same-space index pair
# is created that was not already there.
function _check_families(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, sites::Vector{SiteInfo},
    )
    for s in sites, u in sites
        s === u && continue
        ((s.key[1], s.key[3]) == (u.key[1], u.key[3])) || continue
        _is_wildcard(s.key[2]) || continue
        _is_wildcard(u.key[2]) && throw(
            ArgumentError(
                "this transform covers the family of `$(first(s.generators))` twice, " *
                    "under the indices `$(index_name(s.key[2]))` and " *
                    "`$(index_name(u.key[2]))`; rename one side to a single `Index`"
            )
        )
        throw(
            ArgumentError(
                "this transform covers `$(first(u.generators))` and the family of " *
                    "`$(first(s.generators))` at once; an operator family is either " *
                    "covered as a whole or site by site, not both"
            )
        )
    end
    for d in (rules, inv_rules), (g, r) in d
        _is_wildcard(g.index) || continue
        k = _site_key(g)
        for (term, _) in r, o in term.ops
            _site_key(o) == k || throw(
                ArgumentError(
                    "a rule of the indexed family `$g` maps it onto `$o` on another site. " *
                        "An indexed family is the same local transform at every site, so " *
                        "its rules must stay on the site they act on"
                )
            )
        end
    end
    return nothing
end

# `change_index` is a `*`-algebra isomorphism between the sites: operator relations are
# index-blind and renaming the coefficient index is a bijection on the parameters. The family
# map commutes with it by construction, which is why residuals at the wildcard certify every
# instance and only the rules have to be moved.
function _instantiate_rules(rules::Dict{Op, QAdd}, targets::Vector{Tuple{Index, Index}})
    out = copy(rules)
    for (w, j) in targets, (g, r) in rules
        (g.index == w && g.index != j) || continue
        out[_with_index(g, j)] = change_index(r, w, j)
    end
    return out
end

change_index(r::ParamRelation, from::Index, to::Index) =
    ParamRelation(change_index(r.hi, from, to), change_index(r.lo, from, to), r.sign)

# A relation on an `IndexedVariable` parameter holds per site, so it has to move with the
# rules. Without this a family squeeze declares `cosh(r(i))^2 - sinh(r(i))^2 = 1` and then
# false-negatives at every index but `i`.
function _instantiate_reductions(
        rels::Vector{ParamRelation}, targets::Vector{Tuple{Index, Index}},
    )
    isempty(rels) && return rels
    out = rels
    for (w, j) in targets
        w == j && continue
        out = _merge_reductions(out, ParamRelation[change_index(r, w, j) for r in rels])
    end
    return out
end

_site_info(U::UnitaryTransform, k::_SiteKey) = _site_info(U.sites, k)
function _site_info(sites::Vector{SiteInfo}, k::_SiteKey)
    for s in sites
        s.key == k && return s
    end
    return nothing
end

# Every operator generating the same site as `o`: the ladder pair for a Fock mode, the three
# axes of a Pauli or spin triple, all `n^2` matrix units of an N-level transition. Empty when
# the site cannot describe itself from one member, which the callers handle.
function _site_generators(o::Op)
    k = o.kind
    if k === OP_DESTROY || k === OP_CREATE
        return Op[
            Op(OP_DESTROY, o.name_id, o.space_index, o.index, 0, 0, 0, 0),
            Op(OP_CREATE, o.name_id, o.space_index, o.index, 0, 0, 0, 0),
        ]
    elseif k === OP_PAULI || k === OP_SPIN
        return Op[Op(k, o.name_id, o.space_index, o.index, Int32(ax), 0, 0, 0) for ax in 1:3]
    elseif k === OP_TRANSITION
        n = Int(o.nlev)
        return Op[
            Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)
                for i in 1:n for j in 1:n
        ]
    end
    # `Position`/`Momentum` and `CollectiveTransition` carry the rest of their generating set
    # off the `Op` (a second quadrature name, a level count), so it cannot be rebuilt from one
    # member. Those sites are validated against the rules the caller supplied instead.
    return _EMPTY_OPS
end

_is_phase_space(o::Op) = is_position(o) || is_momentum(o)
_is_fock(o::Op) = is_destroy(o) || is_create(o)

# The levels a collective rule set or generator list actually touches, sorted.
function _collective_levels(gens::Vector{Op})
    levels = Int[]
    for o in gens
        for l in (Int(o.l1), Int(o.l2))
            l in levels || push!(levels, l)
        end
    end
    return sort!(levels)
end

# A rule set covering part of a site would silently mix frames, so completeness is an
# invariant of the type rather than a check at the application site. Checked once per site
# against that site's own generators, not once per rule: for `n` levels the rule set has
# `n^2` keys and each `_site_generators` call builds `n^2` operators.
function _check_complete(sites::Vector{SiteInfo})
    for s in sites
        g = first(s.generators)
        expected = _site_generators(g)
        if isempty(expected)
            # Self-description is unavailable, so check the invariant that is.
            if is_collective_transition(g)
                # The level count is on the space, not the operator, so the rule keys are
                # the only witness: whatever levels they touch must form a full square. A
                # proper subset is admitted deliberately, for a block-diagonal unitary that
                # is the identity elsewhere; operators outside it are then simply uncovered.
                levels = _collective_levels(s.generators)
                for i in levels, j in levels
                    any(o -> o.l1 == i && o.l2 == j, s.generators) || throw(
                        ArgumentError(
                            "incomplete rule set: `$g` is covered on levels " *
                                "$(join(levels, ", ")) but `S^{$i$j}` is not; a collective " *
                                "rule set must cover the full square of the levels it touches"
                        )
                    )
                end
                continue
            end
            _is_phase_space(g) || continue
            (
                any(o -> o.kind === OP_POSITION, s.generators) &&
                    any(o -> o.kind === OP_MOMENTUM, s.generators)
            ) || throw(
                ArgumentError(
                    "incomplete rule set: `$g` has no rule for its conjugate variable"
                )
            )
            continue
        end
        for e in expected
            any(o -> isequal(o, e), s.generators) ||
                throw(ArgumentError("incomplete rule set: `$g` is covered but `$e` is not"))
        end
    end
    return nothing
end

_has_site(sites::Vector{SiteInfo}, k::_SiteKey) = any(s -> s.key == k, sites)

# A `PhaseSpace` is one degree of freedom, so a second differently-named `Position` on it is
# the same observable under another name and there is no rule to add for it. Saying "pass
# every generator" there would send the caller after something that cannot exist.
_uncovered_message(o::Op) = _is_phase_space(o) ?
    "`$o` is a second quadrature name on a `PhaseSpace` this transform already covers; a " *
    "`PhaseSpace` carries exactly one canonical pair, so use `⊗` for a second mode" :
    is_collective_transition(o) ?
    "`$o` acts on a site this transform covers but has no rule. A collective operator does " *
    "not carry the level count of its space, so the levels covered were inferred from the " *
    "generator; name this one there too" :
    "`$o` acts on a site this transform covers but has no rule; pass every generator of " *
    "the site to the constructor"

# Eager canonicalization stores `σ12*σ21` as the leaf `σ11`, so a rule set covering only
# the off-diagonals leaves it untouched and returns a wrong answer with no diagnostic.
# Re-walking per `conjugate` measured 1.8% of it, so the linear `_has_site` stays.
_restrict_to_site(t::QTerm, k::_SiteKey) = Op[o for o in t.ops if _site_key(o) == k]

# `U'XU = X` for every `X` commuting with the exponent, and for every `X` that does not the
# image is `X` times a displacement operator, which is no polynomial in the generators. The
# amplitude's own sites are therefore neither covered nor free, and letting them through
# untouched is the wrong answer with no diagnostic. Tested per sub-product rather than per
# operator: `a` and `a'` each fail against an `a'*a` amplitude while `a'*a` passes, and
# `a'*a` is the whole point. Memoized per call, so a repeated site pattern costs one hash.
function _check_blockers(q::QAdd, U::UnitaryTransform)
    memo = Dict{Tuple{Int, Vector{Op}}, Bool}()
    for (t, _) in q, (bi, b) in enumerate(U.blockers), k in b.keys
        sub = _restrict_to_site(t, k)
        isempty(sub) && continue
        ok = get!(memo, (bi, sub)) do
            iszero(commutator(_single_qadd(_CNUM_ONE, sub), b.amp))
        end
        ok || throw(
            ArgumentError(
                "`$(_single_qadd(_CNUM_ONE, sub))` acts on the site of the displacement " *
                    "amplitude `$(b.amp)` and does not commute with it, so its image is " *
                    "itself times a displacement operator, which is not a polynomial in " *
                    "the generators. Only combinations commuting with the amplitude can " *
                    "be transformed"
            )
        )
    end
    return nothing
end

# Validates, and collects the `(wildcard, target)` pairs the expression needs instantiated.
function _validate_coverage(q::QAdd, U::UnitaryTransform)
    isempty(U.blockers) || _check_blockers(q, U)
    targets = Tuple{Index, Index}[]
    for (t, _) in q
        for o in t.ops
            k = _site_key(o)
            if _has_site(U.sites, k)
                haskey(U.rules, o) || throw(ArgumentError(_uncovered_message(o)))
                continue
            end
            w = isempty(U.wildcards) ? nothing : _wildcard_for(U, o)
            if w !== nothing
                # The unindexed operator is a different operator from any member of the
                # family, not the family's representative, so it is not covered by one.
                has_index(o.index) || throw(
                    ArgumentError(
                        "`$o` carries no index, but this transform covers the indexed " *
                            "family `$(_with_index(o, w))`; they are different operators"
                    )
                )
                # Completeness at the wildcard transports to the instance by renaming, so a
                # missing member here is a missing member of the family.
                haskey(U.rules, _with_index(o, w)) ||
                    throw(ArgumentError(_uncovered_message(o)))
                any(p -> p == (w, o.index), targets) || push!(targets, (w, o.index))
                continue
            end
            # Rules key on the exact `Op`, index included, so any other member of a covered
            # `(space, name)` family (`a` vs `a_i`, or `a_i` vs `a_j`) misses every rule and
            # would come back untransformed. Reaching here with the space and name equal
            # means the index differs, since an exact match was handled above.
            for s in U.sites
                (s.key[1] == k[1] && s.key[3] == k[3]) && throw(
                    ArgumentError(
                        "`$o` belongs to an operator family this transform covers under a " *
                            "different index; a family is covered as a whole only when its " *
                            "rules are keyed on a free index"
                    )
                )
            end
        end
    end
    return targets
end

function _reduce_params(q::QAdd, rels::Vector{ParamRelation}, gated::Bool)
    # `scratch` is reused across the terms of one expression, so a frame conjugation
    # allocates no relation list per term.
    scratch = ParamRelation[]
    return _map_coefficients(c -> _reduce_all(c, rels, gated, scratch), q)
end

_apply_rules(q::QAdd, rules::Dict{Op, QAdd}) = _substitute_op_rules(q, rules)

# `c * q` for a coefficient already in `Coeff` form. `*` handles this too now that `Coeff` is
# a `_CoeffLike`, but this skips its `isone` probe on a value that is never one here.
function _scale_qadd(c::CNum, q::QAdd)
    _iszero_cnum(c) && return _zero_qadd()
    d = QTermDict()
    for (t, ct) in q
        _addto_key!(d, _copy_key(t), _mul_cnum(ct, c))
    end
    return QAdd(d, copy(q.indices))
end

"""
    conjugate(A, U::UnitaryTransform) -> QAdd

The frame change of an observable, `U' * A * U`.

For a Hamiltonian or any other generator of time evolution use [`transform`](@ref) instead:
a time-dependent frame adds a gauge term that `conjugate` deliberately omits. The two agree
exactly when `U` is static.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables α::Real;

julia> conjugate(a, Displace(a, α))
α + a
```

See also [`transform`](@ref), [`UnitaryTransform`](@ref).
"""
function conjugate(q::QAdd, U::UnitaryTransform)
    targets = _validate_coverage(q, U)
    isempty(targets) &&
        return _reduce_params(_apply_rules(q, U.rules), U.reductions, true)
    return _reduce_params(
        _apply_rules(q, _instantiate_rules(U.rules, targets)),
        _instantiate_reductions(U.reductions, targets), true,
    )
end
conjugate(o::QSym, U::UnitaryTransform) = conjugate(_single_qadd(_CNUM_ONE, Op[o]), U)

"""
    transform(H, U::UnitaryTransform) -> QAdd

The frame change of a Hamiltonian, `U' * H * U + gauge_term(U)`.

The gauge term `im * (∂ₜU') * U` is what keeps the transformed Hamiltonian the generator of
motion in the new frame; it vanishes for a static `U`, where this agrees with
[`conjugate`](@ref). Use `conjugate` for observables.

!!! note
    QuantumOptics and QuantumOpticsBase export a `transform` of their own (the numeric change
    of basis). A session that loads either one must call this one as
    `SecondQuantizedAlgebra.transform`.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables ω t;

julia> iszero(transform(ω * a' * a, RotatingFrame(ω * a' * a, t)))
true
```

See also [`conjugate`](@ref), [`gauge_term`](@ref), [`RotatingFrame`](@ref).
"""
transform(q::QAdd, U::UnitaryTransform) =
    iszero(U.gauge) ? conjugate(q, U) : conjugate(q, U) + U.gauge
transform(o::QSym, U::UnitaryTransform) = transform(_single_qadd(_CNUM_ONE, Op[o]), U)

"""
    gauge_term(U::UnitaryTransform) -> QAdd

The `im * (∂ₜU') * U` term [`transform`](@ref) adds and [`conjugate`](@ref) omits. Zero for
a static transform.
"""
gauge_term(U::UnitaryTransform) = U.gauge

"""
    generators(U::UnitaryTransform) -> Vector{Op}

The operators `U` has a rule for, in canonical order.
"""
generators(U::UnitaryTransform) = sort!(collect(keys(U.rules)))

"""
    constraints(U::UnitaryTransform) -> Vector{Num}

Expressions that must vanish for `U` to be unitary. [`Bogoliubov`](@ref) returns
`u^2 - v^2 - 1`; a parametrization that is unitary by construction returns nothing.
"""
constraints(U::UnitaryTransform) = U.constraints

function Base.inv(U::UnitaryTransform)
    # `gauge(U') = -U*gauge(U)*U'`, which collapses to `-gauge(U)` only when the two
    # commute; routing through the inverse rules is exact for every constructor.
    gauge = iszero(U.gauge) ? U.gauge :
        -_reduce_params(_apply_rules(U.gauge, U.inv_rules), U.reductions, true)
    # The blocked sites are a property of the exponent, which inversion only negates.
    return _unchecked_transform(
        copy(U.inv_rules), copy(U.rules), gauge, U.time,
        copy(U.reductions), copy(U.constraints), copy(U.blockers),
    )
end

Base.adjoint(U::UnitaryTransform) = inv(U)

function _merge_rules(outer::Dict{Op, QAdd}, inner::Dict{Op, QAdd})
    out = Dict{Op, QAdd}()
    for (g, r) in outer
        out[g] = _apply_rules(r, inner)
    end
    for (g, r) in inner
        haskey(out, g) || (out[g] = r)
    end
    return out
end

_merge_unique_num(a::Vector{Num}, b::Vector{Num}) =
    _merge_unique(a, b, (x, ys) -> any(y -> isequal(y, x), ys))

_merge_reductions(a::Vector{ParamRelation}, b::Vector{ParamRelation}) = _merge_unique(
    a, b, (r, qs) -> any(q -> isequal(q.hi, r.hi) && isequal(q.lo, r.lo), qs),
)

# A static factor was built with a zero gauge. If its rules move with the time variable it is
# about to adopt, that zero is wrong rather than merely absent, and the composed transform
# would generate the wrong motion with no diagnostic.
function _check_static_gauge(U::UnitaryTransform, t::Num)
    (_iszero_num(U.time) && !_iszero_num(t)) || return nothing
    for r in values(U.rules), (_, c) in r
        _depends_on_time(c, t) && throw(
            ArgumentError(
                "a static transform whose rules depend on `$t` cannot be composed with a " *
                    "time-dependent one: its gauge term was built as zero. Pass `$t` to " *
                    "its constructor instead"
            )
        )
    end
    return nothing
end

# A static transform adopts the other factor's time variable, so `Displace(a, α)` can be
# composed with a rotating frame.
function _common_time(U1::UnitaryTransform, U2::UnitaryTransform)
    _iszero_num(U1.time) && return U2.time
    _iszero_num(U2.time) && return U1.time
    isequal(U1.time, U2.time) && return U1.time
    return throw(
        ArgumentError(
            "cannot compose transforms with different time variables " *
                "(`$(U1.time)` and `$(U2.time)`)"
        )
    )
end

function Base.:*(U1::UnitaryTransform, U2::UnitaryTransform)
    time = _common_time(U1, U2)
    _check_static_gauge(U1, time)
    _check_static_gauge(U2, time)
    reductions = _merge_reductions(U1.reductions, U2.reductions)
    rules = _merge_rules(U1.rules, U2.rules)
    # `conjugate(g, inv(U1*U2)) = U1*(U2 g U2')*U1'`, so `U2` is the *inner* map here while
    # `U1` is the inner one for the forward rules. Writing both in the same order silently
    # breaks `inv` for non-commuting factors.
    inv_rules = _merge_rules(U2.inv_rules, U1.inv_rules)
    gauge = if iszero(U1.gauge)
        U2.gauge
    else
        _reduce_params(_apply_rules(U1.gauge, U2.rules), reductions, true) + U2.gauge
    end
    # A union of two complete key sets is complete on every site either one covered.
    return _unchecked_transform(
        rules, inv_rules, gauge, time, reductions,
        _merge_unique_num(U1.constraints, U2.constraints),
        _merge_unique(
            U1.blockers, U2.blockers, (b, bs) -> any(o -> isequal(o.amp, b.amp), bs),
        ),
    )
end

# === Rule construction helpers ===

# A rule value: a sum of `(coefficient, operators)` terms. Every rule built here is a linear
# combination of single operators plus a constant, so no canonicalization is needed.
function _rule_qadd(pairs::Vector{Tuple{CNum, Vector{Op}}})
    d = QTermDict()
    for (c, ops) in pairs
        _iszero_cnum(c) && continue
        _addto!(d, ops, c)
    end
    return QAdd(d, _EMPTY_INDICES)
end
# Splatting a runtime-length vector into the `Vararg` form forces a fresh specialization per
# length, so callers building a variable number of terms use the `Vector` method directly.
_rule_qadd(pairs::Vararg{Tuple{CNum, Vector{Op}}}) =
    _rule_qadd(Tuple{CNum, Vector{Op}}[pairs...])

_scaled(c::CNum, o::Op) = _rule_qadd((c, Op[o]))

# `U'a†U = (U'aU)†`, so the adjoint rule never has to be spelled out and cannot drift out of
# step with the forward one.
_with_adjoint(g::Op, r::QAdd) = Dict{Op, QAdd}(g => r, adjoint(g) => adjoint(r))

_lowering(a::Op) = is_create(a) ? adjoint(a) : a

_fock_or_throw(a::Op, what::AbstractString) = _is_fock(a) ? _lowering(a) :
    throw(ArgumentError("$what expects a Fock ladder operator, got $(a.kind)"))

# A unit phase is one atom, so conjugating it negates an exponent and `p*conj(p)` folds to 1
# with nothing to declare. Euler-splitting it into `cos`/`sin` would need a relation and
# still never close `exp(im*ϕ)*exp(-im*ϕ)`, since `cos(-ϕ)` is not folded to `cos(ϕ)`.
_phase(ϕ) = _phase_coeff(_as_num(ϕ))
_conj_phase(ϕ) = _conj_cnum(_phase(ϕ))

# A composite argument (`cos(ω*t)` in a rotating frame) is not an atom, so the coefficient
# never reaches the polynomial tier where the pair would be found automatically. Declaring it
# lets the reduction close the residual through the transient swap.
_trig_rel(θ) = ParamRelation(cos(θ), sin(θ), -1)
_hyp_rel(r) = ParamRelation(cosh(r), sinh(r), 1)

# `exp(y)` for a real `y`, oriented onto one atom per rate. `exp(-y)` is a *different* atom
# from `exp(y)` and the two never cancel, the same trap as a Euler-split phase; `1/exp(y)` is
# the one atom at exponent -1 and cancels natively.
function _real_scale(y::Num)
    e = expand(y)
    return _leading_sign(e) < 0 ? _to_cnum(1 / exp(expand(-e))) : _to_cnum(exp(e))
end

_static(
    rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd},
    rels::Vector{ParamRelation} = ParamRelation[],
    blockers::Vector{BlockedAmplitude} = BlockedAmplitude[],
) = UnitaryTransform(rules, inv_rules, _zero_qadd(), _NUM_ZERO, rels, Num[], blockers)

_with_gauge(U::UnitaryTransform, gauge::QAdd, t::Num) = _unchecked_transform(
    U.rules, U.inv_rules, _reduce_params(_family_gauge(U, gauge), U.reductions, true), t,
    U.reductions, U.constraints, U.blockers,
)

# `U = ⊗ᵢUᵢ` gives `im*(∂ₜU')U = Σᵢ im*(∂ₜUᵢ')Uᵢ`, but every timed constructor builds the
# single-site gauge. Summing it here covers all of them at once; `RotatingFrame` bypasses
# this, passing the already-summed `-H0`.
function _family_gauge(U::UnitaryTransform, gauge::QAdd)
    isempty(U.wildcards) && return gauge
    out = gauge
    for w in U.wildcards
        (_any_depends_on_index(out, w) && !any(i -> i == w, out.indices)) || continue
        iszero(index_range(w)) && throw(
            ArgumentError(
                "a time-dependent transform of the family indexed by " *
                    "`$(index_name(w))` needs a range on that index: its gauge term is " *
                    "the sum over the family"
            )
        )
        out = Σ(out, w)
    end
    return out
end

# === Time dependence ===

# `Differential` is defined on real `Num`s, and a complex amplitude arrives as a coefficient
# with real and imaginary parts, so differentiate componentwise.
_dt_num(x::Num, t::Num) = Symbolics.expand_derivatives(Symbolics.Differential(t)(x))
_dt(c::CNum, t::Num) = _cnum(_dt_num(real(c), t), _dt_num(imag(c), t))
_dt(x::_ScalarLike, t::Num) = _dt(_to_cnum(x), t)

function _depends_on_time(x::_CoeffLike, t::Num)
    tv = SymbolicUtils.unwrap(t)
    c = _to_cnum(x)
    for part in (real(c), imag(c))
        for v in Symbolics.get_variables(part)
            isequal(v, tv) && return true
        end
    end
    return false
end

# `im*(∂ₜU')U = -∂ₜθ*G` for `U = exp(-im*θ*G)`, exactly when `[G, ∂ₜG] = 0`. Every generator
# reaching this is time-independent, so the condition holds; the constructors whose generator
# does move in time (`Displace`, a phase-driven `Squeeze`) carry their own closed form.
_gauge(G::QAdd, θ::_ScalarLike, t::Num) = _scale_qadd(_neg_cnum(_dt(θ, t)), G)

# === Fock constructors ===

"""
    Displace(a::Op, α[, t]) -> UnitaryTransform

Displacement of a Fock mode, `D(α) = exp(α*a' - conj(α)*a)`, acting as `a -> a + α`.

`α` may be real or complex, numeric or symbolic; declare a complex amplitude as
`@variables α::Number`. Passing a time variable `t` makes the transform time-dependent, so
[`transform`](@ref) adds the corresponding gauge term.

!!! warning
    A time-dependent `α` needs `t` passed as well. Without it the transform is static by
    construction, its gauge term is zero, and [`transform`](@ref) agrees with
    [`conjugate`](@ref) even though the frame is moving. A static transform carries no time
    variable, so this cannot be detected until the transform is composed with one that does.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables α::Real;

julia> conjugate(a' * a, Displace(a, α))
α^2 + α * a + α * a' + a' * a
```

See also [`Rotation`](@ref), [`Squeeze`](@ref), [`conjugate`](@ref).
"""
function Displace(a::Op, α::_ScalarLike)
    d = _fock_or_throw(a, "`Displace`")
    c = _to_cnum(α)
    return _static(
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (c, Op[]))),
        _with_adjoint(d, _rule_qadd((_CNUM_ONE, Op[d]), (_neg_cnum(c), Op[]))),
    )
end

_rule_plus(g::Op, A::QAdd) = _single_qadd(_CNUM_ONE, Op[g]) + A

_site_keys(A::QAdd) = unique!(_SiteKey[_site_key(o) for (t, _) in A for o in t.ops])

# `U'gU = g + A` truncates at first order only when `A` commutes with the site it displaces,
# and is a unitary map only when `A` is normal: `[A, A']` would otherwise survive in the
# commutator of the images. Both are silent failures, so both are constructor invariants.
function _blocked_amplitude(A::QAdd, gens::Vector{Op}, what::AbstractString)
    for g in gens
        iszero(commutator(A, g)) || throw(
            ArgumentError(
                "$what needs an amplitude commuting with the site it displaces, but " *
                    "`[$A, $g]` is not zero"
            )
        )
    end
    iszero(commutator(A, adjoint(A))) || throw(
        ArgumentError(
            "$what needs an amplitude commuting with its own adjoint, but `[$A, $A']` is " *
                "not zero; the displacement it generates is not unitary"
        )
    )
    keys = _site_keys(A)
    for g in gens
        _site_key(g) in keys && throw(
            ArgumentError(
                "$what cannot displace `$g` by an amplitude acting on its own site"
            )
        )
    end
    return BlockedAmplitude(A, keys)
end

"""
    Displace(a::Op, A::QAdd[, t]) -> UnitaryTransform

Displacement of a Fock mode by an operator-valued amplitude, `exp(A*a' - A'*a)`, acting as
`a -> a + A`.

`A` must commute with `a` and with its own adjoint, which makes the Hadamard series truncate
at first order exactly as it does for a number. This is the polaron (Lang-Firsov) class: the
amplitude is an operator on another site that the displacement conserves.

Operators on `A`'s own site are neither transformed nor left alone. Those commuting with `A`
are invariant and pass through; the rest have images that are not polynomials in the
generators (`X` times a displacement operator), and [`conjugate`](@ref) throws for them
rather than return them unchanged.

# Examples

```jldoctest
julia> h = FockSpace(:c) ⊗ FockSpace(:m);

julia> a = Destroy(h, :a, 1); b = Destroy(h, :b, 2);

julia> @variables g ωm;

julia> conjugate(b, Displace(b, (-g / ωm) * (a' * a)))
b - (g / ωm) * a' * a
```

See also [`Displace`](@ref), [`conjugate`](@ref), [`is_canonical`](@ref).
"""
function Displace(a::Op, A::QAdd)
    d = _fock_or_throw(a, "`Displace`")
    blocked = _blocked_amplitude(A, Op[d, adjoint(d)], "`Displace`")
    return _static(
        _with_adjoint(d, _rule_plus(d, A)),
        _with_adjoint(d, _rule_plus(d, -A)),
        ParamRelation[], BlockedAmplitude[blocked],
    )
end

_dt_qadd(q::QAdd, t::Num) = _map_coefficients(c -> _dt(c, t), q)

function Displace(a::Op, A::QAdd, t::_ScalarLike)
    d = _fock_or_throw(a, "`Displace`")
    tt = _as_num(t)
    Ad = _dt_qadd(A, tt)
    # The scalar gauge below carries `[A, Ȧ]/2`, which for an operator amplitude is an
    # operator rather than a global phase. It is still the whole story only when `A` and `Ȧ`
    # commute; otherwise `e^{-A}∂ₜe^{A}` does not terminate at that order.
    iszero(commutator(A, Ad)) || throw(
        ArgumentError(
            "`Displace` needs an amplitude commuting with its own time derivative, but " *
                "`[$A, $Ad]` is not zero"
        )
    )
    gauge = -im * (Ad * adjoint(d) - adjoint(Ad) * d) -
        (im // 2) * (adjoint(A) * Ad - A * adjoint(Ad))
    return _with_gauge(Displace(a, A), gauge, tt)
end

function Displace(a::Op, α::_ScalarLike, t::_ScalarLike)
    d = _fock_or_throw(a, "`Displace`")
    tt = _as_num(t)
    c = _to_cnum(α)
    cd = _dt(c, tt)
    # `[A, Ȧ]` is a nonzero c-number here, so the gauge is not `-∂ₜθ*G`:
    # `D'∂ₜD = Ȧ - [A, Ȧ]/2` with `A = α*a' - conj(α)*a`.
    cnum = _mul_cnum(
        _CNUM_NEG_IM,
        _mul_cnum(
            _CNUM_HALF,
            _add_cnum(
                _mul_cnum(_conj_cnum(c), cd), _neg_cnum(_mul_cnum(c, _conj_cnum(cd)))
            ),
        ),
    )
    gauge = _rule_qadd(
        (_mul_cnum(_CNUM_NEG_IM, cd), Op[adjoint(d)]),
        (_mul_cnum(_CNUM_IM, _conj_cnum(cd)), Op[d]),
        (cnum, Op[]),
    )
    return _with_gauge(Displace(a, α), gauge, tt)
end

"""
    Rotation(a::Op, θ[, t]) -> UnitaryTransform

Phase rotation of a Fock mode, `exp(-im*θ*a'*a)`, acting as `a -> exp(-im*θ)*a`.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables θ;

julia> conjugate(a' * a, Rotation(a, θ))
a' * a
```

See also [`Displace`](@ref), [`Squeeze`](@ref), [`RotatingFrame`](@ref).
"""
function Rotation(a::Op, θ::_ScalarLike)
    d = _fock_or_throw(a, "`Rotation`")
    return _static(
        _with_adjoint(d, _scaled(_conj_phase(θ), d)),
        _with_adjoint(d, _scaled(_phase(θ), d)),
    )
end

function _timed_fock_rotation(a::Op, θ::_ScalarLike, t::_ScalarLike)
    d = _fock_or_throw(a, "`Rotation`")
    tt = _as_num(t)
    return _with_gauge(Rotation(a, θ), _gauge(adjoint(d) * d, θ, tt), tt)
end

"""
    Squeeze(a::Op, r, ϕ = 0[, t]) -> UnitaryTransform

Single-mode squeezing, `exp((z*a'^2 - conj(z)*a^2)/2)` with `z = r*exp(im*ϕ)`, acting as
`a -> cosh(r)*a + exp(im*ϕ)*sinh(r)*a'`.

`r` and `ϕ` may both depend on `t`; a moving phase adds `∂ₜϕ*sinh(r)^2*(a'*a + 1/2)` to the
gauge term on top of the usual `-∂ₜr*G`.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables r;

julia> conjugate(a, Squeeze(a, r))
cosh(r) * a + sinh(r) * a'
```

See also [`Displace`](@ref), [`Rotation`](@ref), [`is_canonical`](@ref).
"""
function Squeeze(a::Op, r::_ScalarLike, ϕ::_ScalarLike = 0)
    d = _fock_or_throw(a, "`Squeeze`")
    ch = _to_cnum(cosh(r))
    sh = _mul_cnum(_phase(ϕ), _to_cnum(sinh(r)))
    return _static(
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (sh, Op[adjoint(d)]))),
        _with_adjoint(d, _rule_qadd((ch, Op[d]), (_neg_cnum(sh), Op[adjoint(d)]))),
        ParamRelation[_hyp_rel(r)],
    )
end

# `[A, ∂ₜA]` is not zero once `ϕ` moves, so `-∂ₜr*G` is not the gauge. Solving
# `[K, g̃] = im*∂ₜg̃` on the Bogoliubov coefficients `(cosh r, exp(im*ϕ)*sinh r)` fixes the
# operator part, and su(1,1) closes on `a'*a + 1/2` rather than `a'*a`, which is where the
# c-number comes from. Both reduce to `-∂ₜr*G` at constant `ϕ`.
function _squeeze_gauge(d::Op, r::_ScalarLike, ϕ::_ScalarLike, tt::Num)
    sh = _to_cnum(sinh(r))
    ch = _to_cnum(cosh(r))
    idr = _mul_cnum(_CNUM_IM, _dt(r, tt))
    dϕ = _dt(ϕ, tt)
    w = _mul_cnum(dϕ, _mul_cnum(sh, sh))
    m = _mul_cnum(dϕ, _mul_cnum(sh, ch))
    return _rule_qadd(
        (w, Op[adjoint(d), d]),
        (_mul_cnum(_CNUM_HALF, w), Op[]),
        (
            _mul_cnum(_CNUM_HALF, _mul_cnum(_phase(ϕ), _add_cnum(m, _neg_cnum(idr)))),
            Op[adjoint(d), adjoint(d)],
        ),
        (_mul_cnum(_CNUM_HALF, _mul_cnum(_conj_phase(ϕ), _add_cnum(m, idr))), Op[d, d]),
    )
end

function Squeeze(a::Op, r::_ScalarLike, ϕ::_ScalarLike, t::_ScalarLike)
    d = _fock_or_throw(a, "`Squeeze`")
    tt = _as_num(t)
    return _with_gauge(Squeeze(a, r, ϕ), _squeeze_gauge(d, r, ϕ, tt), tt)
end

# === Two-mode and phase-space constructors ===

function _two_modes(a::Op, b::Op, what::AbstractString)
    x = _fock_or_throw(a, what)
    y = _fock_or_throw(b, what)
    _site_key(x) == _site_key(y) &&
        throw(ArgumentError("$what needs two distinct modes, got `$x` twice"))
    return (x, y)
end

# `Position` and `Momentum` carry different names and hold no reference to their Hilbert
# space, so the conjugate variable has to be supplied explicitly.
function _phase_pair(x::Op, p::Op, what::AbstractString)
    (x.kind === OP_POSITION && p.kind === OP_MOMENTUM) ||
        throw(ArgumentError("$what expects a `(Position, Momentum)` pair in that order"))
    (x.space_index == p.space_index && x.index == p.index) ||
        throw(ArgumentError("$what expects both operators on the same site"))
    return nothing
end

_pair_rules(x::Op, p::Op, fx::QAdd, fp::QAdd) = Dict{Op, QAdd}(x => fx, p => fp)

# Each constructor spells its inverse out instead of deriving it, so that its sign convention
# reads at the call site; `canonicality_residuals` round-trips both directions. See devdocs.
function _beamsplitter(a::Op, b::Op, θ::_ScalarLike)
    x, y = _two_modes(a, b, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    ns = _neg_cnum(s)
    return _static(
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (s, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (ns, Op[x]))),
        ),
        merge(
            _with_adjoint(x, _rule_qadd((c, Op[x]), (ns, Op[y]))),
            _with_adjoint(y, _rule_qadd((c, Op[y]), (s, Op[x]))),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function _quad_rotation(x::Op, p::Op, θ::_ScalarLike)
    _phase_pair(x, p, "`Rotation`")
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    ns = _neg_cnum(s)
    return _static(
        _pair_rules(x, p, _rule_qadd((c, Op[x]), (s, Op[p])), _rule_qadd((c, Op[p]), (ns, Op[x]))),
        _pair_rules(x, p, _rule_qadd((c, Op[x]), (ns, Op[p])), _rule_qadd((c, Op[p]), (s, Op[x]))),
        ParamRelation[_trig_rel(θ)],
    )
end

"""
    Rotation(a::Op, b::Op, θ[, t]) -> UnitaryTransform

Rotation mixing two operators.

For two Fock modes this is the beamsplitter `exp(-im*θ*im*(a'*b - b'*a))`, acting as
`a -> cos(θ)*a + sin(θ)*b`. For a `(Position, Momentum)` pair it is the quadrature rotation
`exp(-im*θ*(x^2 + p^2)/2)`, acting as `x -> cos(θ)*x + sin(θ)*p`, which matches
[`Rotation`](@ref) under `a = (x + im*p)/sqrt(2)`.

# Examples

```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> a = Destroy(h, :a, 1); b = Destroy(h, :b, 2);

julia> @variables θ;

julia> conjugate(a, Rotation(a, b, θ))
cos(θ) * a + sin(θ) * b
```

See also [`Squeeze`](@ref), [`is_canonical`](@ref).
"""
Rotation(a::Op, b::Op, θ::_ScalarLike) =
    _is_phase_space(a) ? _quad_rotation(a, b, θ) : _beamsplitter(a, b, θ)

function Rotation(a::Op, b::Op, θ::_ScalarLike, t::_ScalarLike)
    tt = _as_num(t)
    U = Rotation(a, b, θ)
    G = if _is_phase_space(a)
        (a * a + b * b) * (1 // 2)
    else
        im * (adjoint(a) * b - adjoint(b) * a)
    end
    return _with_gauge(U, _gauge(G, θ, tt), tt)
end

function _two_mode_squeeze(
        a::Op, b::Op, u::CNum, v::CNum, rels::Vector{ParamRelation}, cons::Vector{Num},
        what::AbstractString,
    )
    x, y = _two_modes(a, b, what)
    nv = _neg_cnum(v)
    rules = merge(
        _with_adjoint(x, _rule_qadd((u, Op[x]), (v, Op[adjoint(y)]))),
        _with_adjoint(y, _rule_qadd((u, Op[y]), (v, Op[adjoint(x)]))),
    )
    inv_rules = merge(
        _with_adjoint(x, _rule_qadd((u, Op[x]), (nv, Op[adjoint(y)]))),
        _with_adjoint(y, _rule_qadd((u, Op[y]), (nv, Op[adjoint(x)]))),
    )
    return UnitaryTransform(rules, inv_rules, _zero_qadd(), _NUM_ZERO, rels, cons)
end

function _quad_squeeze(x::Op, p::Op, r::_ScalarLike)
    _phase_pair(x, p, "`Squeeze`")
    # `exp(-r)` has a composite argument, so it is not an atom and would never cancel
    # against `exp(r)`. `1/exp(r)` recognizes as the same atom at exponent `-1`, so
    # `_merge_factors` cancels the pair natively and no relation is needed.
    up = _to_cnum(exp(r))
    dn = _to_cnum(1 / exp(r))
    return _static(
        _pair_rules(x, p, _scaled(up, x), _scaled(dn, p)),
        _pair_rules(x, p, _scaled(dn, x), _scaled(up, p)),
        ParamRelation[],
    )
end

"""
    Squeeze(a::Op, b::Op, r[, t]) -> UnitaryTransform

Squeezing of two operators.

For two Fock modes this is the two-mode squeeze `exp(r*(a'*b' - b*a))`, acting as
`a -> cosh(r)*a + sinh(r)*b'`. For a `(Position, Momentum)` pair it is
`exp(-im*r*(x*p + p*x)/2)`, acting as `x -> exp(r)*x` and `p -> exp(-r)*p`, which matches
[`Squeeze`](@ref) at `ϕ = 0`.

# Examples

```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> a = Destroy(h, :a, 1); b = Destroy(h, :b, 2);

julia> @variables r;

julia> conjugate(a, Squeeze(a, b, r))
cosh(r) * a + sinh(r) * b'
```

See also [`Bogoliubov`](@ref), [`Rotation`](@ref).
"""
function Squeeze(a::Op, b::Op, r::_ScalarLike)
    _is_phase_space(a) && return _quad_squeeze(a, b, r)
    return _two_mode_squeeze(
        a, b, _to_cnum(cosh(r)), _to_cnum(sinh(r)), ParamRelation[_hyp_rel(r)], Num[],
        "`Squeeze`",
    )
end

function Squeeze(a::Op, b::Op, r::_ScalarLike, t::_ScalarLike)
    tt = _as_num(t)
    U = Squeeze(a, b, r)
    G = if _is_phase_space(a)
        # (x*p + p*x)/2
        (a * b + b * a) * (1 // 2)
    else
        im * (adjoint(a) * adjoint(b) - b * a)
    end
    return _with_gauge(U, _gauge(G, r, tt), tt)
end

"""
    Bogoliubov(a::Op, b::Op, u, v) -> UnitaryTransform

Two-mode Bogoliubov transformation in raw coefficients, `a -> u*a + v*b'`.

Unitarity is not automatic: `u^2 - v^2 == 1` is what [`constraints`](@ref) returns and what
the coefficient reduction uses to close residuals. Prefer
[`Squeeze`](@ref) when the hyperbolic parametrization `u = cosh(r)`,
`v = sinh(r)` is acceptable.

# Examples

```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> a = Destroy(h, :a, 1); b = Destroy(h, :b, 2);

julia> @variables u v;

julia> is_canonical(Bogoliubov(a, b, u, v))
true
```

See also [`Squeeze`](@ref), [`constraints`](@ref).
"""
Bogoliubov(a::Op, b::Op, u::_ScalarLike, v::_ScalarLike) = _two_mode_squeeze(
    a, b, _to_cnum(u), _to_cnum(v),
    ParamRelation[ParamRelation(u, v, 1)],
    Num[_as_num(u)^2 - _as_num(v)^2 - 1],
    "`Bogoliubov`",
)

"""
    Displace(x::Op, p::Op, dx, dp[, t]) -> UnitaryTransform

Phase-space displacement, `exp(im*(dp*x - dx*p))`, acting as `x -> x + dx` and
`p -> p + dp`.

# Examples

```jldoctest
julia> h = PhaseSpace(:q); x = Position(h, :x); p = Momentum(h, :p);

julia> @variables dx dp;

julia> conjugate(p, Displace(x, p, dx, dp))
dp + p
```

See also [`Displace`](@ref), [`Rotation`](@ref).
"""
function Displace(x::Op, p::Op, dx::_ScalarLike, dp::_ScalarLike)
    _phase_pair(x, p, "`Displace`")
    cx = _to_cnum(dx)
    cp = _to_cnum(dp)
    return _static(
        _pair_rules(
            x, p,
            _rule_qadd((_CNUM_ONE, Op[x]), (cx, Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (cp, Op[])),
        ),
        _pair_rules(
            x, p,
            _rule_qadd((_CNUM_ONE, Op[x]), (_neg_cnum(cx), Op[])),
            _rule_qadd((_CNUM_ONE, Op[p]), (_neg_cnum(cp), Op[])),
        ),
    )
end

"""
    Displace(x::Op, p::Op, DX::QAdd, DP::QAdd) -> UnitaryTransform

Phase-space displacement by operator-valued amplitudes, acting as `x -> x + DX` and
`p -> p + DP`.

Both amplitudes must be Hermitian, commute with the pair they displace, and commute with each
other. As in the Fock case, operators on their sites are transformed only where they commute
with the amplitude, and [`conjugate`](@ref) throws otherwise.

See also [`Displace`](@ref).
"""
function Displace(x::Op, p::Op, DX::QAdd, DP::QAdd)
    _phase_pair(x, p, "`Displace`")
    gens = Op[x, p]
    # `x` and `p` are Hermitian and a displacement must keep them so; the CCR alone does not
    # pin that, since `x -> x + im` preserves it and no unitary induces it.
    for A in (DX, DP)
        isequal(adjoint(A), A) || throw(
            ArgumentError(
                "`Displace` needs Hermitian amplitudes on a phase-space pair, but `$A` is " *
                    "not; the displaced quadrature would not be an observable"
            )
        )
    end
    iszero(commutator(DX, DP)) || throw(
        ArgumentError("`Displace` needs commuting amplitudes, but `[$DX, $DP]` is not zero")
    )
    bx = _blocked_amplitude(DX, gens, "`Displace`")
    bp = _blocked_amplitude(DP, gens, "`Displace`")
    return _static(
        _pair_rules(x, p, _rule_plus(x, DX), _rule_plus(p, DP)),
        _pair_rules(x, p, _rule_plus(x, -DX), _rule_plus(p, -DP)),
        ParamRelation[], BlockedAmplitude[bx, bp],
    )
end

function Displace(x::Op, p::Op, dx::_ScalarLike, dp::_ScalarLike, t::_ScalarLike)
    _phase_pair(x, p, "`Displace`")
    tt = _as_num(t)
    cx = _to_cnum(dx)
    cp = _to_cnum(dp)
    dxd = _dt(cx, tt)
    dpd = _dt(cp, tt)
    # Same defect as the Fock displacement: `[A, Ȧ] = im*(dx*∂ₜdp - dp*∂ₜdx)` is a nonzero
    # c-number, so the gauge carries a scalar beyond `-∂ₜθ*G`.
    cnum = _mul_cnum(
        _neg_cnum(_CNUM_HALF),
        _add_cnum(_mul_cnum(cp, dxd), _neg_cnum(_mul_cnum(cx, dpd))),
    )
    gauge = _rule_qadd((dpd, Op[x]), (_neg_cnum(dxd), Op[p]), (cnum, Op[]))
    return _with_gauge(Displace(x, p, dx, dp), gauge, tt)
end

# === Pauli / spin constructors ===

_axis_or_throw(ax::Integer) =
    1 <= ax <= 3 ? Int32(ax) : throw(ArgumentError("axis must be 1, 2 or 3, got $ax"))
# Reached when a Pauli/Spin rotation is given a symbolic or fractional second argument,
# which for a triple can only have been meant as an axis.
_axis_or_throw(ax) = throw(ArgumentError("axis must be 1, 2 or 3, got $ax"))

_axis_op(o::Op, ax::Integer) = Op(o.kind, o.name_id, o.space_index, o.index, Int32(ax), 0, 0, 0)

_triple_or_throw(S::Op, what::AbstractString) =
    (S.kind === OP_PAULI || S.kind === OP_SPIN) ? nothing :
    throw(ArgumentError("$what expects a Pauli or Spin operator, got $(S.kind)"))

"""
    Rotation(S::Op, axis::Integer, θ[, t]) -> UnitaryTransform

Rotation of a Pauli or spin triple by `θ` about `axis` (`1`, `2` or `3`).

This is `exp(-im*θ*S_axis)` for a [`Spin`](@ref) and `exp(-im*θ*σ_axis/2)` for a
[`Pauli`](@ref), so that in both cases the triple rotates by the angle `θ`.

# Examples

```jldoctest
julia> h = SpinSpace(:S); Sx = Spin(h, :S, 1);

julia> @variables θ;

julia> conjugate(Sx, Rotation(Sx, 3, θ))
cos(θ) * Sx - sin(θ) * Sy
```

See also [`is_canonical`](@ref), [`RotatingFrame`](@ref).
"""
function _axis_rotation(S::Op, axis::_ScalarLike, θ::_ScalarLike)
    _triple_or_throw(S, "`Rotation`")
    m = Int(_axis_or_throw(axis))
    ou, ov, om = _axis_op(S, mod1(m + 1, 3)), _axis_op(S, mod1(m + 2, 3)), _axis_op(S, m)
    c = _to_cnum(cos(θ))
    s = _to_cnum(sin(θ))
    ns = _neg_cnum(s)
    return _static(
        Dict{Op, QAdd}(
            ou => _rule_qadd((c, Op[ou]), (ns, Op[ov])),
            ov => _rule_qadd((c, Op[ov]), (s, Op[ou])),
            om => _scaled(_CNUM_ONE, om),
        ),
        Dict{Op, QAdd}(
            ou => _rule_qadd((c, Op[ou]), (s, Op[ov])),
            ov => _rule_qadd((c, Op[ov]), (ns, Op[ou])),
            om => _scaled(_CNUM_ONE, om),
        ),
        ParamRelation[_trig_rel(θ)],
    )
end

function Rotation(S::Op, axis::Integer, θ::_ScalarLike, t::_ScalarLike)
    _triple_or_throw(S, "`Rotation`")
    tt = _as_num(t)
    om = _axis_op(S, Int(_axis_or_throw(axis)))
    # The Pauli generator carries the factor 1/2 that makes `θ` the rotation angle.
    G = S.kind === OP_PAULI ? _scaled(_CNUM_HALF, om) : _scaled(_CNUM_ONE, om)
    return _with_gauge(_axis_rotation(S, axis, θ), _gauge(G, θ, tt), tt)
end

# `Rotation(S, axis, θ)` and `Rotation(a, θ, t)` are both `(Op, scalar, scalar)`. Dispatching
# on `axis::Integer` picked the triple method whenever the angle happened to be written as an
# integer literal, so `Rotation(a, 2, t)` and `Rotation(a, 2.0, t)` built different things.
# The operator is what actually distinguishes the two, so branch on that.
function Rotation(o::Op, x::_ScalarLike, y::_ScalarLike)
    (o.kind === OP_PAULI || o.kind === OP_SPIN) && return _axis_rotation(o, x, y)
    return _timed_fock_rotation(o, x, y)
end

# === N-level constructor ===

_transition_op(o::Op, i::Integer, j::Integer) =
    Op(OP_TRANSITION, o.name_id, o.space_index, o.index, Int32(i), Int32(j), o.g, o.nlev)

# A collective `Op` carries no level count (`nlev == 0`): it lives on the space, which the
# operator does not reference. The builder below therefore takes `n` from the caller.
_collective_op(o::Op, i::Integer, j::Integer) =
    Op(OP_COLLECTIVE_TRANSITION, o.name_id, o.space_index, NO_INDEX, Int32(i), Int32(j), 0, 0)

# `U'|i⟩⟨j|U = Σ_kl conj(W[i,k])*W[j,l]*|k⟩⟨l|`: the contraction runs over W's *row* index.
# It holds verbatim for a collective `S^{ij} = Σ_a |i⟩_a⟨j|` under `⊗_a U`, because
# conjugation is linear in the atom sum. Only linearity is used, never the matrix-unit
# product law, which is why the collective failure of `S^{ij}S^{kl} = δ_{jk}S^{il}` is
# irrelevant here and matters only to the residuals.
function _matrix_unit_rules(σ::Op, W::AbstractMatrix, n::Int, mk::F) where {F}
    out = Dict{Op, QAdd}()
    # Hoisted out of the innermost loop: each entry was otherwise recognized `n^2` times, and
    # for a symbolic `W` every `_to_cnum` walks a SymbolicUtils tree.
    Wc = CNum[_to_cnum(W[i, j]) for i in 1:n, j in 1:n]
    Wk = CNum[_conj_cnum(Wc[i, j]) for i in 1:n, j in 1:n]
    ops = Op[mk(k, l) for k in 1:n, l in 1:n]
    for i in 1:n, j in 1:n
        d = QTermDict()
        for k in 1:n, l in 1:n
            c = _fold_roots(_mul_cnum(Wk[i, k], Wc[j, l]))
            _iszero_cnum(c) && continue
            _addto!(d, Op[ops[k, l]], c)
        end
        out[mk(i, j)] = QAdd(d, _EMPTY_INDICES)
    end
    return out
end

_nlevel_rules(σ::Op, W::AbstractMatrix) =
    _matrix_unit_rules(σ, W, Int(σ.nlev), (i, j) -> _transition_op(σ, i, j))
_collective_rules(σ::Op, W::AbstractMatrix) =
    _matrix_unit_rules(σ, W, size(W, 1), (i, j) -> _collective_op(σ, i, j))

# `W†[i,j] = conj(W[j,i])`, conjugated through the coefficient so a `Number`-symtype entry
# gets the same `conj` wrapper the reduction knows how to see through.
_dagger(W::AbstractMatrix) =
    [to_num(_conj_cnum(_to_cnum(W[j, i]))) for i in axes(W, 2), j in axes(W, 1)]

function _nlevel_or_throw(σ::Op, W::AbstractMatrix, what::AbstractString)
    σ.kind === OP_TRANSITION ||
        throw(ArgumentError("$what expects a `Transition` operator, got $(σ.kind)"))
    n = Int(σ.nlev)
    size(W) == (n, n) ||
        throw(ArgumentError("`W` must be $n×$n for a $n-level space, got $(size(W))"))
    return nothing
end

"""
    Rotation(σ::Op, W::AbstractMatrix) -> UnitaryTransform
    Rotation(σ::Op, W::AbstractMatrix, t) -> UnitaryTransform
    Rotation(σ::Op, W::AbstractMatrix, t, gauge::QAdd) -> UnitaryTransform

Basis rotation of an N-level system by the unitary matrix `W`, acting as
`σ^{ij} -> sum_kl conj(W[i,k])*W[j,l]*σ^{kl}`.

A [`CollectiveTransition`](@ref) is accepted too, under the same contraction: it is linear in
the atom sum, so it never uses the matrix-unit product law that collective transitions do not
obey. There the level count comes from `size(W)`, since a collective operator does not carry
one.

`W` is the matrix of the unitary itself, `W[k,l] = ⟨k|U|l⟩`. Unitarity of `W` is not
checked on construction; [`is_canonical`](@ref) tests it through the transformed algebra.

With a time variable `t`, the gauge term is `im*Ẇ'W`, obtained by differentiating `W`
entrywise. The four-argument form replaces that with the `gauge` given, for a `W` whose time
dependence is not recoverable by entrywise differentiation. It is not checked: an incorrect
`gauge` makes [`transform`](@ref) return something that is not the generator of motion in the
new frame, and [`is_canonical`](@ref) does not test the gauge.

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, 2); σ = Transition(h, :σ, 1, 2);

julia> conjugate(σ, Rotation(σ, [0 1; 1 0]))
σ₂₁
```

See also [`RotatingFrame`](@ref), [`is_canonical`](@ref).
"""
function Rotation(σ::Op, W::AbstractMatrix)
    is_collective_transition(σ) && return _collective_rotation(σ, W)
    _nlevel_or_throw(σ, W, "`Rotation`")
    return _static(_nlevel_rules(σ, W), _nlevel_rules(σ, _dagger(W)))
end

function _collective_or_throw(σ::Op, W::AbstractMatrix, what::AbstractString)
    size(W, 1) == size(W, 2) ||
        throw(ArgumentError("$what needs a square `W` on a collective site, got $(size(W))"))
    n = size(W, 1)
    (1 <= σ.l1 <= n && 1 <= σ.l2 <= n) || throw(
        ArgumentError("`$σ` uses a level outside `1:$n`, the size of `W`")
    )
    return n
end

function _collective_rotation(σ::Op, W::AbstractMatrix)
    _collective_or_throw(σ, W, "`Rotation`")
    return _static(_collective_rules(σ, W), _collective_rules(σ, _dagger(W)))
end

# `U = Σ_ij W[i,j]*σ^{ij}`, so `im*(∂ₜU')U` is the operator of the matrix `im*Ẇ'W`. That is
# Hermitian for free: `Ẇ'W + W'Ẇ = ∂ₜ(W'W) = 0`. A statement about `W` alone, so the
# collective site shares it and only the operator builder differs.
function _matrix_unit_gauge(W::AbstractMatrix, n::Int, mk::F, tt::Num) where {F}
    Wc = CNum[_to_cnum(W[k, l]) for k in 1:n, l in 1:n]
    Wd = CNum[_conj_cnum(_dt(Wc[k, l], tt)) for k in 1:n, l in 1:n]
    d = QTermDict()
    for j in 1:n, l in 1:n
        c = _CNUM_ZERO
        for k in 1:n
            c = _add_cnum(c, _mul_cnum(Wd[k, j], Wc[k, l]))
        end
        c = _mul_cnum(_CNUM_IM, c)
        _iszero_cnum(c) && continue
        _addto!(d, Op[mk(j, l)], c)
    end
    return QAdd(d, _EMPTY_INDICES)
end

_nlevel_gauge(σ::Op, W::AbstractMatrix, tt::Num) =
    _matrix_unit_gauge(W, Int(σ.nlev), (i, j) -> _transition_op(σ, i, j), tt)

function Rotation(σ::Op, W::AbstractMatrix, t::_ScalarLike)
    tt = _as_num(t)
    if is_collective_transition(σ)
        n = _collective_or_throw(σ, W, "`Rotation`")
        gauge = _matrix_unit_gauge(W, n, (i, j) -> _collective_op(σ, i, j), tt)
        return _with_gauge(_collective_rotation(σ, W), gauge, tt)
    end
    _nlevel_or_throw(σ, W, "`Rotation`")
    return _with_gauge(Rotation(σ, W), _nlevel_gauge(σ, W, tt), tt)
end

function Rotation(σ::Op, W::AbstractMatrix, t::_ScalarLike, gauge::QAdd)
    return _with_gauge(Rotation(σ, W), gauge, _as_num(t))
end

# === Canonicality self-test ===

# Ungated: only "is it zero" matters here, and the gate makes reduction non-monotone.
# See "`is_canonical` needs no representation of `U`" in the devdocs for why preserving the
# algebra and Hermiticity is sufficient, not merely necessary.
_residual(q::QAdd, U::UnitaryTransform) = _reduce_params(q, U.reductions, false)

function _fock_residuals!(out::Vector{QAdd}, d::Op, U::UnitaryTransform)
    lo = conjugate(d, U)
    hi = conjugate(adjoint(d), U)
    push!(out, _residual(commutator(lo, hi) - 1, U))
    push!(out, _residual(adjoint(lo) - hi, U))
    return nothing
end

# The `(Position, Momentum)` pair of a set of generators, or `nothing`. One helper for the
# three places that used to scan for it: coverage, residuals, and the rotating-frame basis.
function _quadrature_pair(gens)
    x = nothing
    p = nothing
    for o in gens
        o.kind === OP_POSITION && (x = o)
        o.kind === OP_MOMENTUM && (p = o)
    end
    return (x === nothing || p === nothing) ? nothing : (x, p)
end

function _phase_residuals!(out::Vector{QAdd}, site::SiteInfo, U::UnitaryTransform)
    pair = _quadrature_pair(site.generators)
    pair === nothing && return nothing
    x, p = pair
    xt = conjugate(x, U)
    pt = conjugate(p, U)
    push!(out, _residual(commutator(pt, xt) + im, U))
    # `x` and `p` are Hermitian, and the CCR alone does not pin that: `x -> x + im`
    # preserves it while no unitary induces it.
    push!(out, _residual(adjoint(xt) - xt, U))
    push!(out, _residual(adjoint(pt) - pt, U))
    return nothing
end

function _triple_residuals!(out::Vector{QAdd}, S::Op, U::UnitaryTransform)
    img = QAdd[conjugate(_axis_op(S, ax), U) for ax in 1:3]
    for j in 1:3
        push!(out, _residual(adjoint(img[j]) - img[j], U))
        for k in 1:3
            l = 6 - j - k
            if S.kind === OP_PAULI
                # σⱼσₖ = δⱼₖ + im*εⱼₖₗ*σₗ
                r = j == k ? img[j] * img[k] - 1 :
                    img[j] * img[k] - im * _levi_civita[j][k] * img[l]
                push!(out, _residual(r, U))
            elseif j < k
                # [Sⱼ, Sₖ] = im*εⱼₖₗ*Sₗ
                r = commutator(img[j], img[k]) - im * _levi_civita[j][k] * img[l]
                push!(out, _residual(r, U))
            end
        end
    end
    return nothing
end

function _nlevel_residuals!(out::Vector{QAdd}, σ::Op, U::UnitaryTransform)
    n = Int(σ.nlev)
    img = QAdd[conjugate(_transition_op(σ, i, j), U) for i in 1:n, j in 1:n]
    for i in 1:n, j in 1:n
        push!(out, _residual(adjoint(img[i, j]) - img[j, i], U))
    end
    # σ^{ij}σ^{kl} = δⱼₖ σ^{il}. Testing the products against the first row and column is
    # enough: σ^{ij} = σ^{i1}σ^{1j}, so associativity collapses the general product,
    # ẽ_{ij}ẽ_{kl} = ẽ_{i1}(ẽ_{1j}ẽ_{k1})ẽ_{1l} = δⱼₖ ẽ_{i1}ẽ_{1l} = δⱼₖ ẽ_{il}. The full
    # table costs n^4 products for no extra coverage.
    for i in 1:n, j in 1:n
        push!(out, _residual(img[i, 1] * img[1, j] - img[i, j], U))
        r = img[1, i] * img[j, 1]
        push!(out, _residual(i == j ? r - img[1, 1] : r, U))
    end
    completeness = sum(QAdd[img[i, i] for i in 1:n])
    push!(out, _residual(expand_completeness(completeness) - 1, U))
    return nothing
end

# A collective site satisfies the bracket law `[S^{ij}, S^{kl}] = δ_{jk}S^{il} - δ_{li}S^{kj}`
# and *not* the matrix-unit product law, so `_nlevel_residuals!` does not transfer. Three
# things pin it, and the third is not a consolation prize for the missing completeness
# identity: the contragredient map `S^{ij} -> -S^{ji}` satisfies the bracket exactly, is
# Hermitian-consistent exactly, and is involutive so it round-trips exactly, yet no unitary
# induces it. Only the invariance of the atom number `Σ_i S^{ii}` rejects it.
function _collective_residuals!(out::Vector{QAdd}, site::SiteInfo, U::UnitaryTransform)
    g = first(site.generators)
    levels = _collective_levels(site.generators)
    img = Dict{Tuple{Int, Int}, QAdd}(
        (i, j) => conjugate(_collective_op(g, i, j), U) for i in levels for j in levels
    )
    for i in levels, j in levels
        push!(out, _residual(adjoint(img[(i, j)]) - img[(j, i)], U))
    end
    # Brackets against a pivot row and column only. `F = {E_pj} ∪ {E_ip}` generates the whole
    # algebra (`E_ij = [E_ip, E_pj] + δ_ij E_pp`), so Jacobi carries preservation against `F`
    # to every pair: `2|L|-1` probes instead of `|L|^2`.
    p = first(levels)
    probes = unique!(Tuple{Int, Int}[(p, j) for j in levels])
    append!(probes, Tuple{Int, Int}[(i, p) for i in levels if i != p])
    for (f1, f2) in probes, k in levels, l in levels
        want = _zero_qadd()
        f2 == k && (want = want + img[(f1, l)])
        l == f1 && (want = want - img[(k, f2)])
        push!(out, _residual(commutator(img[(f1, f2)], img[(k, l)]) - want, U))
    end
    number = sum(QAdd[img[(i, i)] for i in levels])
    have = sum(QAdd[_single_qadd(_CNUM_ONE, Op[_collective_op(g, i, i)]) for i in levels])
    push!(out, _residual(number - have, U))
    return nothing
end

# Whether any rule sends a generator off its own site. A map that does not cannot break
# cross-site commutativity, so the quadratic sweep below is skipped for it entirely.
function _mixes_sites(U::UnitaryTransform)
    for (g, r) in U.rules
        k = _site_key(g)
        for (term, _) in r, o in term.ops
            _site_key(o) == k || return true
        end
    end
    return false
end

# Operators on different sites commute, and preserving the algebra of each site separately
# does not imply that. A mode-mixing map is exactly where it can fail: the whole content of a
# multimode frame is which combination each generator becomes, and a non-orthonormal one
# passes every same-site test.
function _cross_site_residuals!(out::Vector{QAdd}, U::UnitaryTransform, gs::Vector{Op})
    _mixes_sites(U) || return nothing
    img = QAdd[conjugate(g, U) for g in gs]
    for j in eachindex(gs), k in (j + 1):lastindex(gs)
        _site_key(gs[j]) == _site_key(gs[k]) && continue
        push!(out, _residual(commutator(img[j], img[k]), U))
    end
    return nothing
end

"""
    canonicality_residuals(U::UnitaryTransform) -> Vector{QAdd}

Expressions that vanish exactly when `U` maps its generators the way a unitary does: the
defining algebra of every site it touches, Hermiticity of the images, agreement between the
forward map and the stored inverse, and, for a map that mixes sites, commutativity between
the images of different ones.

Use [`is_canonical`](@ref) for the yes/no answer.
"""
function canonicality_residuals(U::UnitaryTransform)
    out = QAdd[]
    # A transform stores both directions, and preserving the algebra says nothing about the
    # second one: a wrong `inv_rules` passes every test above while `inv` returns garbage.
    Ui = inv(U)
    gs = generators(U)
    for g in gs
        push!(out, _residual(conjugate(conjugate(g, U), Ui) - g, U))
    end
    for site in U.sites
        g = first(site.generators)
        kind = g.kind
        if _is_fock(g)
            _fock_residuals!(out, _lowering(g), U)
        elseif _is_phase_space(g)
            _phase_residuals!(out, site, U)
        elseif kind === OP_PAULI || kind === OP_SPIN
            _triple_residuals!(out, g, U)
        elseif kind === OP_TRANSITION
            _nlevel_residuals!(out, g, U)
        elseif kind === OP_COLLECTIVE_TRANSITION
            _collective_residuals!(out, site, U)
        else
            throw(ArgumentError("no canonicality test for $(kind)"))
        end
    end
    _cross_site_residuals!(out, U, gs)
    return out
end

# An upper bound on `|x|` over all real arguments, or `nothing` when a factor is unbounded.
# `cos` and `sin` are the only transcendentals admitted: a frame phase is bounded by one, so
# a tiny scalar in front of it really does bound the term, while a tiny scalar in front of a
# `sinh` says nothing at all.
# A numeric literal reaches here wrapped, so `isa Number` alone would miss it.
_literal(x) = _const_value(SymbolicUtils.unwrap(x))

function _sup_norm(x)
    n = _literal(x)
    n === nothing || return Float64(abs(n))
    v = SymbolicUtils.unwrap(x)
    (v isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(v)) || return nothing
    op = SymbolicUtils.operation(v)
    args = SymbolicUtils.arguments(v)
    # A phase is unimodular. `_phase_poly_bound` covers it on the polynomial tier; this is
    # the same fact for a phase that reached the `Complex{Num}` tail.
    (length(args) == 1 && (op === cos || op === sin || op === expim)) && return 1.0
    # A root of a nonnegative literal is itself a constant. One survives whenever two exact
    # eigenvector entries with different radicands multiply.
    if op === ssqrt && length(args) == 1
        b = _literal(args[1])
        (b === nothing || !iszero(imag(b)) || real(b) < 0) && return nothing
        return sqrt(Float64(real(b)))
    end
    if op === (+) || op === (*)
        acc = op === (+) ? 0.0 : 1.0
        for a in args
            b = _sup_norm(a)
            b === nothing && return nothing
            acc = op === (+) ? acc + b : acc * b
        end
        return acc
    end
    if op === (^) && length(args) == 2
        e = _literal(args[2])
        e isa Integer || return nothing
        if e >= 0
            b = _sup_norm(args[1])
            return b === nothing ? nothing : b^e
        end
        # A negative power has no bound from a bound on the base, which can approach zero.
        # A nonzero literal base is the exception, and it is the common one: dividing by a
        # rate leaves `√κ^2` behind as a constant once the root has been undone.
        base = _literal(args[1])
        (base === nothing || iszero(base)) && return nothing
        return Float64(abs(base))^e
    end
    return nothing
end

# Bound of a polynomial whose every factor is a unit phase: the monomials cannot grow, so
# the sum of their scalars bounds the whole thing. `nothing` when some factor is not a
# phase, since an unconstrained parameter admits no bound.
function _phase_poly_bound(p::Poly)
    acc = 0.0
    for m in p.terms
        all(_is_phase, m.syms) || return nothing
        acc += abs(m.scalar)
    end
    return acc
end

# A residual built from floating-point coefficients cancels only to rounding, so a term
# bounded by `atol` counts as zero. Symbolic residuals are unaffected: an unbounded factor
# makes the bound `nothing`, and the parameters of an exact transform are never bounded.
function _bounded_by(c::CNum, atol::Real)
    _is_native(c) && return abs(c.z) <= atol
    # A phase polynomial is bounded term by term without ever lowering to SymbolicUtils.
    if c.tail isa Poly
        b = _phase_poly_bound(c.tail)
        b === nothing || return abs(c.z) + b <= atol
    end
    re, im = _realimag(c)
    br = _sup_norm(expand(re))
    br === nothing && return false
    bi = _sup_norm(expand(im))
    return bi !== nothing && br + bi <= atol
end

# Last resort, tried only once the coefficient tier has already failed to cancel. An
# involutive block divides by `√κ`, and `κ` is a sum, so its image carries `Δ^2/κ` and
# `4g^2/κ` as separate monomials that the polynomial tier never recombines into `1`. That is
# rational-function cancellation, which only SymbolicUtils does, and it is exact when it
# fires, so it has to be asked before any tolerance is.
#
# `ssqrt` is a square root, so `ssqrt(x)^2 = x` on either branch, but SymbolicUtils will not
# apply that on its own: `ssqrt` is opaque to it precisely so that a negative argument stays
# unevaluated. Supplying the one identity is what lets an involutive frame keep both halves,
# an exact certificate and a substitution that works below threshold.
const _SSQRT_SQUARE = SymbolicUtils.@rule ssqrt(~x)^2 => ~x

_unroot(x::Num) = _rewrite_with(x, _SSQRT_SQUARE)

# Same identities the polynomial tier applies, restated for a coefficient that reached the
# `Complex{Num}` tail. `_reduce_all` cannot find the pair once the terms sit inside divisions,
# and without it a rounded `0.3cos² + 0.3000000000000004sin²` stays a sum of two separately
# bounded terms that `_sup_norm` can only bound by their sum rather than by the rounding.
const _TRIG_SQUARE = SymbolicUtils.@rule cos(~x)^2 => 1 - sin(~x)^2
const _HYP_SQUARE = SymbolicUtils.@rule cosh(~x)^2 => 1 + sinh(~x)^2

_rewrite_with(x::Num, rules...) = Num(
    SymbolicUtils.Postwalk(
        SymbolicUtils.PassThrough(SymbolicUtils.Chain(collect(rules)))
    )(SymbolicUtils.unwrap(expand(x)))
)

# The bound of last resort: collapse the Pythagorean pairs first, so what is left is the
# rounding itself rather than two bounded terms that happen to cancel.
function _bounded_after_rewrite(c::CNum, atol::Real)
    re, im = _realimag(c)
    tot = 0.0
    for part in (re, im)
        y = _rewrite_with(part, _SSQRT_SQUARE, _TRIG_SQUARE, _HYP_SQUARE)
        b = _sup_norm(expand(simplify(y)))
        b === nothing && return false
        tot += b
    end
    return tot <= atol
end

function _simplifies_to_zero(c::CNum)
    _is_native(c) && return false
    re, im = _realimag(c)
    return _zero_after_simplify(re) && _zero_after_simplify(im)
end

# `simplify_fractions` before `simplify`: undoing the roots leaves the same term written two
# ways, once as a rational coefficient and once as a division, and only fraction handling puts
# them over a common denominator where they cancel.
function _zero_after_simplify(x::Num)
    y = _unroot(x)
    iszero(simplify(y)) && return true
    return iszero(simplify(expand(simplify_fractions(simplify(y)))))
end

# Per coefficient, not per tier: one residual can mix a rounded native term with a symbolic
# one that only cancels under `simplify`, and demanding that a single tier settle all of them
# would fail a residual whose every term is zero.
#
# The `cos`/`sin` pair is reduced first, ungated. A phase atom cancels on its own by adding
# exponents, but `cos` and `sin` are two unrelated factors, so `0.3cos² + 0.3sin²` stays a
# sum of two bounded terms and `_sup_norm` can only bound it by `0.6`. Applying the identity
# collapses it to the rounding it actually is, which `atol` then covers.
function _coeff_iszero(c::CNum, atol::Real, scratch::Vector{ParamRelation})
    r = _reduce_all(c, _NO_RELATIONS, false, scratch)
    _iszero_cnum(r) && return true
    _simplifies_to_zero(r) && return true
    return atol > 0 && (_bounded_by(r, atol) || _bounded_after_rewrite(r, atol))
end

const _NO_RELATIONS = ParamRelation[]

function _residual_iszero(q::QAdd, atol::Real)
    iszero(q) && return true
    scratch = ParamRelation[]
    return all(p -> _coeff_iszero(p.second, atol, scratch), q.arguments)
end

"""
    is_canonical(U::UnitaryTransform; atol = 1e-12) -> Bool

Whether `U` really is induced by a unitary: `true` iff every
[`canonicality_residuals`](@ref) entry vanishes.

The test runs on the transformed algebra, so it needs no matrix representation and covers
Fock and phase-space transforms too. A transform built from floating-point coefficients
cancels only to rounding, hence `atol`; pass `atol = 0` to demand exact cancellation.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables r;

julia> is_canonical(Squeeze(a, r))
true
```
"""
# Not worth short-circuiting on the first nonzero residual: building them dwarfs testing
# them, and a `true` answer has to build and test every one anyway.
is_canonical(U::UnitaryTransform; atol::Real = 1.0e-12) =
    all(q -> _residual_iszero(q, atol), canonicality_residuals(U))

# === Rotating frames ===
#
# The map is derived from the adjoint action, not by classifying term shapes. For
# `U = exp(-im*H0*t)` every generator obeys `dg̃/dt = im*[H0, g̃]`, so when `ad_{H0}` closes
# on the generators of the sites `H0` touches it is an `m x m` matrix `A` there and the
# whole rule map is `exp(im*A*t)`. A commutator term that is not one of those generators is
# exactly what puts a frame outside the exactly solvable class, and it is named as such.

# A site that cannot describe itself from one member takes its generators from `H0`, which is
# where the missing information (the other quadrature name) actually appears.
function _frame_site_generators(H0::QAdd, o::Op)
    gens = _site_generators(o)
    isempty(gens) || return gens
    k = _site_key(o)
    present = Op[]
    for (term, _) in H0, q in term.ops
        (_site_key(q) == k && !any(g -> isequal(g, q), present)) && push!(present, q)
    end
    if is_collective_transition(o)
        # `[S^{ij}, S^{kl}]` stays inside `L×L` for `i,j,k,l ∈ L`, so any full square is a
        # subalgebra and `ad` closes on it. The square has to be guessed, though: the level
        # count lives on the space and a collective `Op` does not reference it. `1:maximum`
        # is the guess, matching what the level-diagonal path infers, so an `H0` naming only
        # `S^{22}` still covers `S^{12}`. A level above every one `H0` names is invisible and
        # comes back uncovered.
        n = maximum(_collective_levels(present))
        return Op[_collective_op(o, i, j) for i in 1:n for j in 1:n]
    end
    _is_phase_space(o) || throw(
        ArgumentError("no frame basis for a $(o.kind) site; `[H0, ·]` cannot be closed on it")
    )
    # A `PhaseSpace` stores its own name and not its quadratures', so the conjugate's name
    # exists nowhere but in `H0`. A zero-coefficient term will not carry it either: that is
    # dropped at construction.
    _quadrature_pair(present) === nothing && throw(
        ArgumentError(
            "`H0` uses only one member of the phase-space pair on the site of `$o`; a " *
                "frame needs both the position and the momentum, and the conjugate's name " *
                "is recoverable only from `H0` itself"
        )
    )
    x, p = _quadrature_pair(present)
    return Op[x, p]
end

# The generators of every site `H0` touches, in a fixed order.
function _frame_basis(H0::QAdd)
    basis = Op[]
    seen = _SiteKey[]
    for (term, _) in H0, o in term.ops
        k = _site_key(o)
        (k in seen) && continue
        push!(seen, k)
        append!(basis, _frame_site_generators(H0, o))
    end
    isempty(basis) && throw(
        ArgumentError("`H0` has no operator term; there is no frame to rotate into")
    )
    return basis
end

_basis_index(pos::Dict{Op, Int}, t::QTerm) =
    (length(t.ops) == 1 && isempty(t.ne)) ? get(pos, t.ops[1], nothing) : nothing

# `[G, g_k] = Σ_j A[j,k] g_j + c_k`. The commutator is built through `*`, so it is already
# canonical and `normal_order` on it would be a second full pass for the same form. The
# inhomogeneous `c` is what a drive produces, and only the Heisenberg algebra produces one:
# no commutator of transitions, Pauli or spin operators is a bare c-number.
function _ad_affine(G::QAdd, basis::Vector{Op})
    m = length(basis)
    A = fill(_CNUM_ZERO, m, m)
    c = fill(_CNUM_ZERO, m)
    pos = Dict{Op, Int}(g => j for (j, g) in enumerate(basis))
    for (k, g) in enumerate(basis)
        for (term, coeff) in commutator(G, g)
            if isempty(term.ops) && isempty(term.ne)
                c[k] = _add_cnum(c[k], coeff)
                continue
            end
            j = _basis_index(pos, term)
            j === nothing && throw(
                ArgumentError(
                    "`H0` does not generate an exactly solvable frame: `[H0, $g]` " *
                        "produced `$(_single_qadd(_CNUM_ONE, term.ops, term.ne))`, which " *
                        "is not one of the generators `H0` acts on"
                )
            )
            A[j, k] = _add_cnum(A[j, k], coeff)
        end
    end
    return (A, c)
end

function _frame_rate(λ::CNum, g::Op)
    _iszero_num(imag(λ)) || throw(
        ArgumentError(
            "the frame rate of `$g` is `$(to_num(λ))`, which is not real; " *
                "`exp(-im*H0*t)` is unitary only for a Hermitian `H0`"
        )
    )
    return real(λ)
end

# Connected components of `A`'s sparsity graph: the invariant subspaces the exponential
# factors over.
function _ad_components(A::Matrix{CNum})
    m = size(A, 1)
    parent = collect(1:m)
    root(x) = (
        while parent[x] != x
            x = parent[x]
        end; x
    )
    for j in 1:m, k in 1:m
        (j == k || _iszero_cnum(A[j, k])) && continue
        a, b = root(j), root(k)
        a == b || (parent[a] = b)
    end
    groups = Dict{Int, Vector{Int}}()
    for j in 1:m
        push!(get!(() -> Int[], groups, root(j)), j)
    end
    return sort!(collect(values(groups)); by = first)
end

# `maximum(abs, C - C')` would allocate two more `m x m` matrices to answer this.
function _is_hermitian(C::Matrix{ComplexF64})
    scale = 0.0
    for v in C
        scale = max(scale, abs(v))
    end
    tol = 1.0e-10 * max(1.0, scale)
    for u in axes(C, 1), v in axes(C, 2)
        abs(C[u, v] - conj(C[v, u])) <= tol || return false
    end
    return true
end

# A rate times a numeric unit, or `nothing` when the entry is neither real nor imaginary.
function _split_scalar(c::CNum)
    re, im = real(c), imag(c)
    _iszero_num(im) && return (re, ComplexF64(1, 0))
    _iszero_num(re) && return (im, ComplexF64(0, 1))
    return nothing
end

# A zero-diagonal 2x2 block is `w * B` with numeric `B`, so the rate is read off instead of
# square-rooted out of the determinant. `B^2 = I` exponentiates to cos/sin and `B^2 = -I`
# to cosh/sinh, which is the whole difference between a rotation and a squeeze.
function _block_factor(A::Matrix{CNum}, j::Int, k::Int)
    (_iszero_cnum(A[j, j]) && _iszero_cnum(A[k, k])) || return nothing
    sβ = _split_scalar(A[j, k])
    sγ = _split_scalar(A[k, j])
    (sβ === nothing || sγ === nothing) && return nothing
    w, bjk = sβ
    q = _to_cnum(sγ[1] / w)
    _is_native(q) || return nothing
    bkj = sγ[2] * _to_complex(q)
    s = bjk * bkj
    (s ≈ 1 || s ≈ -1) || return nothing
    # `w*B` is invariant under flipping both, so prefer the form without a leading minus.
    if _is_real_negative(_to_cnum(w))
        w, bjk, bkj = -w, -bjk, -bkj
    end
    return (w, bjk, bkj, real(s) > 0)
end

# One phase per level, not one per transition. `σ^{ij}` rotates at `E_i - E_j`, so a frame
# built per generator mints an atom for every difference and `E₁-E₃ = (E₁-E₂) + (E₂-E₃)`
# stops being visible: the residuals would then need the angle-addition identity. Carrying
# `P_i * conj(P_j)` instead makes every such relation exponent arithmetic on `n` atoms.
function _level_phase_rules(M::Matrix{Complex{Num}}, n::Int, mk::F, θ::Num) where {F}
    # Every other frame path rejects a non-Hermitian `G` (`_frame_rate`, `_eigen_or_throw`,
    # `_dressed_symbolic`). Without this one, `real(M[k, k])` silently drops the imaginary
    # part and the result reports itself canonical.
    all(k -> _iszero_num(imag(M[k, k])), 1:n) || throw(
        ArgumentError(
            "an exactly solvable frame needs a Hermitian generator; its level energies " *
                "are not real"
        )
    )
    P = [_phase_coeff(real(M[k, k]) * θ) for k in 1:n]
    rules = Dict{Op, QAdd}()
    inv_rules = Dict{Op, QAdd}()
    for i in 1:n, j in 1:n
        g = mk(i, j)
        c = _mul_cnum(P[i], _conj_cnum(P[j]))
        rules[g] = _scaled(c, g)
        inv_rules[g] = _scaled(_conj_cnum(c), g)
    end
    return (rules, inv_rules)
end

_level_phase_rules(σ::Op, M::Matrix{Complex{Num}}, θ::Num) =
    _level_phase_rules(M, Int(σ.nlev), (i, j) -> _transition_op(σ, i, j), θ)

function _frame_diagonal_rule!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd},
        g::Op, λ::CNum, θ::Num,
    )
    if _iszero_cnum(λ)
        merge!(rules, _with_adjoint(g, _scaled(_CNUM_ONE, g)))
        merge!(inv_rules, _with_adjoint(g, _scaled(_CNUM_ONE, g)))
        return nothing
    end
    if isequal(adjoint(g), g)
        # A Hermitian generator that stays Hermitian can carry no phase at all. An imaginary
        # rate is instead a real scaling, `dg̃/dθ = im*λ*g̃` with `im*λ` real, which is the
        # squeeze of a quadrature. Only a quadrature reaches here: a projector is never an
        # `ad` eigenvector (`[E_kl, E_ii]` cannot return `E_ii`) and a spin axis always
        # couples to the other two, so both land on a block instead of a singleton.
        _iszero_num(real(λ)) || throw(
            ArgumentError(
                "`$g` is Hermitian but `[G, $g]` is a real multiple of it, which no " *
                    "unitary can produce"
            )
        )
        y = -imag(λ) * θ
        merge!(rules, _with_adjoint(g, _scaled(_real_scale(y), g)))
        merge!(inv_rules, _with_adjoint(g, _scaled(_real_scale(-y), g)))
        return nothing
    end
    # One phase atom per rate. The inverse is its conjugate, which is the same atom at the
    # opposite exponent, so the two cancel on multiplication with no relation to declare.
    # `_frame_rate` rejects a non-real rate here: a real scaling of a ladder operator scales
    # its commutator by the square and is never unitary.
    p = _phase_coeff(_frame_rate(λ, g) * θ)
    merge!(rules, _with_adjoint(g, _scaled(p, g)))
    merge!(inv_rules, _with_adjoint(g, _scaled(_conj_cnum(p), g)))
    return nothing
end

function _frame_block_rule!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, rels::Vector{ParamRelation},
        basis::Vector{Op}, A::Matrix{CNum}, j::Int, k::Int, tt::Num,
    )
    f = _block_factor(A, j, k)
    f === nothing && return false
    w, bjk, bkj, rotates = f
    θ = w * tt
    # `exp(im*w*t*B) = dg*I + im*og*B`; the inverse is even in `dg` and odd in `og`, so it
    # reuses the same atoms rather than negating the argument. A rotating block has
    # eigenphases `exp(±im*θ)`, so it is written on one phase atom (`cos θ = (p + p⁻¹)/2`,
    # `im*sin θ = (p - p⁻¹)/2`) and needs no relation. A hyperbolic block is not a phase.
    rotates || push!(rels, _hyp_rel(θ))
    dg, img = if rotates
        p = _phase_coeff(θ)
        pc = _conj_cnum(p)
        (
            _mul_cnum(_CNUM_HALF, _add_cnum(p, pc)),
            _mul_cnum(_CNUM_HALF, _add_cnum(p, _neg_cnum(pc))),
        )
    else
        (_to_cnum(cosh(θ)), _mul_cnum(_CNUM_IM, _to_cnum(sinh(θ))))
    end
    for (u, v, b) in ((j, k, bkj), (k, j, bjk))
        c = _mul_cnum(img, _to_cnum(b))
        rules[basis[u]] = _rule_qadd((dg, Op[basis[u]]), (c, Op[basis[v]]))
        inv_rules[basis[u]] = _rule_qadd((dg, Op[basis[u]]), (_neg_cnum(c), Op[basis[v]]))
    end
    return true
end

# `A[comp, comp]` as `d*I + w*C` with a real rate `w`, a numeric `C`, and a shared diagonal
# `d`, or `nothing` when no such split exists. See "Any number of generators may couple" in
# the devdocs for why the diagonal comes out first. Hermiticity of `C` is what the `eigen`
# path needs and the involutive one does not, so it is the caller's question.
function _common_factor(A::Matrix{CNum}, comp::Vector{Int})
    m = length(comp)
    dg = A[comp[1], comp[1]]
    d = all(j -> _iszero_cnum(_add_cnum(A[j, j], _neg_cnum(dg))), comp) ? dg : _CNUM_ZERO
    B = Matrix{CNum}(undef, m, m)
    for u in 1:m, v in 1:m
        e = A[comp[u], comp[v]]
        B[u, v] = u == v ? _add_cnum(e, _neg_cnum(d)) : e
    end
    ref = nothing
    for u in 1:m, v in 1:m
        (u == v || _iszero_cnum(B[u, v])) && continue
        ref = _split_scalar(B[u, v])
        break
    end
    ref === nothing && return nothing
    w = ref[1]
    C = Matrix{ComplexF64}(undef, m, m)
    for u in 1:m, v in 1:m
        s = _split_scalar(B[u, v])
        s === nothing && return nothing
        q = _to_cnum(s[1] / w)
        _is_native(q) || return nothing
        C[u, v] = s[2] * _to_complex(q)
    end
    return (w, C, d)
end

function _common_factor_hermitian(A::Matrix{CNum}, comp::Vector{Int})
    f = _common_factor(A, comp)
    (f === nothing || !_is_hermitian(f[2])) && return nothing
    return f
end

# `A[comp, comp]` split as `d*I + C` with the shared diagonal `d` pulled out, or a zero `d`
# when the diagonal is not shared. Unlike `_common_factor` this never fails and never reaches
# for a `ComplexF64`: the involution test below is a polynomial identity, so it needs neither
# a rate divided out nor numeric entries.
function _shared_diagonal(A::Matrix{CNum}, comp::Vector{Int})
    m = length(comp)
    dg = A[comp[1], comp[1]]
    d = all(j -> _iszero_cnum(_add_cnum(A[j, j], _neg_cnum(dg))), comp) ? dg : _CNUM_ZERO
    C = Matrix{CNum}(undef, m, m)
    for u in 1:m, v in 1:m
        e = A[comp[u], comp[v]]
        C[u, v] = u == v ? _add_cnum(e, _neg_cnum(d)) : e
    end
    return (C, d)
end

# `κ` with `C^2 = κ*I`, or `nothing`. An involution exponentiates in closed form at any size
# and with any `κ`, where `eigen` needs a Hermitian `C` and returns numeric eigenvectors.
# `κ` is whatever the coefficients make it, symbolic included, and the test is exact
# cancellation rather than a tolerance: a float reaches here only if the user wrote one.
function _involution_scale(C::Matrix{CNum})
    m = size(C, 1)
    κ = _CNUM_ZERO
    for u in 1:m, v in 1:m
        s = _CNUM_ZERO
        for k in 1:m
            s = _add_cnum(s, _mul_cnum(C[u, k], C[k, v]))
        end
        if u == v
            u == 1 ? (κ = s) : (_iszero_cnum(_add_cnum(s, _neg_cnum(κ))) || return nothing)
        else
            _iszero_cnum(s) || return nothing
        end
    end
    return κ
end

# `√κ` as an exact, real, rounding-free quantity, or `nothing`. A single monomial always has
# one: the polynomial tier carries `Rational` exponents, so halving them is exact and only the
# scalar has to come out whole. This is what keeps `√(g^2)` equal to `g` rather than an opaque
# radical, and with it a rotation stays on one phase atom instead of splitting into `cos`/`sin`.
function _exact_root(κ::CNum)
    if _is_native(κ)
        iszero(imag(κ.z)) && real(κ.z) >= 0 || return nothing
        r = sqrt(real(κ.z))
        return isinteger(r) ? Num(Int(r)) : nothing
    end
    # An exact rational `κ` is deliberately *not* rooted here even when it could be: a
    # rational rate makes `_phase_coeff` fold only part of the pair back to `cos`/`sin`, so
    # the two images end up spelled differently and the residual stops closing. It keeps the
    # radical and certifies exactly, which is the property that matters.
    _is_poly(κ) || return nothing
    ts = κ.tail.terms
    length(ts) == 1 || return nothing
    m = ts[1]
    (iszero(imag(m.scalar)) && real(m.scalar) > 0) || return nothing
    rs = sqrt(real(m.scalar))
    isinteger(rs) || return nothing
    return real(to_num(_from_poly(Monomial[Monomial(complex(rs), m.syms, m.exps .// 2)])))
end

# `√κ` as an unevaluated node. Built with `term` rather than by calling `ssqrt`, for two
# reasons: `ssqrt` returns `nothing` on a `Num`-wrapped constant, and on a bare positive one it
# would fold to a `Float64` and put the rounding back. Left unevaluated it stays exact through
# the residuals and still evaluates, on either branch, once parameters are substituted.
_safe_root(x::Num) = Num(SymbolicUtils.term(ssqrt, SymbolicUtils.unwrap(x)))

# The pair `(dg, og)` in `exp(im*θ*C) = dg*I + og*C` for `C^2 = κ*I`, from
# `cos(√κ θ)` and `im*sin(√κ θ)/√κ`. Both are entire in `κ`, so one formula covers every
# sign; the branches below only pick the spelling that reads best and closes cheapest.
function _involution_pair(κ::CNum, θ::Num, rels::Vector{ParamRelation})
    # Nilpotent: the series stops after one term, so this is a polynomial with no atom at all
    # and it cancels exactly. A parametric amplifier driven exactly at threshold (`Δ = 2g`)
    # lands here, which is why the closed form is worth having.
    _iszero_cnum(κ) && return (_CNUM_ONE, _mul_cnum(_CNUM_IM, _to_cnum(θ)))
    μ = _exact_root(κ)
    if μ !== nothing
        # A real root means a genuine phase: one atom whose conjugate is itself at the
        # opposite exponent, so the two cancel with no relation to declare. Worth taking
        # whenever it is available, since `cos`/`sin` would need one and read worse.
        p = _phase_coeff(μ * θ)
        pc = _conj_cnum(p)
        return (
            _mul_cnum(_CNUM_HALF, _add_cnum(p, pc)),
            _cnum_div(_mul_cnum(_CNUM_HALF, _add_cnum(p, _neg_cnum(pc))), _to_cnum(μ)),
        )
    end
    ν = _exact_root(_neg_cnum(κ))
    if ν !== nothing
        # `√κ` is imaginary, so the phase spelling would be a lie: `expim` declares a unit
        # modulus and `conj` would negate the wrong thing. Real `cosh`/`sinh` instead. Only
        # on an exact root: `ssqrt(35)` folds to a `Float64` where `ssqrt(-35)` does not, so
        # an irrational one is better served by the unevaluated `√κ` in the branch below.
        y = ν * θ
        push!(rels, _hyp_rel(y))
        return (
            _to_cnum(cosh(y)),
            _cnum_div(_mul_cnum(_CNUM_IM, _to_cnum(sinh(y))), _to_cnum(ν)),
        )
    end
    # Sign unknown, or known but with an irrational root. `cos(√κ θ)` and `sin(√κ θ)/√κ` are
    # even in `√κ`, hence single-valued in `κ`, so the branch is irrelevant to the value and
    # `cos² + sin² = 1` closes on either. `ssqrt` rather than `sqrt` so that substituting
    # parameters that make `κ` negative evaluates instead of raising a `DomainError`.
    kn = to_num(κ)
    # A Hermitian generator gives a real `κ`; anything else is not a frame this path owns.
    _iszero_num(imag(kn)) || return nothing
    μ = _safe_root(real(kn))
    x = μ * θ
    push!(rels, _trig_rel(x))
    return (
        _to_cnum(cos(x)),
        _cnum_div(_mul_cnum(_CNUM_IM, _to_cnum(sin(x))), _to_cnum(μ)),
    )
end

# `exp(im*θ*C)` for an involutive `C`, on any number of coupled generators. `C^2 = κ*I` makes
# the series close on `I` and `C`, so the answer is exact and stays symbolic where
# `_frame_multi_rule!` would need a Hermitian `C` and numeric eigenvectors. No rate is divided
# out: `κ` absorbs any scale, which is what lets `C` keep the coefficients it was built from.
function _frame_involutive_rule!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, rels::Vector{ParamRelation},
        basis::Vector{Op}, A::Matrix{CNum}, comp::Vector{Int}, θ::Num,
    )
    C, d = _shared_diagonal(A, comp)
    κ = _involution_scale(C)
    κ === nothing && return false
    pair = _involution_pair(κ, θ, rels)
    pair === nothing && return false
    m = length(comp)
    dg, og = pair
    pd = _iszero_cnum(d) ? _CNUM_ONE : _phase_coeff(_frame_rate(d, basis[comp[1]]) * θ)
    pdc = _conj_cnum(pd)
    # `ad` acts on coefficient vectors, so the image of a generator is a column.
    for u in 1:m
        fwd = Tuple{CNum, Vector{Op}}[]
        bwd = Tuple{CNum, Vector{Op}}[]
        for v in 1:m
            off = _mul_cnum(og, _to_cnum(C[v, u]))
            # `dg` is even in the rate and `og` odd, so the inverse negates only the latter
            # and reuses both atoms rather than exponentiating a second time.
            e = u == v ? _add_cnum(dg, off) : off
            ei = u == v ? _add_cnum(dg, _neg_cnum(off)) : _neg_cnum(off)
            g = Op[basis[comp[v]]]
            push!(fwd, (_mul_cnum(pd, e), g))
            push!(bwd, (_mul_cnum(pdc, ei), g))
        end
        rules[basis[comp[u]]] = _rule_qadd(fwd)
        inv_rules[basis[comp[u]]] = _rule_qadd(bwd)
    end
    return true
end

# === Structured eigendecompositions ===

# An exact real as `q*√r` with `r` squarefree, kept in that normal form throughout: `Symbolics`
# combines neither `√(1//2)*√2` nor `√2*√2`, so leaving the arithmetic to it would build
# eigenvector entries that are correct and unreadable.
const _RootVal = Tuple{Rational{Int}, Int}

# Every square factor of `r` pulled into the rational part.
function _root_reduce(q::Rational{Int}, r::Int)
    iszero(q) && return (zero(Rational{Int}), 1)
    for p in 2:isqrt(r)
        while r % (p * p) == 0
            r ÷= p * p
            q *= p
        end
    end
    return (q, r)
end

_root_mul(a::_RootVal, b::_RootVal) = _root_reduce(a[1] * b[1], a[2] * b[2])

# A root of a positive literal, typed `Real` so `_conj_atom` folds it. Without the symtype the
# conjugation cannot tell a radical from an unknown complex atom and leaves a `conj(...)`
# wrapper that never cancels, which breaks every residual built from these entries.
_real_root(x::Num) = Num(SymbolicUtils.term(ssqrt, SymbolicUtils.unwrap(x); type = Real))

_root_to_num(v::_RootVal) =
    iszero(v[1]) ? _NUM_ZERO : (v[2] == 1 ? Num(v[1]) : v[1] * _real_root(Num(v[2])))

# `sinpi(r)` as a `_RootVal`, or `nothing`. Only the rational multiples whose sine is a single
# rational multiple of a single root are listed. `sinpi(1//12)` and `sinpi(1//10)` are sums of
# two roots and `sinpi(1//8)` is a nested one, and neither fits the normal form; a matrix
# needing them falls back to `eigen` exactly as before. Denominators 1, 2, 3, 4 and 6 covered.
function _sinpi_exact(r::Rational{Int})
    s = mod(r, 2)
    neg = s > 1
    neg && (s -= 1)
    s > 1 // 2 && (s = 1 - s)
    v = if s == 0
        (zero(Rational{Int}), 1)
    elseif s == 1 // 2
        (1 // 1, 1)
    elseif s == 1 // 6
        (1 // 2, 1)
    elseif s == 1 // 4
        (1 // 2, 2)
    elseif s == 1 // 3
        (1 // 2, 3)
    else
        return nothing
    end
    return neg ? (-v[1], v[2]) : v
end

_cospi_exact(r::Rational{Int}) = _sinpi_exact(1 // 2 - r)

# A real matrix entry as a `Num`, keeping an integral value an `Int`. Same atom-identity reason
# as `_scaled_rate`: a `2.0` scale would mint a second phase atom beside the one `2` builds.
_exact_real(x::ComplexF64) = (r = real(x); isinteger(r) ? Num(Int(r)) : Num(r))

# `(a, b)` with `C == a*I + b*(shift + shift')`, or `nothing`. Entries compare exactly: a
# uniform coupling divided out by `_common_factor` lands on a literal, and a tolerance here
# would answer exactly for a matrix that is not the structure.
function _tridiag_toeplitz(C::Matrix{ComplexF64})
    n = size(C, 1)
    n >= 2 || return nothing
    a, b = C[1, 1], C[1, 2]
    iszero(b) && return nothing
    for u in 1:n, v in 1:n
        d = abs(u - v)
        want = d == 0 ? a : (d == 1 ? b : zero(ComplexF64))
        C[u, v] == want || return nothing
    end
    return (a, b)
end

# The first row of `C` when every row is a cyclic shift of it, or `nothing`. Any circulant is
# diagonalized by the DFT, not only the nearest-neighbour band, so a ring with longer-range
# hopping qualifies too. Below size three a shift lands back on the band and the class is
# degenerate.
function _circulant_row(C::Matrix{ComplexF64})
    n = size(C, 1)
    n >= 3 || return nothing
    for u in 2:n, v in 1:n
        C[u, v] == C[1, mod1(v - u + 1, n)] || return nothing
    end
    return ComplexF64[C[1, v] for v in 1:n]
end

# The discrete sine transform: `V[j, k] ∝ sin(jkπ/(n+1))`, `E[k] = a + 2b*cos(kπ/(n+1))`. These
# are the Chebyshev polynomials of the second kind evaluated on the spectrum.
function _sine_eigen(n::Int, a::ComplexF64, b::ComplexF64)
    (iszero(imag(a)) && iszero(imag(b))) || return nothing
    an, bn = _exact_real(a), _exact_real(b)
    # `√(2/(n+1))` in normal form.
    nrm = _root_reduce(1 // (n + 1), 2 * (n + 1))
    V = Matrix{Complex{Num}}(undef, n, n)
    E = Vector{Num}(undef, n)
    for k in 1:n
        c = _cospi_exact(k // (n + 1))
        c === nothing && return nothing
        E[k] = an + bn * _root_to_num(_root_mul((2 // 1, 1), c))
        for j in 1:n
            s = _sinpi_exact((j * k) // (n + 1))
            s === nothing && return nothing
            V[j, k] = Complex(_root_to_num(_root_mul(nrm, s)), _NUM_ZERO)
        end
    end
    return (V, E)
end

# The DFT: `V[j, k] ∝ exp(2πi*jk/n)` whatever the row, and `E[k] = Σⱼ cⱼ*exp(2πi*jk/n)`. The
# row must be real, and its sine sum must vanish, or the eigenvalues are not real and the
# callers (all of which need a Hermitian block) would be handed a different problem.
function _dft_eigen(n::Int, row::Vector{ComplexF64})
    all(c -> iszero(imag(c)), row) || return nothing
    cs = Num[_exact_real(c) for c in row]
    nrm = _root_reduce(1 // n, n)
    V = Matrix{Complex{Num}}(undef, n, n)
    E = Vector{Num}(undef, n)
    for k in 1:n
        # Wavenumbers run from zero, unlike the sine transform's modes.
        p = 2 * (k - 1)
        re = _NUM_ZERO
        im = _NUM_ZERO
        for j in 0:(n - 1)
            rc = _cospi_exact((j * p) // n)
            rs = _sinpi_exact((j * p) // n)
            (rc === nothing || rs === nothing) && return nothing
            re += cs[j + 1] * _root_to_num(rc)
            im += cs[j + 1] * _root_to_num(rs)
        end
        _iszero_num(_unroot(expand(im))) || return nothing
        E[k] = re
        for j in 1:n
            rc = _cospi_exact((j * p) // n)
            rs = _sinpi_exact((j * p) // n)
            (rc === nothing || rs === nothing) && return nothing
            V[j, k] = Complex(
                _root_to_num(_root_mul(nrm, rc)), _root_to_num(_root_mul(nrm, rs))
            )
        end
    end
    return (V, E)
end

# Exact eigenvectors and eigenvalues of `C`, or `nothing` when it is neither structure. The
# eigenvector matrix of a uniform band depends on the size alone, so it is exact whatever `a`
# and `b` are; only the eigenvalues read them, and they stay as exact as what the user wrote.
# A uniform nearest-neighbour chain reduces to the open band and the same chain with periodic
# boundaries to the closed one, which is what makes those two cases worth detecting.
function _structured_eigen(C::Matrix{ComplexF64})
    f = _tridiag_toeplitz(C)
    f === nothing || return _sine_eigen(size(C, 1), f[1], f[2])
    r = _circulant_row(C)
    r === nothing || return _dft_eigen(size(C, 1), r)
    return nothing
end

# Whether every weight `V[u,p]*conj(V[v,p])` folds to a plain number. True for a size whose
# eigenvector entries share one radicand (the squares then cancel), false as soon as two
# different roots meet.
# Whether every entry itself is a plain number. Stricter than `_rational_weights`, and the
# right test when the entries will be multiplied pairwise in every combination.
function _all_rational(V::Matrix{Complex{Num}})
    for e in V
        _is_native(_to_cnum(e)) || return false
    end
    return true
end

function _rational_weights(V::Matrix{Complex{Num}})
    m = size(V, 1)
    for u in 1:m, v in 1:m, p in 1:m
        c = _eigen_weight(V, u, v, p)
        (c === nothing || _is_native(c)) || return false
    end
    return true
end

# One term of a spectral outer product, or `nothing` when it drops out. Numeric eigenvectors
# carry rounding noise, so a tiny weight there is a structural zero the floats missed; the
# exact ones hit zero on the nose and need no threshold.
function _eigen_weight(V::Matrix{ComplexF64}, u::Int, v::Int, p::Int)
    ν = V[u, p] * conj(V[v, p])
    return abs(ν) <= 1.0e-12 ? nothing : _to_cnum(ν)
end
function _eigen_weight(V::Matrix{Complex{Num}}, u::Int, v::Int, p::Int)
    ν = _fold_roots(_mul_cnum(_to_cnum(V[u, p]), _conj_cnum(_to_cnum(V[v, p]))))
    return _iszero_cnum(ν) ? nothing : ν
end

# `√x*√x -> x` inside a coefficient. Two eigenvector entries multiplied square the
# normalization back to a rational, and that is the one identity `SymbolicUtils` will not
# apply on its own. A native coefficient carries no root, so the common path pays one branch.
function _fold_roots(c::CNum)
    _is_native(c) && return c
    re, im = _realimag(c)
    return _to_cnum(Complex(_unroot(re), _unroot(im)))
end

# Snap the spectrum onto shared magnitudes: `eigen` returns `±√2` as two floats that are not
# exact negatives, which would mint two unrelated atoms instead of one at opposite exponents.
function _snap_spectrum!(μ::Vector{Float64})
    tol = 1.0e-10 * max(1.0, maximum(abs, μ))
    reps = Float64[]
    for p in eachindex(μ)
        a = abs(μ[p])
        if a <= tol
            μ[p] = 0.0
            continue
        end
        # An integral magnitude has to land exactly on its integer, so that `_scaled_rate`
        # below can keep the scale an `Int`. A circulant ring returns `1.0000000000000002`.
        abs(a - round(a)) <= tol && (a = round(a))
        i = 0
        for n in eachindex(reps)
            if abs(reps[n] - a) <= tol
                i = n
                break
            end
        end
        i == 0 ? push!(reps, a) : (a = reps[i])
        μ[p] = sign(μ[p]) * a
    end
    return μ
end

# Atom identity is `===`, and `(-g)*(-1.0)` is equal to `g` and hashes alike without being the
# same object, so a `Float64` scale of magnitude one mints a second `expim(g*t)` that never
# merges with the first. An integral eigenvalue therefore scales the rate as an `Int`.
_scaled_rate(w::Num, μ::Float64) = isinteger(μ) ? Int(μ) * w : μ * w

# `exp(im*(d + w*C)*θ)` on a component, from the spectral decomposition of `C`: each eigenvalue
# gets one phase atom and the image of a generator is a column of `V*diag(P)*V'`.
function _emit_eigen_rules!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, basis::Vector{Op},
        comp::Vector{Int}, V::AbstractMatrix, rates::Vector{Num}, pd::CNum, tt::Num,
    )
    m = length(comp)
    P = CNum[
        _iszero_num(rates[p]) ? pd : _mul_cnum(pd, _phase_coeff(rates[p] * tt)) for p in 1:m
    ]
    M = Matrix{CNum}(undef, m, m)
    for u in 1:m, v in 1:m
        c = _CNUM_ZERO
        for p in 1:m
            ν = _eigen_weight(V, u, v, p)
            ν === nothing && continue
            c = _add_cnum(c, _mul_cnum(ν, P[p]))
        end
        M[u, v] = c
    end
    # `ad` acts on coefficient vectors, so the image of a generator is a column of `M`, and
    # the inverse map is `M'` because `M` is unitary.
    for u in 1:m
        fwd = Tuple{CNum, Vector{Op}}[]
        bwd = Tuple{CNum, Vector{Op}}[]
        for v in 1:m
            push!(fwd, (M[v, u], Op[basis[comp[v]]]))
            push!(bwd, (_conj_cnum(M[u, v]), Op[basis[comp[v]]]))
        end
        rules[basis[comp[u]]] = _rule_qadd(fwd)
        inv_rules[basis[comp[u]]] = _rule_qadd(bwd)
    end
    return nothing
end

# `exp(im*A*t)` on a component of any size: diagonalize `A = d*I + w*C`. A structured `C`
# answers exactly; otherwise the eigenvectors are numbers and these residuals cancel to
# rounding rather than exactly, which is what `is_canonical`'s `atol` is for.
function _frame_multi_rule!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd},
        basis::Vector{Op}, A::Matrix{CNum}, comp::Vector{Int}, tt::Num,
    )
    f = _common_factor_hermitian(A, comp)
    f === nothing && throw(
        ArgumentError(
            "`H0` couples $(join(string.(basis[comp]), ", ")) in a way this frame cannot " *
                "exponentiate in closed form: the block is neither an involution nor a " *
                "real rate times a numeric Hermitian matrix plus a shared diagonal. Give " *
                "numeric couplings, or build the transform with an explicit constructor"
        )
    )
    w, C, d = f
    pd = _iszero_cnum(d) ? _CNUM_ONE : _phase_coeff(_frame_rate(d, basis[comp[1]]) * tt)
    s = _structured_eigen(C)
    # A radical weight carries the coefficient off the polynomial tier, and that tier is where
    # a phase cancels against its conjugate. Exactness and phase cancellation are not both
    # available, so a frame takes the exact spectrum only when its weights are rational.
    s === nothing || _rational_weights(s[1]) || (s = nothing)
    if s === nothing
        F = eigen(Hermitian(C))
        μ = _snap_spectrum!(F.values)
        rates = Num[_scaled_rate(w, e) for e in μ]
        _emit_eigen_rules!(rules, inv_rules, basis, comp, F.vectors, rates, pd, tt)
    else
        V, λ = s
        _emit_eigen_rules!(rules, inv_rules, basis, comp, V, Num[w * e for e in λ], pd, tt)
    end
    return nothing
end

# === Dressed frames ===

# Whether `_site_matrix` applies, without paying for the matrix or committing to its throws.
# The representative transition of a single-N-level-site `H0`, or `nothing`. `_site_matrix`
# is the throwing wrapper; both used to walk `H0` with the same three predicates.
# The representative collective operator of a single-collective-site `H0`, or `nothing`.
function _collective_site_op(H0::QAdd)
    σ = nothing
    for (term, _) in H0
        isempty(term.ops) && continue
        (length(term.ops) == 1 && is_collective_transition(term.ops[1]) && isempty(term.ne)) ||
            return nothing
        σ === nothing && (σ = term.ops[1])
        _site_key(term.ops[1]) == _site_key(σ) || return nothing
    end
    return σ
end

# The level matrix of a collective `H0`, over the levels it names. `n` comes from the highest
# level present, since the `Op` does not carry the space's level count.
function _collective_matrix(H0::QAdd, σ::Op)
    n = 0
    for (term, _) in H0
        isempty(term.ops) && continue
        n = max(n, Int(term.ops[1].l1), Int(term.ops[1].l2))
    end
    M = fill(Complex(_NUM_ZERO, _NUM_ZERO), n, n)
    for (term, c) in H0
        isempty(term.ops) && continue
        o = term.ops[1]
        M[Int(o.l1), Int(o.l2)] += to_num(c)
    end
    return M
end

function _nlevel_site_op(H0::QAdd)
    σ = nothing
    for (term, _) in H0
        isempty(term.ops) && continue
        (length(term.ops) == 1 && is_transition(term.ops[1]) && isempty(term.ne)) ||
            return nothing
        σ === nothing && (σ = term.ops[1])
        _site_key(term.ops[1]) == _site_key(σ) || return nothing
    end
    return σ
end

# Exactly zero for a symbolic entry, but a numeric one only has to be negligible against
# the level scale: conjugating an `H0` into its own dressed basis leaves rounding on the
# off-diagonal, and an exact test would send it back down the off-diagonal path.
function _off_diag_iszero(x::Complex{Num}, scale::Float64)
    (_iszero_num(real(x)) && _iszero_num(imag(x))) && return true
    v = _literal(real(x))
    w = _literal(imag(x))
    (v === nothing || w === nothing) && return false
    return abs(v) + abs(w) <= 1.0e-10 * max(scale, 1.0)
end

function _is_level_diagonal(M::Matrix{Complex{Num}})
    scale = 0.0
    for i in axes(M, 1)
        v = _literal(real(M[i, i]))
        v === nothing || (scale = max(scale, abs(v)))
    end
    return all(
        i -> all(j -> i == j || _off_diag_iszero(M[i, j], scale), axes(M, 2)),
        axes(M, 1),
    )
end

# The `n x n` matrix of an N-level `H0`, plus a representative operator of its site.
function _site_matrix(H0::QAdd)
    σ = nothing
    for (term, _) in H0
        isempty(term.ops) && continue
        (length(term.ops) == 1 && is_transition(term.ops[1])) || throw(
            ArgumentError(
                "`DressedFrame` needs an `H0` built from single N-level transitions; " *
                    "`$(_single_qadd(_CNUM_ONE, term.ops, term.ne))` is not one"
            )
        )
        o = term.ops[1]
        σ === nothing && (σ = o)
        _site_key(o) == _site_key(σ) || throw(
            ArgumentError(
                "`H0` spans more than one site; a dressed frame diagonalizes one site at " *
                    "a time"
            )
        )
    end
    σ === nothing && throw(ArgumentError("`H0` has no N-level term to diagonalize"))
    n = Int(σ.nlev)
    M = fill(Complex(_NUM_ZERO, _NUM_ZERO), n, n)
    for (term, c) in H0
        isempty(term.ops) && continue
        o = term.ops[1]
        M[Int(o.l1), Int(o.l2)] += to_num(c)
    end
    return σ, M
end

function _numeric_matrix(M::Matrix{Complex{Num}})
    n = size(M, 1)
    out = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n, j in 1:n
        c = _to_cnum(M[i, j])
        _is_native(c) || return nothing
        out[i, j] = _to_complex(c)
    end
    return out
end

# `eigen(Hermitian(·))` reads one triangle, so a non-Hermitian `H0` would be silently
# symmetrized into a different problem instead of rejected.
function _eigen_or_throw(Mn::Matrix{ComplexF64}, what::AbstractString)
    scale = max(1.0, maximum(abs, Mn))
    maximum(abs, Mn - Mn') <= 1.0e-10 * scale ||
        throw(ArgumentError("$what needs a Hermitian `H0`; its level matrix is not"))
    # A structured level matrix answers with exact eigenvectors, so a uniform chain or ring
    # dresses without rounding noise even when its energies are numeric.
    s = _structured_eigen(Mn)
    s === nothing || return s
    F = eigen(Hermitian(Mn))
    return F.vectors, Num[Num(e) for e in F.values]
end

# The two-level mixing angle `tan(2θ) = 2*m12/(m11 - m22)`.
function _mixing_angle(m11::Complex{Num}, m12::Complex{Num}, m21::Complex{Num}, m22::Complex{Num})
    (isequal(m12, m21) && _iszero_num(imag(m12))) || throw(
        ArgumentError(
            "`DressedFrame` needs a real, symmetric off-diagonal element, got " *
                "`$m12` and `$m21`"
        )
    )
    return atan(2 * real(m12), real(m11) - real(m22)) / 2
end

# Abel-Ruffini blocks a general symbolic eigendecomposition, but not a structured one: a
# block-diagonal `H0` diagonalizes block by block, and a two-level block has the closed form
# above. Only blocks of three or more need numbers.
# A symbolic block of three or more levels that is `d*I + w*C` with a numeric Hermitian `C`
# has *numeric* eigenvectors, and the eigenvector matrix is all a dressed basis needs. The
# energies `d + w*μ` stay symbolic. That covers a uniformly coupled chain or ring with a
# fully symbolic rate, which is what a symbolic eigendecomposition would otherwise be needed
# for. `_snap_spectrum!` runs for the same atom-identity reason as in the frame path.
function _dressed_block!(
        W::Matrix{Complex{Num}}, E::Vector{Num}, A::Matrix{CNum}, comp::Vector{Int},
        need_phases::Bool,
    )
    f = _common_factor_hermitian(A, comp)
    f === nothing && return false
    w, C, d = f
    # The shared diagonal is real: `_dressed_symbolic` checked the level energies first.
    dn = real(d)
    s = _structured_eigen(C)
    # See `_frame_multi_rule!`: a frame built on these needs phase cancellation, which a
    # radical entry forfeits. A static dressed basis has no phase and keeps the exact one.
    s === nothing || !need_phases || _all_rational(s[1]) || (s = nothing)
    if s === nothing
        F = eigen(Hermitian(C))
        V, μ = F.vectors, _snap_spectrum!(F.values)
        for (q, k) in enumerate(comp)
            E[k] = dn + _scaled_rate(w, μ[q])
            for (u, j) in enumerate(comp)
                W[j, k] = Complex(Num(real(V[u, q])), Num(imag(V[u, q])))
            end
        end
    else
        V, λ = s
        for (q, k) in enumerate(comp)
            E[k] = dn + w * λ[q]
            for (u, j) in enumerate(comp)
                W[j, k] = V[u, q]
            end
        end
    end
    return true
end

# The dressed basis `W` of a symbolic level matrix, its relations, and the dressed energies
# as `E`, or `nothing` for the latter when some block cannot state them in closed form. Only
# `RotatingFrame` needs the energies, and only for a matrix every block can answer for.
function _dressed_symbolic(M::Matrix{Complex{Num}}, need_phases::Bool)
    n = size(M, 1)
    all(i -> _iszero_num(imag(M[i, i])), 1:n) ||
        throw(ArgumentError("`H0` is not Hermitian: its level energies are not real"))
    A = CNum[_to_cnum(M[i, j]) for i in 1:n, j in 1:n]
    W = Complex{Num}[Complex(Num(i == j ? 1 : 0), _NUM_ZERO) for i in 1:n, j in 1:n]
    E = Num[real(M[i, i]) for i in 1:n]
    have_E = true
    rels = ParamRelation[]
    for comp in _ad_components(A)
        length(comp) == 1 && continue
        if length(comp) == 2
            # The two-level block goes first and stays exact: routing it through
            # `_dressed_block!` would trade its symbolic mixing angle for numeric
            # eigenvectors whenever the coupling ratio happened to fold to a literal.
            j, k = comp
            θ = _mixing_angle(M[j, j], M[j, k], M[k, j], M[k, k])
            W[j, j] = W[k, k] = Complex(cos(θ), _NUM_ZERO)
            W[j, k] = Complex(-sin(θ), _NUM_ZERO)
            W[k, j] = Complex(sin(θ), _NUM_ZERO)
            push!(rels, _trig_rel(θ))
            # Its dressed energies need a `sqrt` no relation closes, so this block cannot
            # supply them and a frame built on them would not certify.
            have_E = false
            continue
        end
        _dressed_block!(W, E, A, comp, need_phases) || throw(
            ArgumentError(
                "`DressedFrame` diagonalizes a symbolic block of three or more levels " *
                    "only when it is one real rate away from a numeric Hermitian matrix " *
                    "plus a shared diagonal, and levels $(join(comp, ", ")) are not. " *
                    "There is no general symbolic eigendecomposition to fall back on; " *
                    "give numeric level energies and couplings for this block"
            )
        )
    end
    return W, rels, (have_E ? E : nothing)
end

"""
    DressedFrame(H0::QAdd) -> UnitaryTransform

The static frame of the eigenstates of an N-level `H0`, i.e. `U` with `U'*H0*U` diagonal.

`H0` must be a Hermitian combination of transitions on a single site. With numeric level
energies and couplings any number of levels works, through an eigendecomposition. Symbolic
ones are diagonalized block by block, where a block is a set of levels `H0` couples: a
one-level block is already diagonal and a two-level one is the rotation by the mixing angle
`atan(2*Ω, Δ)/2`. A symbolic block of three or more levels works when it is one real rate
away from a numeric Hermitian matrix plus a shared diagonal, since its eigenvectors are then
numbers and its energies stay symbolic; a uniformly coupled chain or ring qualifies with both
its spacing and its coupling free. Anything else there needs numbers.

There is no gauge term, so this is a change of basis rather than a change of frame in time.
The dressed energies are `conjugate(H0, DressedFrame(H0))`.

That conjugation comes out diagonal in the numeric case. In the symbolic one the frame is
still exact, and `is_canonical` confirms it, but the off-diagonal element is not reduced to
zero: that needs `cos(θ)^2*Ω - sin(θ)^2*Ω - cos(θ)*sin(θ)*Δ = 0`, an identity between the
mixing angle and its own arguments, which neither relation kind of the reduction layer can
state.

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, 2); σ(i, j) = Transition(h, :σ, i, j);

julia> U = DressedFrame(0.5 * σ(1, 1) - 0.5 * σ(2, 2) + 0.4 * (σ(1, 2) + σ(2, 1)));

julia> is_canonical(U)
true
```

See also [`RotatingFrame`](@ref), [`Rotation`](@ref), [`transform`](@ref).
"""
function DressedFrame(H0::QAdd)
    σ, M = _site_matrix(H0)
    Mn = _numeric_matrix(M)
    Mn === nothing || return Rotation(σ, first(_eigen_or_throw(Mn, "`DressedFrame`")))
    W, rels, _ = _dressed_symbolic(M, false)
    U = Rotation(σ, W)
    # Same rule set, only the relations differ, so it needs no second validation.
    return _unchecked_transform(U.rules, U.inv_rules, U.gauge, U.time, rels, U.constraints)
end

function _level_diagonal(σ::Op, E::Vector{Num})
    d = QTermDict()
    for (k, e) in enumerate(E)
        c = _to_cnum(e)
        _iszero_cnum(c) && continue
        _addto!(d, Op[_transition_op(σ, k, k)], c)
    end
    return QAdd(d, _EMPTY_INDICES)
end

# `exp(-im*θ*G) = V*exp(-im*θ*Λ)*V'`: the dressed basis change wrapped around the frame of
# the dressed energies, built by composition. See "A single N-level site is exactly solvable"
# in the devdocs. Every factor here is static, so the composition never touches a gauge and
# the caller installs the one belonging to `G`.
function _dressed_rules(σ::Op, M::Matrix{Complex{Num}}, θ::Num)
    Mn = _numeric_matrix(M)
    if Mn === nothing
        # Symbolic, but a block that splits as `d*I + w*C` still has numeric eigenvectors
        # and symbolic energies `d + w*μ`, which is everything this needs.
        W, _, E = _dressed_symbolic(M, true)
        E === nothing && throw(
            ArgumentError(
                "an off-diagonal N-level generator with symbolic entries needs every " *
                    "coupled block to be one real rate away from a numeric Hermitian " *
                    "matrix; a two-level block is not, since its dressed energies need a " *
                    "`sqrt` the coefficient reduction cannot close. Give numeric level " *
                    "energies and couplings"
            )
        )
        return _dressed_compose(σ, W, E, θ)
    end
    V, E = _eigen_or_throw(Mn, "`RotatingFrame`")
    return _dressed_compose(σ, V, E, θ)
end

function _dressed_compose(σ::Op, W::AbstractMatrix, E::Vector{Num}, θ::Num)
    D = Rotation(σ, W)
    U = D * _static(_exponentiate(_level_diagonal(σ, E), θ)...) * inv(D)
    return (U.rules, U.inv_rules, U.reductions)
end

"""
    RotatingFrame(H0::QAdd, t) -> UnitaryTransform

The frame `U = exp(-im*H0*t)`, derived from the adjoint action of `H0`.

`H0` must be Hermitian and time-independent, and `[H0, ·]` must close on the generators of
the sites it touches. A generator that only picks up a multiple of itself gets a phase; a
pair that mixes gets the closed-form rotation or squeeze of its `2x2` block. That covers a
number-diagonal `H0` written in any form, a spin or Pauli `Sz`, a beamsplitter
`g*(a'*b + b'*a)`, a harmonic `(x^2 + p^2)/2`, and a parametric `a^2 + a'^2`. For a moving
`H0 = f(t)*G` use `UnitaryTransform(G, θ, t)` with `θ` the integral of `f`.

A drive makes the adjoint action affine rather than linear, and that is covered too: the
generators are shifted by the c-numbers that homogenize it, so `η*(a + a')` gives a
displacement linear in `t` and `ω*a'*a + η*(a + a')` the displaced rotation. A drive with no
such shift is resonant, grows without bound, and throws.

Any number of generators may mix, on either of two closed forms. A block that squares to a
multiple of the identity is exact and keeps its rate symbolic; otherwise the block must be
one real rate away from a numeric Hermitian matrix plus a diagonal they all share, and is
diagonalized. A hopping chain or ring of identical modes qualifies with a fully symbolic
frequency and coupling, and so does a spin rotation about an arbitrary axis. Where the
eigenbasis is numeric those frames cancel to rounding rather than exactly.

A single N-level site closes whatever the shape of `H0`, so an off-diagonal one is covered
too, on any number of levels: it factors as [`DressedFrame`](@ref) around the diagonal
frame of the dressed energies. Its eigenvectors have to be numbers, which they are whenever
the level matrix is numeric, and also when a symbolic block is one rate away from a numeric
Hermitian one.

A [`CollectiveTransition`](@ref) site is covered on the level diagonal, with one phase per
level as below. Off the diagonal there is no dressed factorization to fall back on, and it
throws.

An operator family carrying a free index is covered as a whole: `Rotation(a_i, θ)` is the
same local transform at every site of `i`, and [`conjugate`](@ref) instantiates it wherever
it meets the family, inside a sum or out of it. `RotatingFrame(Σ(ω*a_i'*a_i, i), t)` is such
a family. A per-slot index like `i(3)` names one resolved site instead.

An N-level frame carries one phase per level rather than one per transition, so the
commensurate frequencies of three or more levels (`E₁-E₃ = (E₁-E₂) + (E₂-E₃)`) cancel as
exponent arithmetic.

The gauge term is exactly `-H0`, so `transform(H0, RotatingFrame(H0, t))` is zero. Where the
eigenbasis is numeric it is zero only to rounding, and that rounding scales with the
parameters of `H0`, so no tolerance on the coefficients alone can certify it.

A commutator that leaves the linear span of the generators puts the frame outside the
exactly solvable class; it throws and names the term responsible. That is what rules out
the mixed atom-cavity frames: `[Δ*σ22 + g*(a'*σ12 + a*σ21), σ12]` produces `a*σ22`, and no
finite set of polynomial generators is closed under this `ad`.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables ω η t;

julia> transform(ω * a' * a + η * (a + a'), RotatingFrame(ω * a' * a, t))
exp(-im*t*ω)*η * a + exp(im*t*ω)*η * a'
```

See also [`transform`](@ref), [`gauge_term`](@ref).
"""
# === Affine adjoint action ===
#
# A drive makes `[G, ·]` affine rather than linear: `[η*(a + a'), a] = -η`. Shifting the
# generators by c-numbers absorbs that, since `g_k -> g_k + s_k` is homogeneous exactly when
# `Σ_j A[j,k] s_j = c_k`. The image is then the homogeneous one plus the constant the shift
# leaves behind, which needs no separate exponential: it is read off the rules already built.

_cnum_div(a::CNum, b::CNum) = _to_cnum(to_num(a) / to_num(b))

# `Aᵀ s = c` on one component, or `nothing` when `c` is outside the range of `Aᵀ` and the
# drive is therefore secular. Any particular solution will do: two differ by a kernel
# element, which `M` fixes, so the image below does not depend on the choice.
function _solve_shift(A::Matrix{CNum}, c::Vector{CNum}, comp::Vector{Int})
    m = length(comp)
    M = Matrix{CNum}(undef, m, m + 1)
    for u in 1:m
        for v in 1:m
            M[u, v] = A[comp[v], comp[u]]
        end
        M[u, m + 1] = c[comp[u]]
    end
    for k in 1:m
        # A literal pivot is preferred over a symbolic one, which could be secretly zero.
        pivot = 0
        for u in k:m
            _iszero_cnum(M[u, k]) && continue
            pivot == 0 && (pivot = u)
            if _is_native(M[u, k])
                pivot = u
                break
            end
        end
        pivot == 0 && continue
        if pivot != k
            for v in 1:(m + 1)
                M[k, v], M[pivot, v] = M[pivot, v], M[k, v]
            end
        end
        for u in 1:m
            (u == k || _iszero_cnum(M[u, k])) && continue
            f = _cnum_div(M[u, k], M[k, k])
            for v in k:(m + 1)
                M[u, v] = _add_cnum(M[u, v], _neg_cnum(_mul_cnum(f, M[k, v])))
            end
        end
    end
    s = fill(_CNUM_ZERO, m)
    for u in 1:m
        k = findfirst(v -> !_iszero_cnum(M[u, v]), 1:m)
        # A row that eliminated to `0 = rhs` is consistent only if the rhs went too.
        k === nothing && (_iszero_cnum(M[u, m + 1]) ? continue : return nothing)
        s[k] = _cnum_div(M[u, m + 1], M[u, k])
    end
    return s
end

# The coefficient of `g` in an already-built rule image. Every image is linear in the
# generators by construction, so the homogeneous matrix is read back rather than rebuilt,
# which is what keeps the forward and inverse shifts in step with their own directions.
function _image_coeff(r::QAdd, g::Op)
    for (term, coeff) in r
        (length(term.ops) == 1 && isempty(term.ne) && isequal(term.ops[1], g)) &&
            return coeff
    end
    return _CNUM_ZERO
end

_shift_constant(image::QAdd, basis::Vector{Op}, comp::Vector{Int}, s::Vector{CNum}, u::Int) =
    _add_cnum(
    sum(
        v -> _mul_cnum(_image_coeff(image, basis[comp[v]]), s[v]),
        eachindex(comp); init = _CNUM_ZERO,
    ),
    _neg_cnum(s[u]),
)

_add_constant!(rules::Dict{Op, QAdd}, g::Op, δ::CNum) =
    _iszero_cnum(δ) || (rules[g] = rules[g] + _rule_qadd((δ, Op[])))

function _apply_affine_shifts!(
        rules::Dict{Op, QAdd}, inv_rules::Dict{Op, QAdd}, basis::Vector{Op},
        A::Matrix{CNum}, c::Vector{CNum}, comps::Vector{Vector{Int}}, θ::Num,
    )
    for comp in comps
        # Every existing frame lands here, so the pass costs one scan of `c` and no more.
        all(k -> _iszero_cnum(c[k]), comp) && continue
        if all(u -> all(v -> _iszero_cnum(A[u, v]), comp), comp)
            # An undriven-but-unrotated component: `dg̃/dθ = im*c` integrates to a shift
            # linear in `θ`, which is the displacement a resonant drive produces.
            for k in comp
                δ = _mul_cnum(_CNUM_IM, _mul_cnum(c[k], _to_cnum(θ)))
                _add_constant!(rules, basis[k], δ)
                _add_constant!(inv_rules, basis[k], _neg_cnum(δ))
            end
            continue
        end
        s = _solve_shift(A, c, comp)
        s === nothing && throw(
            ArgumentError(
                "`H0` drives $(join(string.(basis[comp]), ", ")) resonantly: the " *
                    "inhomogeneous part of `[H0, ·]` is outside the range of its " *
                    "homogeneous part, so the frame grows without bound and has no closed " *
                    "form here. Split `H0`, or build the transform explicitly"
            )
        )
        for (u, k) in enumerate(comp)
            _add_constant!(rules, basis[k], _shift_constant(rules[basis[k]], basis, comp, s, u))
            _add_constant!(
                inv_rules, basis[k],
                _shift_constant(inv_rules[basis[k]], basis, comp, s, u),
            )
        end
    end
    return nothing
end

# The scalar `exp(-im*θ*G)` is exponentiated against is arbitrary: a frame passes the time
# variable, a static transform an angle. Nothing below differentiates it, so the two differ
# only in the gauge their callers install.
function _exponentiate(G::QAdd, θ::Num)
    if _nlevel_site_op(G) !== nothing
        σ, M = _site_matrix(G)
        _is_level_diagonal(M) || return _dressed_rules(σ, M, θ)
        return (_level_phase_rules(σ, M, θ)..., ParamRelation[])
    end
    σc = _collective_site_op(G)
    if σc !== nothing
        # One phase per level, not one per transition: `S^{ij}` rotates at `E_i - E_j`, and
        # phasing the pairs would mint an atom per difference, leaving `E₁-E₃` invisible.
        M = _collective_matrix(G, σc)
        _is_level_diagonal(M) || throw(
            ArgumentError(
                "an off-diagonal collective `H0` is not supported: the dressed " *
                    "factorization behind the N-level case has no counterpart here. " *
                    "Diagonalize it yourself and pass the eigenbasis to `Rotation(S, W)`"
            )
        )
        rules, inv_rules = _level_phase_rules(
            M, size(M, 1), (i, j) -> _collective_op(σc, i, j), θ,
        )
        return (rules, inv_rules, ParamRelation[])
    end
    basis = _frame_basis(G)
    A, c = _ad_affine(G, basis)
    rules = Dict{Op, QAdd}()
    inv_rules = Dict{Op, QAdd}()
    rels = ParamRelation[]
    comps = _ad_components(A)
    for comp in comps
        # An adjoint-related component is derived rather than exponentiated again: two
        # independent runs could pick opposite rate signs and break Hermiticity.
        all(j -> haskey(rules, basis[j]), comp) && continue
        if length(comp) == 1
            _frame_diagonal_rule!(rules, inv_rules, basis[comp[1]], A[comp[1], comp[1]], θ)
            # The exact `2x2` form goes first where it applies: it reads the rate straight off
            # the block rather than through a numeric eigendecomposition, and its sign
            # orientation is what the printed two-mode forms are pinned to. The involutive
            # path generalizes it to any size and any scale; `eigen` is the last resort.
        elseif !(
                length(comp) == 2 &&
                    _frame_block_rule!(rules, inv_rules, rels, basis, A, comp[1], comp[2], θ)
            ) && !_frame_involutive_rule!(rules, inv_rules, rels, basis, A, comp, θ)
            _frame_multi_rule!(rules, inv_rules, basis, A, comp, θ)
        end
        for j in comp
            ga = adjoint(basis[j])
            haskey(rules, ga) && continue
            rules[ga] = adjoint(rules[basis[j]])
            inv_rules[ga] = adjoint(inv_rules[basis[j]])
        end
    end
    # After the adjoint derivation, never before: that block copies a rule wholesale, so a
    # shift added earlier would be counted twice for one of each adjoint pair.
    _apply_affine_shifts!(rules, inv_rules, basis, A, c, comps, θ)
    return (rules, inv_rules, rels)
end

"""
    UnitaryTransform(G::QAdd, θ) -> UnitaryTransform
    UnitaryTransform(G::QAdd, θ, t) -> UnitaryTransform

The transform `U = exp(-im*θ*G)`, derived from the adjoint action of `G`.

This is the generic constructor behind [`Rotation`](@ref), [`Squeeze`](@ref) and
[`RotatingFrame`](@ref): any `G` whose commutator closes on the generators of the sites it
touches is admissible, on the same terms [`RotatingFrame`](@ref) documents. Use it for a
generator no named constructor covers.

The two-argument form is static and its gauge term is zero. The three-argument form declares
`θ` a function of `t`, whose gauge term is then `-∂ₜθ*G`; that is exact because `G` itself is
time-independent. `RotatingFrame(H0, t)` is the case `θ = t`.

# Examples

```jldoctest
julia> h = FockSpace(:f); @qnumbers a::Destroy(h);

julia> @variables θ;

julia> conjugate(a, UnitaryTransform(a' * a, θ))
exp(-im*θ) * a
```

See also [`RotatingFrame`](@ref), [`conjugate`](@ref), [`is_canonical`](@ref).
"""
function UnitaryTransform(G::QAdd, θ::_ScalarLike)
    rules, inv_rules, rels = _exponentiate(G, _as_num(θ))
    return UnitaryTransform(rules, inv_rules, _zero_qadd(), _NUM_ZERO, rels, Num[])
end

function UnitaryTransform(G::QAdd, θ::_ScalarLike, t::_ScalarLike)
    tt = _as_num(t)
    return _with_gauge(UnitaryTransform(G, θ), _gauge(G, θ, tt), tt)
end

# `exp(-im*H0*t)` is the propagator of `H0` only for a time-independent one, and the gauge
# `-H0` assumes it too. Without this the wrong answer certifies: `RotatingFrame(ω*cos(t)*a'a)`
# returned `exp(-im*cos(t)*t*ω)*a` where the phase is `exp(-im*ω*sin(t))`.
function _check_static_H0(H0::QAdd, tt::Num)
    for (_, c) in H0
        _depends_on_time(c, tt) && throw(
            ArgumentError(
                "`RotatingFrame` needs a time-independent `H0`: `exp(-im*H0*t)` is not " *
                    "the propagator of a moving `H0`, and its gauge is not `-H0`. For " *
                    "`H0 = f($tt)*G` use `UnitaryTransform(G, θ, $tt)` with `θ` the " *
                    "integral of `f`"
            )
        )
    end
    return nothing
end

# `im*(∂ₜU')U = im*(im*H0)*U'U = -H0` for any time-independent `H0`, with no per-shape gauge
# table and no assumption that per-site pieces add up. That is `_gauge(H0, t, t)`, since
# `∂ₜt = 1`, so the frame is the `θ = t` case of the generic constructor.
function RotatingFrame(H0::QAdd, t::_ScalarLike)
    tt = _as_num(t)
    _check_static_H0(H0, tt)
    return UnitaryTransform(H0, tt, tt)
end
