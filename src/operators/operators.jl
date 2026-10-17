function Base.adjoint(o::Op)
    k = o.kind
    if k === OP_DESTROY
        return with_kind(o, OP_CREATE)
    elseif k === OP_CREATE
        return with_kind(o, OP_DESTROY)
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        return with_levels(o, o.l2, o.l1)
    else
        return o
    end
end

"""
    order_key(op::Op) -> Tuple

Total ordering key of an operator that ties two operators exactly when they are `isequal`.

The tuple is `(space_index, kind, name rank, index, l1, l2, g, nlev)`: operators group by
subspace, then by kind, then alphabetically by name, and the packed level and axis fields
keep the remaining operators apart. Unlike `hash`, the key is reproducible across sessions
and Julia versions.
"""
order_key(o::Op) = (o.space_index, Int(o.kind), name_rank(o.name_id), index_key(o.index), Int(o.l1), Int(o.l2), Int(o.g), Int(o.nlev))

"""
    canonical_isless(a, b) -> Bool

Whether `a` precedes `b` in the canonical order of operators, operator sums, operator
functions, superoperators or dissipation channels.

The order is total and reproducible: for two objects of the same type, exactly one of
`canonical_isless(a, b)`, `canonical_isless(b, a)` and `isequal(a, b)` holds, in every
session and Julia version. Use it to pick a canonical representative, e.g. of an operator
and its adjoint, or to sort with `sort(v; lt = canonical_isless)`.

It is not a method of `Base.isless` on purpose. `a + b < 2a` would read as a comparison of
operator values, and `sort`, `max` and `<` would silently use an order that has no physical
meaning.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: canonical_isless

julia> h = FockSpace(:cavity); a = Destroy(h, :a); b = Destroy(h, :b);

julia> canonical_isless(a, b), canonical_isless(b, a)
(true, false)

julia> sort([b, a]; lt = canonical_isless)
2-element Vector{Op}:
 a
 b
```

See also [`QAdd`](@ref), [`channels`](@ref).
"""
canonical_isless(a::Op, b::Op) = isless(order_key(a), order_key(b))

const SiteKey = Tuple{Int32, Index, Int32}

site_key(o::Op) =
    (o.space_index, o.index, (is_position(o) || is_momentum(o)) ? Int32(0) : o.name_id)

@inline site_family(k::OpKind) =
    (k === OP_DESTROY || k === OP_CREATE) ? 0x00 :
    (k === OP_POSITION || k === OP_MOMENTUM) ? 0x05 :
    UInt8(k)

function site_compare(a::Op, b::Op, ne::Vector{NonEqualPair})
    ka = a.kind; kb = b.kind
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    fa = site_family(ka)
    fa === site_family(kb) || return UInt8(ka) < UInt8(kb) ? Less : Greater
    if !(fa === 0x05 && ka !== kb)
        a.name_id == b.name_id || return name_rank(a.name_id) < name_rank(b.name_id) ? Less : Greater
    end
    a.index == b.index && return Equal
    fa === 0x05 && return Undetermined
    ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

function can_commute(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return false
    elseif ka === OP_TRANSITION && kb === OP_TRANSITION
        return false
    elseif ka === OP_PAULI && kb === OP_PAULI
        return false
    elseif ka === OP_SPIN && kb === OP_SPIN
        return a.l1 <= b.l1
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        a.name_id == b.name_id || return true
        return a.l1 > b.l1 || (a.l1 == b.l1 && a.l2 >= b.l2)
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return false
    else
        return true
    end
end

function commute_pair(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return (CNUM_ONE, a, CNUM_ZERO, a, 0)
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return (CNUM_NEG_IM, a, CNUM_ZERO, a, 0)
    elseif ka === OP_SPIN && kb === OP_SPIN
        contracted, c = axis_product(a, b)
        return (c, contracted, CNUM_ZERO, contracted, 1)
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        c1 = a.l2 == b.l1 ? CNUM_ONE : CNUM_ZERO
        c2 = b.l2 == a.l1 ? CNUM_NEG1 : CNUM_ZERO
        return (c1, with_levels(a, a.l1, b.l2), c2, with_levels(a, b.l1, a.l2), 1)
    else
        return (CNUM_ZERO, a, CNUM_ZERO, a, 0)
    end
end

const levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

@inline levi_civita_coeff(eps::Int) = eps == 1 ? CNUM_IM : (eps == -1 ? CNUM_NEG_IM : CNUM_ZERO)

@inline axis_product(a::Op, b::Op) =
    (with_levels(a, 6 - a.l1 - b.l1, 0), levi_civita_coeff(levi_civita[a.l1][b.l1]))

@inline same_site(a::Op, b::Op) =
    a.name_id == b.name_id && a.space_index == b.space_index && a.index == b.index

may_reduce(a::Op, b::Op) =
    (a.kind === OP_TRANSITION && b.kind === OP_TRANSITION) ||
    (a.kind === OP_PAULI && b.kind === OP_PAULI)

function reduce_pair(a::Op, b::Op)
    ka = a.kind
    if ka === OP_TRANSITION && b.kind === OP_TRANSITION
        same_site(a, b) || return (NoReduction, a, CNUM_ZERO)
        a.l2 == b.l1 || return (ScalarReduction, a, CNUM_ZERO)
        return (OpReduction, with_levels(a, a.l1, b.l2), CNUM_ONE)
    elseif ka === OP_PAULI && b.kind === OP_PAULI
        same_site(a, b) || return (NoReduction, a, CNUM_ZERO)
        a.l1 == b.l1 && return (ScalarReduction, a, CNUM_ONE)
        product, c = axis_product(a, b)
        return (OpReduction, product, c)
    else
        return (NoReduction, a, CNUM_ZERO)
    end
end

function ground_state_expand(op::Op)
    op.kind === OP_TRANSITION || return (false, 0, 0, 0)
    (op.l1 == op.g && op.l2 == op.g) || return (false, 0, 0, 0)
    return (true, Int(op.g), Int(op.nlev), Int(op.space_index))
end
