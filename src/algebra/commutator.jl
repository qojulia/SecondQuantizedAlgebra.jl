"""
    commutator(a, b)

Commutator of `a` and `b`.

Returns ``[a, b] = ab - ba`` in canonical form. If either argument is a formal
operator function, the result stays formal until lowered with [`taylor`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> commutator(a, a')
1

julia> commutator(a' * a, a)
-a
```

See also [`anticommutator`](@ref), [`normal_order`](@ref), [`taylor`](@ref).
"""
function commutator end

commutator(::Number, ::Number) = zero_qadd()
commutator(::Number, ::QField) = zero_qadd()
commutator(::QField, ::Number) = zero_qadd()

function commutator(a::QSym, b::QSym)
    site_compare(a, b, EMPTY_NE) === Equal || return zero_qadd()
    isequal(a, b) && return zero_qadd()
    forward = can_commute(a, b)
    reverse = can_commute(b, a)
    if !forward && reverse
        c1, r1, c2, r2, n = commute_pair(a, b)
    elseif forward && !reverse
        c1, r1, c2, r2, n = commute_pair(b, a)
        c1, c2 = -c1, -c2
    else
        return a * b - b * a
    end
    d = QTermDict()
    iszero(c1) || addto!(d, n == 0 ? EMPTY_OPS : Op[r1], c1)
    iszero(c2) || addto!(d, n == 0 ? EMPTY_OPS : Op[r2], c2)
    return QAdd(d, EMPTY_INDICES)
end

function commutator(a::QAdd, b::QSym)
    term = QTerm(Op[b], EMPTY_NE)
    return commute_terms(a, [(term, CNUM_ONE, acts_on(term))], EMPTY_INDICES) do _, _
        b
    end
end
commutator(a::QSym, b::QAdd) = -commutator(b, a)

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return zero_qadd()
    bside = [(tb, cb, acts_on(tb)) for (tb, cb) in b.arguments]
    return commute_terms(a, bside, b.indices) do tb, cb
        QAdd(QTermDict(copy_key(tb) => cb), b.indices)
    end
end

function commute_terms(
        bfactor::F, a::QAdd, bside::Vector{Tuple{QTerm, CNum, Vector{Int}}},
        b_indices::Vector{Index},
    ) where {F}
    indexed = !(isempty(a.indices) && isempty(b_indices))
    d = QTermDict()
    indices = EMPTY_INDICES
    for (ta, ca) in a.arguments
        aon_a = acts_on(ta)
        qa = indexed ? QAdd(QTermDict(copy_key(ta) => ca), a.indices) : ZERO_QADD
        for (tb, cb, aon_b) in bside
            isdisjoint(aon_a, aon_b) && continue
            if indexed
                qb = bfactor(tb, cb)
                indices = merge_unique(indices, merge_into!(d, qa * qb - qb * qa))
            else
                emit_product!(d, ta.ops, ca, ta.ne, tb.ops, cb, tb.ne, EMPTY_INDICES, false)
                emit_product!(d, tb.ops, -cb, tb.ne, ta.ops, ca, ta.ne, EMPTY_INDICES, false)
            end
        end
    end
    return QAdd(d, drop_unused_indices(d, indices))
end

"""
    anticommutator(a, b)

Anticommutator of `a` and `b`.

Returns ``\\{a, b\\} = ab + ba`` in canonical form.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> anticommutator(a, a')
1 + 2 * a' * a
```

See also [`commutator`](@ref).
"""
anticommutator(a, b) = a * b + b * a
