"""
    SuperTerm

One sandwich `ρ ↦ left * ρ * right` of a [`SuperOperator`](@ref), with the pairwise index
inequality constraints shared by both sides.

# Fields
- `left::Vector{Op}`: the ordered operator product acting from the left
- `right::Vector{Op}`: the ordered operator product acting from the right
- `ne::Vector{NonEqualPair}`: pairwise `α ≠ β` constraints (canonicalized)

`SuperTerm` is the dict key of a `SuperOperator`; iterating one yields
`Pair{SuperTerm, CNum}`. The fields are owned by the stored key and must not be mutated.
Build keys with [`superterm`](@ref).
"""
struct SuperTerm
    left::Vector{Op}
    right::Vector{Op}
    ne::Vector{NonEqualPair}
    hash::UInt
    SuperTerm(left::Vector{Op}, right::Vector{Op}, ne::Vector{NonEqualPair}) =
        new(left, right, ne, superterm_hash(left, right, ne))
    SuperTerm(left::Vector{Op}, right::Vector{Op}, ne::Vector{NonEqualPair}, h::UInt) =
        new(left, right, ne, h)
end

@inline superterm_hash(left::Vector{Op}, right::Vector{Op}, ne::Vector{NonEqualPair}) =
    hash(:SuperTerm, hash(left, hash(right, hash(ne, zero(UInt)))))

Base.isequal(a::SuperTerm, b::SuperTerm) = a.hash == b.hash && isequal(a.left, b.left) &&
    isequal(a.right, b.right) && isequal(a.ne, b.ne)
Base.:(==)(a::SuperTerm, b::SuperTerm) = isequal(a, b)
Base.hash(t::SuperTerm, h::UInt) = hash(t.hash, h)

const SuperTermDict = Dict{SuperTerm, CNum}

copy_key(t::SuperTerm) = SuperTerm(copy(t.left), copy(t.right), copy_ne(t.ne), t.hash)

"""
    superterm(left, right, ne = Tuple{Index, Index}[]) -> SuperTerm

Build the storage key of the sandwich `ρ ↦ left * ρ * right`. The operator vectors are
copied and `ne` is canonicalized; the sides are used as given, so a lookup with
`get(S, superterm(...), 0)` matches only sides already in canonical order, such as those
read from an iterated term.
"""
superterm(left::AbstractVector{<:QSym}, right::AbstractVector{<:QSym}, ne::Vector{NonEqualPair} = EMPTY_NE) =
    SuperTerm(Op[o for o in left], Op[o for o in right], canonical_ne(ne))

function addto_key!(d::SuperTermDict, term::SuperTerm, c::CNum)
    existing = get(d, term, CNUM_ZERO)
    if existing === CNUM_ZERO
        iszero_cnum(c) && return d
        d[term] = c
        return d
    end
    new_c = add_cnum(existing, c)
    if iszero_cnum(new_c)
        delete!(d, term)
    else
        d[term] = new_c
    end
    return d
end

function copy_args(d::SuperTermDict)
    out = SuperTermDict()
    sizehint!(out, length(d))
    for (term, c) in d
        out[copy_key(term)] = c
    end
    return out
end

# Joint dependence of a sandwich on `idx`, without concatenating the two sides.
function depends_on_sandwich(c::CNum, left::Vector{Op}, right::Vector{Op}, idx::Index)
    for op in right
        op.index == idx && return true
    end
    return depends_on_index_term(c, left, idx)
end

function depends_on_scope_sandwich(
        c::CNum, left::Vector{Op}, right::Vector{Op}, ne::Vector{NonEqualPair}, idx::Index,
    )
    depends_on_sandwich(c, left, right, idx) && return true
    for (α, β) in ne
        (α == idx || β == idx) && return true
    end
    return false
end

bound_indices(t::SuperTerm, c::CNum, indices::Vector{Index}) =
    Index[idx for idx in indices if depends_on_scope_sandwich(c, t.left, t.right, t.ne, idx)]

@inline function pair_referenced(
        p::NonEqualPair, t::SuperTerm, c::CNum, scope::Vector{Index},
    )
    α, β = p
    α in scope && return true
    β in scope && return true
    depends_on_sandwich(c, t.left, t.right, α) && return true
    depends_on_sandwich(c, t.left, t.right, β) && return true
    return false
end

function prune_dead_ne(args::SuperTermDict, indices::Vector{Index})
    dead((term, c)) = any(p -> !pair_referenced(p, term, c, indices), term.ne)
    any(dead, args) || return args
    out = SuperTermDict()
    sizehint!(out, length(args))
    for (term, c) in args
        kept = NonEqualPair[]
        for p in term.ne
            pair_referenced(p, term, c, indices) && push_ne_unique!(kept, p)
        end
        kept_ne = isempty(kept) ? EMPTY_NE : insertion_sort!(kept, isless)
        addto_key!(out, SuperTerm(copy(term.left), copy(term.right), kept_ne), c)
    end
    return out
end

function drop_unused_indices(d::SuperTermDict, indices::Vector{Index})
    isempty(indices) && return indices
    used = Index[]
    for idx in indices
        for (term, c) in d
            if depends_on_scope_sandwich(c, term.left, term.right, term.ne, idx)
                push!(used, idx)
                break
            end
        end
    end
    length(used) == length(indices) && return indices
    return used
end

function constraint_pairs(args::SuperTermDict)
    out = NonEqualPair[]
    for term in keys(args), p in term.ne
        push_ne_unique!(out, p)
    end
    return out
end

"""
    canonical_sandwich!(out, left, right, c, ne)

Canonicalize each side under the shared constraints and add the outer product of the two
canonical sides into `out`. Both sides are copied.
"""
function canonical_sandwich!(
        out::SuperTermDict, left::Vector{Op}, right::Vector{Op}, c::CNum,
        ne::Vector{NonEqualPair},
    )
    iszero_cnum(c) && return out
    cne = canonical_ne(ne)
    lt = QTermDict()
    canonicalize!(lt, copy(left), CNUM_ONE, cne)
    rt = QTermDict()
    canonicalize!(rt, copy(right), CNUM_ONE, cne)
    for (a, ca) in lt, (b, cb) in rt
        addto_key!(out, SuperTerm(a.ops, b.ops, cne), mul_cnum(c, mul_cnum(ca, cb)))
    end
    return out
end

"""
    settle_scope!(out, temp::SuperTermDict, scope)

The two-sided analogue of the `QAdd` method: every scope index that a canonical sandwich of
`temp` no longer mentions through either side or its coefficient is eliminated with its
constrained count, or kept with its constraints when the count is unknown.
"""
function settle_scope!(out::SuperTermDict, temp::SuperTermDict, scope::Vector{Index})
    for (term, cv) in temp
        mentions(c, idx) = depends_on_sandwich(c, term.left, term.right, idx)
        coef, ne = eliminate_binders(mentions, cv, term.ne, scope)
        key = ne === term.ne ? term : SuperTerm(term.left, term.right, ne)
        addto_key!(out, key, coef)
    end
    return out
end

function settle_sandwich!(
        out::SuperTermDict, left::Vector{Op}, right::Vector{Op}, c::CNum,
        ne::Vector{NonEqualPair}, scope::Vector{Index},
    )
    isempty(scope) && return canonical_sandwich!(out, left, right, c, ne)
    temp = SuperTermDict()
    canonical_sandwich!(temp, left, right, c, ne)
    return settle_scope!(out, temp, scope)
end

# The indices of `ops` that the summation index `i` may coincide with under `ne`.
coinciding_pairs(i::Index, ops::Vector{Op}, ne::Vector{NonEqualPair}) =
    Tuple{Index, Index}[(i, ext) for ext in distinct_op_indices(ops) if may_coincide(i, ext, ne)]

"""
    emit_sandwich!(out, left, right, c, ne, scope)

Emit one joint sandwich bound by `scope` (rules I5 and I6): partition the bound indices by
their coincidences with the indices of either side ([`partition_binders!`](@ref)), with the
diagonal substitution applied to the unsorted sides, then canonicalize both sides and count
every binder the result no longer mentions.
"""
function emit_sandwich!(
        out::SuperTermDict, left::Vector{Op}, right::Vector{Op}, c::CNum,
        ne::Vector{NonEqualPair}, scope::Vector{Index},
    )
    isempty(scope) && return canonical_sandwich!(out, left, right, c, ne)
    partition_binders!((left, right), c, ne, scope, scope) do (l, r), pc, pne, pscope, _
        settle_sandwich!(out, l, r, pc, pne, pscope)
    end
    return out
end

"""
    emit_applied!(out::QTermDict, ops, c, ne, scope)

The one-sided version of [`emit_sandwich!`](@ref), used when a map is applied to an operator.
"""
function emit_applied!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
        scope::Vector{Index},
    )
    isempty(scope) && return canonicalize!(out, ops, c, ne)
    partition_binders!((ops,), c, ne, scope, scope) do (part,), pc, pne, pscope, _
        temp = QTermDict()
        canonicalize!(temp, copy(part), pc, pne)
        settle_scope!(out, temp, pscope)
    end
    return out
end

reversed_adjoint(ops::Vector{Op}) = Op[adjoint(o) for o in Iterators.reverse(ops)]
