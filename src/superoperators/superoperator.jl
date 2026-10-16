"""Marker for the trusted `SuperOperator` constructor: the terms are already canonical."""
struct CanonicalTerms end

"""
    SuperOperator
    SuperOperator(L::Lindbladian)
    SuperOperator(pairs, indices)

A polynomial linear map on operators.

The map acts as `ρ ↦ Σ c * A * ρ * B`. Build maps with [`sprepost`](@ref), [`spre`](@ref) and [`spost`](@ref), or lower a
[`Lindbladian`](@ref) with `SuperOperator(L)`. `SuperOperator(pairs, indices)` rebuilds a
map from `SuperTerm => coefficient` pairs and a summation scope. Maps compose with `∘` or
`*`, the right factor acting first, and `S(ρ)` applies a map. `transpose(S)` satisfies
`Tr(A * S(ρ)) == Tr(transpose(S)(A) * ρ)` and is the Heisenberg-picture map, while
`adjoint(S)` is the Hilbert–Schmidt adjoint. Iterating a map gives its
`SuperTerm => coefficient` pairs.

A map is not an operator, so it does not enter operator products or averages. The name is
public but not exported, since the numeric backends export a `SuperOperator` of their own.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: SuperOperator

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> SuperOperator(dissipator(a))
-1//2 * (·) ⋅ a' * a + a ⋅ (·) ⋅ a' - 1//2 * a' * a ⋅ (·)
```

See also [`Lindbladian`](@ref), [`SuperExpr`](@ref), [`compose`](@ref).
"""
struct SuperOperator
    arguments::SuperTermDict
    indices::Vector{Index}
    function SuperOperator(::CanonicalTerms, arguments::SuperTermDict, indices::Vector{Index})
        return new(prune_dead_ne(arguments, indices), indices)
    end
end

superop(args::SuperTermDict, indices::Vector{Index}) =
    SuperOperator(CanonicalTerms(), args, indices)
superop_pruned(args::SuperTermDict, indices::Vector{Index}) =
    superop(args, drop_unused_indices(args, indices))

function SuperOperator(pairs, indices::Vector{Index})
    out = SuperTermDict()
    for (t, c) in pairs
        cc = to_cnum(c)
        emit_sandwich!(out, t.left, t.right, cc, t.ne, bound_indices(t, cc, indices))
    end
    return superop_pruned(out, copy(indices))
end

Base.zero(::Type{SuperOperator}) = superop(SuperTermDict(), Index[])
Base.zero(::SuperOperator) = zero(SuperOperator)
Base.one(::Type{SuperOperator}) =
    superop(SuperTermDict(SuperTerm(Op[], Op[], EMPTY_NE) => CNUM_ONE), Index[])
Base.one(::SuperOperator) = one(SuperOperator)
Base.iszero(S::SuperOperator) = isempty(S.arguments)
function Base.isone(S::SuperOperator)
    length(S.arguments) == 1 || return false
    t, c = first(S.arguments)
    return isempty(t.left) && isempty(t.right) && isequal(c, CNUM_ONE)
end

Base.length(S::SuperOperator) = length(S.arguments)
Base.iterate(S::SuperOperator) = iterate(S.arguments)
Base.iterate(S::SuperOperator, state) = iterate(S.arguments, state)
Base.eltype(::Type{SuperOperator}) = Pair{SuperTerm, CNum}
Base.keys(S::SuperOperator) = keys(S.arguments)
Base.haskey(S::SuperOperator, t::SuperTerm) = haskey(S.arguments, t)
Base.get(S::SuperOperator, t::SuperTerm, default) = get(S.arguments, t, default)

function Base.isequal(a::SuperOperator, b::SuperOperator)
    isequal(a.arguments, b.arguments) || return false
    return same_index_set(a.indices, b.indices)
end
Base.:(==)(a::SuperOperator, b::SuperOperator) = isequal(a, b)
function Base.hash(S::SuperOperator, h::UInt)
    ih = foldl((acc, idx) -> acc + hash(idx), S.indices; init = zero(UInt))
    return hash(:SuperOperator, hash(S.arguments, hash(ih, hash(length(S.indices), h))))
end

"""
    left_operator(t::SuperTerm)
    right_operator(t::SuperTerm)

The operator acting from the left or from the right in the sandwich `t`.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: sprepost, left_operator, right_operator

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> t, c = first(sprepost(a, a')); left_operator(t), right_operator(t)
(a, a')
```

See also [`SuperTerm`](@ref).
"""
left_operator(t::SuperTerm) = side_qadd(t.left, t.ne)
right_operator(t::SuperTerm) = side_qadd(t.right, t.ne)
@doc (@doc left_operator) right_operator

side_qadd(ops::Vector{Op}, ne::Vector{NonEqualPair}) =
    QAdd(QTermDict(QTerm(copy(ops), copy_ne(ne)) => CNUM_ONE), Index[])

get_sum_indices(S::SuperOperator) = copy(S.indices)

constraint_pairs(S::SuperOperator) = constraint_pairs(S.arguments)

function get_prefactor(S::SuperOperator)
    length(S.arguments) == 1 || throw(
        ArgumentError("get_prefactor requires a single-term map, got $(length(S.arguments)) terms")
    )
    return to_num(first(values(S.arguments)))
end

"""
    coefficient(S::SuperOperator, sandwich::SuperOperator) -> Coeff

The coefficient of the single sandwich `sandwich`, such as `sprepost(a, a')` or `spre(a' * a)`,
in the map `S`, or zero when `S` has no such sandwich. The rules for operator products apply
to both sides: the sandwich must be one canonical term with unit coefficient, index
constraints must match, and summed maps are refused.
"""
function coefficient(S::SuperOperator, sandwich::SuperOperator)
    isempty(S.indices) || summed_coefficient_error()
    isempty(sandwich.indices) || summed_coefficient_error()
    length(sandwich.arguments) == 1 || throw(
        ArgumentError(
            "`coefficient` needs a single sandwich with unit coefficient, but the map has " *
                "$(length(sandwich.arguments)) canonical terms"
        )
    )
    key, kc = first(sandwich.arguments)
    isone(kc) || throw(
        ArgumentError(
            "`coefficient` needs a single sandwich with unit coefficient, but the map has " *
                "coefficient $(to_num(kc))"
        )
    )
    found = CNUM_ZERO
    matched = false
    for (t, c) in S.arguments
        (isequal(t.left, key.left) && isequal(t.right, key.right)) || continue
        (matched || !isequal(t.ne, key.ne)) && ambiguous_coefficient_error()
        found = c
        matched = true
    end
    return found
end

@noinline function shared_binder_error(idx::Index)
    throw(
        ArgumentError(
            "summation index $(index_name(idx)) is bound by both operands. Bound indices " *
                "of composed maps, or of a map and the operator it acts on, must be " *
                "distinct; rename one side with `change_index`.",
        ),
    )
end

@noinline function capture_error(idx::Index)
    throw(
        ArgumentError(
            "summation index $(index_name(idx)) is bound by one operand and free in the " *
                "other, so the sum would capture the free index. Rename the bound index " *
                "with `change_index`.",
        ),
    )
end

@noinline function summed_side_error(what::AbstractString)
    throw(
        ArgumentError(
            "$what cannot take an operator with its own bound sum: `J ρ J'` with " *
                "`J = Σᵢ Jᵢ` needs two independent copies of the summation index. Use " *
                "`dissipator(Jᵢ, Γ, i => j)` for collective decay or `Σ(dissipator(Jᵢ), i)` " *
                "for independent decay, and `Σ(sprepost(Aᵢ, Bᵢ), i)` for a summed sandwich.",
        ),
    )
end

function mentions_index(S::SuperOperator, idx::Index)
    for (t, c) in S.arguments
        depends_on_scope_sandwich(c, t.left, t.right, t.ne, idx) && return true
    end
    return false
end

function check_scopes(a_indices::Vector{Index}, b_indices::Vector{Index})
    (isempty(a_indices) || isempty(b_indices)) && return nothing
    for idx in a_indices
        idx in b_indices && shared_binder_error(idx)
    end
    return nothing
end

function check_capture(binders::Vector{Index}, operand::Union{QAdd, SuperOperator})
    for idx in binders
        idx in operand.indices && continue
        mentions_index(operand, idx) && capture_error(idx)
    end
    return nothing
end

function check_operands(a::Union{QAdd, SuperOperator}, b::Union{QAdd, SuperOperator})
    check_scopes(a.indices, b.indices)
    check_capture(a.indices, b)
    check_capture(b.indices, a)
    return nothing
end

as_side(x::QAdd) = x
as_side(x::QSym) = to_qadd(x)
as_side(x::Coefficient) = single_qadd(to_cnum(x), Op[])

"""
    sprepost(A, B)

The map that multiplies an operator by `A` from the left and by `B` from the right.

It acts as `ρ ↦ A * ρ * B`. `A` and `B` are operators or scalars. A side with a formal operator function gives a
[`SuperExpr`](@ref). A side that carries its own `Σ` throws; sum the map instead,
`Σ(sprepost(Aᵢ, Bᵢ), i)`.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: sprepost

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> S = sprepost(a, a')
a ⋅ (·) ⋅ a'

julia> S(a' * a)
1 + 3 * a' * a + a' * a' * a * a
```

See also [`spre`](@ref), [`spost`](@ref), [`SuperOperator`](@ref).
"""
function sprepost(A::Union{QField, Coefficient}, B::Union{QField, Coefficient})
    return sandwich_of(as_part(A), as_part(B))
end

function sprepost_poly(A::QAdd, B::QAdd)
    isempty(A.indices) || summed_side_error("`sprepost`")
    isempty(B.indices) || summed_side_error("`sprepost`")
    out = SuperTermDict()
    for (ta, ca) in A.arguments, (tb, cb) in B.arguments
        canonical_sandwich!(out, ta.ops, tb.ops, ca * cb, merge_ne(ta.ne, tb.ne))
    end
    return superop(out, Index[])
end

"""
    spre(A)

The map that multiplies an operator by `A` from the left, `ρ ↦ A * ρ`.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: spre

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> spre(a)(a')
1 + a' * a
```

See also [`spost`](@ref), [`sprepost`](@ref).
"""
spre(A::Union{QField, Coefficient}) = sprepost(A, 1)

"""
    spost(B)

The map that multiplies an operator by `B` from the right, `ρ ↦ ρ * B`.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: spost

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> spost(a)(a')
a' * a
```

See also [`spre`](@ref), [`sprepost`](@ref).
"""
spost(B::Union{QField, Coefficient}) = sprepost(1, B)

function Base.:+(a::SuperOperator, b::SuperOperator)
    d = copy_args(a.arguments)
    for (t, c) in b.arguments
        addto_key!(d, copy_key(t), c)
    end
    return superop_pruned(d, merge_unique(copy(a.indices), b.indices))
end
Base.:+(a::SuperOperator) = a

function Base.:-(a::SuperOperator)
    d = SuperTermDict()
    sizehint!(d, length(a.arguments))
    for (t, c) in a.arguments
        d[copy_key(t)] = -c
    end
    return superop(d, copy(a.indices))
end
Base.:-(a::SuperOperator, b::SuperOperator) = a + (-b)

function scale(S::SuperOperator, c::CNum)
    iszero(c) && return zero(SuperOperator)
    isequal(c, CNUM_ONE) && return S
    d = SuperTermDict()
    sizehint!(d, length(S.arguments))
    for (t, v) in S.arguments
        nv = v * c
        iszero(nv) || addto_key!(d, copy_key(t), nv)
    end
    return superop_pruned(d, copy(S.indices))
end

Base.:*(S::SuperOperator, c::Coefficient) = scale(S, to_cnum(c))
Base.:*(c::Coefficient, S::SuperOperator) = scale(S, to_cnum(c))
Base.:/(S::SuperOperator, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? S * (1 // b) : S * inv(b)
Base.:/(S::SuperOperator, b::Coefficient) = S * inv(b)
Base.://(S::SuperOperator, b::Integer) = S * (1 // b)
Base.://(S::SuperOperator, b::Coefficient) = S / b

Base.:+(S::SuperOperator, c::Coefficient) = S + scale(one(SuperOperator), to_cnum(c))
Base.:+(c::Coefficient, S::SuperOperator) = S + c
Base.:-(S::SuperOperator, c::Coefficient) = S + (-to_cnum(c))
Base.:-(c::Coefficient, S::SuperOperator) = (-S) + c

"""
    compose(A, B)
    A ∘ B
    A * B

The composition of two maps, applying `B` first.

Composing maps that bind the same summation index throws; rename one with
[`change_index`](@ref).

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: spre, spost, compose

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> compose(spre(a), spost(a'))
a ⋅ (·) ⋅ a'
```

See also [`SuperOperator`](@ref).
"""
function compose(A::SuperOperator, B::SuperOperator)
    check_operands(A, B)
    out = SuperTermDict()
    for (ta, ca) in A.arguments
        ba = bound_indices(ta, ca, A.indices)
        for (tb, cb) in B.arguments
            bb = bound_indices(tb, cb, B.indices)
            emit_sandwich!(
                out, vcat(ta.left, tb.left), vcat(tb.right, ta.right), ca * cb,
                merge_ne(ta.ne, tb.ne), merge_unique(ba, bb),
            )
        end
    end
    return superop_pruned(out, merge_unique(copy(A.indices), B.indices))
end

Base.:*(A::SuperOperator, B::SuperOperator) = compose(A, B)
Base.:∘(A::SuperOperator, B::SuperOperator) = compose(A, B)

function apply_map(S::SuperOperator, ρ::QAdd)
    check_operands(S, ρ)
    out = QTermDict()
    for (t, c) in S.arguments
        bt = bound_indices(t, c, S.indices)
        for (r, cr) in ρ.arguments
            br = bound_indices(r, cr, ρ.indices)
            emit_applied!(
                out, vcat(t.left, r.ops, t.right), c * cr, merge_ne(t.ne, r.ne),
                merge_unique(bt, br),
            )
        end
    end
    return QAdd(out, drop_unused_indices(out, merge_unique(copy(S.indices), ρ.indices)))
end

(S::SuperOperator)(ρ::QAdd) = apply_map(S, ρ)
(S::SuperOperator)(ρ::QSym) = apply_map(S, to_qadd(ρ))
(S::SuperOperator)(ρ::Coefficient) = apply_map(S, single_qadd(to_cnum(ρ), Op[]))

@noinline function negative_power_error()
    throw(ArgumentError("negative powers of a superoperator are not supported"))
end

@noinline function scoped_power_error()
    throw(
        ArgumentError(
            "powers of a map with bound summation indices need independent copies of " *
                "those indices; compose explicitly with a renamed copy, " *
                "`S ∘ change_index(S, i, k)`",
        ),
    )
end

Base.inv(::SuperOperator) = negative_power_error()

function Base.:^(S::SuperOperator, n::Integer)
    n >= 0 || negative_power_error()
    n == 0 && return one(SuperOperator)
    n >= 2 && !isempty(S.indices) && scoped_power_error()
    result = S
    for _ in 2:n
        result = compose(result, S)
    end
    return result
end

commutator(A::SuperOperator, B::SuperOperator) = compose(A, B) - compose(B, A)
anticommutator(A::SuperOperator, B::SuperOperator) = compose(A, B) + compose(B, A)

function acts_on(S::SuperOperator)
    aon = Int[]
    for t in keys(S.arguments)
        for o in t.left
            push!(aon, o.space_index)
        end
        for o in t.right
            push!(aon, o.space_index)
        end
    end
    unique!(sort!(aon))
    return aon
end

"""
    ⊗(A::SuperOperator, B::SuperOperator)

The map that acts as `A` and `B` on two disjoint subsystems.

It equals `A ∘ B`, and throws when the two maps act on a common subspace.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: SuperOperator

julia> h = FockSpace(:a) ⊗ FockSpace(:b); a = Destroy(h, :a, 1); b = Destroy(h, :b, 2);

julia> SuperOperator(dissipator(a)) ⊗ SuperOperator(dissipator(b)) ==
           SuperOperator(dissipator(a)) ∘ SuperOperator(dissipator(b))
true
```

See also [`compose`](@ref).
"""
function ⊗(A::SuperOperator, B::SuperOperator)
    isdisjoint(acts_on(A), acts_on(B)) || throw(
        ArgumentError("`⊗` of superoperators requires maps on disjoint subspaces; use `∘`"),
    )
    return compose(A, B)
end

"""
    adjoint(S::SuperOperator)

The adjoint of a map with respect to the Hilbert–Schmidt inner product.

It satisfies `Tr(X' * S(ρ)) == Tr(S'(X)' * ρ)`, and turns each term `c * A ρ B` into
`conj(c) * A' X B'`. It is the Heisenberg-picture map only for maps that preserve
Hermiticity; `transpose` is the Heisenberg-picture map in general.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: spre

julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> adjoint(spre(a))(a' * a)
a' * a' * a
```

See also `transpose`.
"""
function Base.adjoint(S::SuperOperator)
    out = SuperTermDict()
    for (t, c) in S.arguments
        emit_sandwich!(
            out, reversed_adjoint(t.left), reversed_adjoint(t.right), conj(c), t.ne,
            bound_indices(t, c, S.indices),
        )
    end
    return superop_pruned(out, copy(S.indices))
end

"""
    transpose(S::SuperOperator)

The Heisenberg-picture map of a map on operators.

It satisfies `Tr(A * S(ρ)) == Tr(transpose(S)(A) * ρ)`, and turns each term `c * L ρ R`
into `c * R A L`. For a generator `L`, `transpose(L)(A)` is the drift of the observable
`A`, and `average(transpose(L)(A))` its equation of motion.

# Examples
```jldoctest
julia> h = FockSpace(:c); a = Destroy(h, :a);

julia> transpose(2 * dissipator(a))(a' * a)
-2 * a' * a
```

See also [`Lindbladian`](@ref).
"""
function Base.transpose(S::SuperOperator)
    out = SuperTermDict()
    sizehint!(out, length(S.arguments))
    for (t, c) in S.arguments
        out[SuperTerm(copy(t.right), copy(t.left), copy_ne(t.ne))] = c
    end
    return superop(out, copy(S.indices))
end

dagger(S::SuperOperator) = adjoint(S)
qadjoint(S::SuperOperator) = adjoint(S)
