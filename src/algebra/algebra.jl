# Scalar values accepted as coefficients in operator expressions. `Coeff` is included for
# internal composition, while public inputs ordinarily arrive as numbers or symbolic values.
const Coefficient = Union{Number, SymbolicUtils.BasicSymbolic, Coeff}

Base.:*(a::QSym, b::Coefficient) = single_qadd(to_cnum(b), Op[a])
Base.:*(b::Coefficient, a::QSym) = a * b

function Base.:*(a::QAdd, b::Coefficient)
    b isa Number && isone(b) && return a
    cb = to_cnum(b)
    d = QTermDict()
    for (term, c) in a.arguments
        new_c = mul_cnum(c, cb)
        iszero_cnum(new_c) && continue
        d[copy_key(term)] = new_c
    end
    return QAdd(d, copy(a.indices))
end
Base.:*(a::Coefficient, b::QAdd) = b * a

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    addto!(d, Op[a], CNUM_ONE)
    addto!(d, Op[b], CNUM_ONE)
    return QAdd(d, EMPTY_INDICES)
end

function Base.:+(a::QAdd, b::QSym)
    d = copy_args(a.arguments)
    addto!(d, Op[b], CNUM_ONE)
    return QAdd(d, drop_unused_indices(d, a.indices))
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = copy_args(a.arguments)
    for (term, c) in b.arguments
        addto_key!(d, copy_key(term), c)
    end
    return QAdd(d, drop_unused_indices(d, merge_unique(a.indices, b.indices)))
end

function Base.:+(a::QSym, b::Coefficient)
    d = QTermDict()
    addto!(d, Op[a], CNUM_ONE)
    addto!(d, EMPTY_OPS, to_cnum(b))
    return QAdd(d, EMPTY_INDICES)
end
Base.:+(a::Coefficient, b::QSym) = b + a

function Base.:+(a::QAdd, b::Coefficient)
    # `iszero` on a `BasicSymbolic` is itself symbolic, so guard the shortcut on `Number`.
    b isa Number && iszero(b) && return a
    d = copy_args(a.arguments)
    addto!(d, EMPTY_OPS, to_cnum(b))
    return QAdd(d, copy(a.indices))
end
Base.:+(a::Coefficient, b::QAdd) = b + a

Base.zero(::Type{QAdd}) = zero_qadd()
Base.zero(::QAdd) = zero_qadd()

Base.:+(a::QSym) = single_qadd(CNUM_ONE, Op[a])
Base.:+(a::QAdd) = a

Base.:-(a::QSym) = single_qadd(CNUM_NEG1, Op[a])

function Base.:-(a::QAdd)
    d = QTermDict()
    for (term, c) in a.arguments
        d[copy_key(term)] = neg_cnum(c)
    end
    return QAdd(d, copy(a.indices))
end

Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Coefficient) = a + (-b)
Base.:-(a::Coefficient, b::QField) = a + (-b)

Base.:/(a::QSym, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? a * (1 // b) : a * inv(b)
Base.:/(a::QAdd, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? a * (1 // b) : a * inv(b)
Base.:/(a::QSym, b::Coefficient) = a * inv(b)
Base.:/(a::QAdd, b::Coefficient) = a * inv(b)

Base.://(a::QSym, b::Integer) = a * (1 // b)
Base.://(a::QAdd, b::Integer) = a * (1 // b)
# A symbolic or non-integer denominator cannot be the denominator of Julia's
# rational literal, but operator division supports all scalar coefficients.
Base.://(a::QSym, b::Coefficient) = a / b
Base.://(a::QAdd, b::Coefficient) = a / b

function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return single_qadd(CNUM_ONE, EMPTY_OPS)
    out = QTermDict()
    ops = Op[a for _ in 1:n]
    canonicalize!(out, ops, CNUM_ONE, EMPTY_NE)
    return QAdd(out, Index[])
end

function Base.:^(a::QAdd, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return single_qadd(CNUM_ONE, EMPTY_OPS)
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

"""
    normal_order(expr::QField) -> QAdd

Route every term of `expr` through the canonicalization pipeline.

In practice this is the identity on anything built through public arithmetic:
`*`, `+`, `-`, `^`, [`commutator`](@ref), [`Σ`](@ref), [`substitute`](@ref),
and `adjoint` all canonicalize eagerly, so the result of any such call is
already normal-ordered. Reach for `normal_order` explicitly only when an
expression was assembled through low-level internals that bypass the
arithmetic, or when interfacing with code that expects a finalizer call.
[`simplify`](@ref) uses it internally before simplifying coefficients.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> a * a'
1 + a' * a

julia> normal_order(a * a')
1 + a' * a
```

See also [`simplify`](@ref), [`expand_completeness`](@ref),
[`normal_to_symmetric`](@ref), [`symmetric_to_normal`](@ref).
"""
normal_order(op::QSym) = single_qadd(CNUM_ONE, Op[op])

function normal_order(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        stream!(out, copy(t.ops), c, t.ne)
    end
    return QAdd(out, copy(q.indices))
end

function simplify_raw_component(x; kwargs...)
    return Num(SymbolicUtils.simplify(SymbolicUtils.expand(x); kwargs...))
end

function contains_division(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(x) || return false
    SymbolicUtils.operation(x) === (/) && return true
    for argument in SymbolicUtils.arguments(x)
        contains_division(argument) && return true
    end
    return false
end

function simplify_prefactor(x::CNum; kwargs...)
    is_native(x) && return x                    # already simplest product form
    is_poly(x) && return reduce_trig(x)        # the CAS is not reached on this tier
    tail = x.tail
    if !contains_phase(tail.expr) && !contains_division(tail.expr)
        re, im = raw_realimag(tail.expr)
        return cnum(
            simplify_raw_component(re; kwargs...), simplify_raw_component(im; kwargs...),
        )
    end
    expanded = SymbolicUtils.expand(tail.expr)
    return to_cnum(simplify_raw(expanded; kwargs...))
end

SymbolicUtils.simplify(c::Coeff; kwargs...) = simplify_prefactor(c; kwargs...)
function drop_unused_indices(d::QTermDict, indices::Vector{Index})
    isempty(indices) && return indices
    used = Index[]
    for idx in indices
        for (term, c) in d
            if depends_on_index_term(c, term.ops, idx)
                push!(used, idx)
                break
            end
        end
    end
    length(used) == length(indices) && return indices
    return used
end

"""
    simplify(expr::QField) -> QAdd

Normal-order `expr`, then simplify each coefficient symbolically and drop
summation indices that no surviving term depends on.

The operator-level work (commutation, same-site composition, like-term
collection) all happens inside [`normal_order`](@ref). What `simplify` adds
on top is purely at the symbolic-coefficient layer: `Symbolics.expand`
followed by `SymbolicUtils.simplify` runs on each surviving prefactor, and
any term whose coefficient simplifies to zero is dropped. A final pass
removes summation indices that no remaining term references.

That symbolic step is expensive, so reach for `simplify` as a finalizer
when cancellations or accumulated symbolic factors need to be folded; use
`normal_order` for intermediate steps.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x y;

julia> expr = (x^2 + 2x*y + y^2) * a' * a - (x + y)^2 * a' * a;

julia> simplify(expr)
0
```

See also [`normal_order`](@ref), [`expand`](@ref), [`expand_completeness`](@ref).
"""
SymbolicUtils.simplify(op::QSym; kwargs...) = single_qadd(CNUM_ONE, Op[op])

SymbolicUtils.simplify(q::QAdd; kwargs...) =
    map_coefficients(simplify_prefactor, normal_order(q))

# Restoring the `QAdd` contract (no zero entries, every advertised index has a live carrier)
# after a coefficient rewrite, once here rather than at each caller.
function map_coefficients(f::F, q::QAdd) where {F}
    out = QTermDict()
    for (term, c) in q.arguments
        new_c = f(c)
        iszero_cnum(new_c) && continue
        addto_key!(out, copy_key(term), new_c)
    end
    return QAdd(out, drop_unused_indices(out, q.indices))
end

exponential_form(op::QSym) =
    exponential_form(single_qadd(CNUM_ONE, Op[op]))
exponential_form(q::QAdd) = map_coefficients(exponential_form, q)

trigonometric_form(op::QSym) =
    trigonometric_form(single_qadd(CNUM_ONE, Op[op]))
trigonometric_form(q::QAdd) = map_coefficients(trigonometric_form, q)

"""
    expand(expr::QField) -> QAdd

Expand the symbolic prefactor of each term via `Symbolics.expand`.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x y;

julia> expand((x + y)^2 * a)
(x^2 + 2x*y + y^2) * a
```

See also [`simplify`](@ref).
"""
function Symbolics.expand(s::QAdd; kwargs...)
    d = QTermDict()
    for (term, c) in s.arguments
        new_c = expand_prefactor(c; kwargs...)
        iszero_cnum(new_c) && continue
        addto!(d, term.ops, new_c, term.ne)
    end
    return QAdd(d, copy(s.indices))
end
Symbolics.expand(op::QSym; kwargs...) = single_qadd(CNUM_ONE, Op[op])

expand_prefactor(x::CNum; kwargs...) = (is_native(x) || is_poly(x) || iszero_cnum(x)) ? x : cnum(Symbolics.expand(real(x)), Symbolics.expand(imag(x)))

"""
    expand_completeness(q) -> QAdd

Rewrite every ground-state projector ``\\sigma^{gg}`` in `q` via the
completeness relation ``\\sigma^{gg} = 1 - \\sum_{k \\neq g} \\sigma^{kk}``.

`*` keeps ``\\sigma^{gg}`` atomic; reach for `expand_completeness` when
downstream code needs the projector eliminated, e.g. before converting to a
numeric basis where the ``\\sigma^{kk}`` for ``k \\neq g`` form the
independent degrees of freedom.

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, 2);

julia> σ11 = Transition(h, :σ, 1, 1);

julia> expand_completeness(σ11)
1 - σ₂₂
```

See also [`assume_distinct_index`](@ref), [`normal_order`](@ref).
"""
function expand_completeness(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        depended = Index[idx for idx in q.indices if depends_on_index_term(c, t.ops, idx)]
        if isempty(depended)
            expand_gs_ops(copy(t.ops), c) do ops1, c1
                stream!(out, ops1, c1, t.ne)
            end
        else
            expand_gs_ops(copy(t.ops), c) do ops1, c1
                emit_scaled_by_scope!(out, ops1, c1, t.ne, depended)
            end
        end
    end
    return QAdd(out, copy(q.indices))
end

expand_completeness(op::QSym) = expand_completeness(single_qadd(CNUM_ONE, Op[op]))

"""
    assume_distinct_index(q::QAdd, pairs::Vector{Tuple{Index, Index}}) -> QAdd

Re-canonicalize `q` under the declared pairwise `≠` constraints on free
indices, then run [`expand_completeness`](@ref) so any ground-state
projectors that emerge from same-site composition are expanded.

Use this when two free indices semantically range over distinct atoms or
modes but no `Σ` supplies the constraint. SQA cannot infer "different symbol
implies different site" automatically: two operators carrying different free
indices on the same Hilbert subspace remain in their physical order, and no
same-site collapse fires between them. Declaring the pair distinct here lets
the canonical sort resolve their relationship and triggers any composition or
completeness rewriting it unlocks.

# Examples

```jldoctest
julia> h = NLevelSpace(:atom, 2);

julia> j = Index(h, :j, 5, h); k = Index(h, :k, 5, h);

julia> σ(i, m, idx) = IndexedOperator(Transition(h, :σ, i, m), idx);

julia> q = σ(2, 1, k) * σ(1, 2, j);

julia> assume_distinct_index(q, [(j, k)])
σ_j₁₂ * σ_k₂₁
```

See also [`expand_completeness`](@ref).
"""
function assume_distinct_index(q::QAdd, pairs::Vector{Tuple{Index, Index}})
    out = QTermDict()
    for (t, c) in q
        ne = t.ne
        for (a, b) in pairs
            ne = merge_ne_pair(ne, a, b)
        end
        canonicalize!(out, copy(t.ops), c, ne)
    end
    return expand_completeness(QAdd(out, copy(q.indices)))
end

const ZERO_QADD = QAdd(QTermDict(), Index[])
zero_qadd() = ZERO_QADD

"""
    commutator(a, b) -> QAdd

Return the commutator ``[a, b] = a\\,b - b\\,a`` as a [`QAdd`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> commutator(a, a')
1
```

See also [`anticommutator`](@ref), [`normal_order`](@ref).
"""
function commutator end

commutator(::Number, ::Number) = zero_qadd()
commutator(::Number, ::QField) = zero_qadd()
commutator(::QField, ::Number) = zero_qadd()

function commutator(a::QSym, b::QSym)
    # `can_commute`/`commute_pair` below are only defined on provably-same-site
    # pairs, which is exactly `Equal`. Anything else commutes or is undetermined.
    site_compare(a, b, EMPTY_NE) === Equal || return zero_qadd()
    isequal(a, b) && return zero_qadd()
    # If exactly one direction needs a swap, [a, b] is the swap residual.
    forward = can_commute(a, b)
    reverse = can_commute(b, a)
    if !forward && reverse
        _, _, c1, ops1, c2, ops2 = commute_pair(a, b)
        return single_qadd(c1, isempty(ops1) ? EMPTY_OPS : copy(ops1)) +
            single_qadd(c2, isempty(ops2) ? EMPTY_OPS : copy(ops2))
    elseif forward && !reverse
        _, _, c1, ops1, c2, ops2 = commute_pair(b, a)
        return single_qadd(neg_cnum(c1), isempty(ops1) ? EMPTY_OPS : copy(ops1)) +
            single_qadd(neg_cnum(c2), isempty(ops2) ? EMPTY_OPS : copy(ops2))
    end
    return a * b - b * a
end

# Skip term pairs on disjoint subspaces: [aₖ,bₗ]=0, so their products only cancel.
function commutator(a::QAdd, b::QSym)
    sb = b.space_index
    d = QTermDict()
    if isempty(a.indices)
        for (ta, ca) in a.arguments
            sb in acts_on(ta) || continue
            emit_product!(
                d, ta.ops, ca, ta.ne, Op[b], CNUM_ONE, EMPTY_NE, EMPTY_INDICES, false,
            )
            emit_product!(
                d, Op[b], CNUM_NEG1, EMPTY_NE, ta.ops, ca, ta.ne, EMPTY_INDICES, false,
            )
        end
        return QAdd(d, EMPTY_INDICES)
    end
    indices = EMPTY_INDICES
    for (ta, ca) in a.arguments
        sb in acts_on(ta) || continue
        qa = QAdd(QTermDict(copy_key(ta) => ca), a.indices)
        indices = merge_unique(indices, merge_into!(d, qa * b - b * qa))
    end
    return QAdd(d, drop_unused_indices(d, indices))
end
commutator(a::QSym, b::QAdd) = -commutator(b, a)

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return zero_qadd()
    bside = [(tb, cb, acts_on(tb)) for (tb, cb) in b.arguments]
    d = QTermDict()
    if isempty(a.indices) && isempty(b.indices)
        for (ta, ca) in a.arguments
            aon_a = acts_on(ta)
            for (tb, cb, aon_b) in bside
                isdisjoint(aon_a, aon_b) && continue
                emit_product!(
                    d, ta.ops, ca, ta.ne, tb.ops, cb, tb.ne, EMPTY_INDICES, false,
                )
                emit_product!(
                    d, tb.ops, neg_cnum(cb), tb.ne, ta.ops, ca, ta.ne, EMPTY_INDICES, false,
                )
            end
        end
        return QAdd(d, EMPTY_INDICES)
    end
    indices = EMPTY_INDICES
    for (ta, ca) in a.arguments
        aon_a = acts_on(ta)
        qa = QAdd(QTermDict(copy_key(ta) => ca), a.indices)
        for (tb, cb, aon_b) in bside
            isdisjoint(aon_a, aon_b) && continue
            qb = QAdd(QTermDict(copy_key(tb) => cb), b.indices)
            indices = merge_unique(indices, merge_into!(d, qa * qb - qb * qa))
        end
    end
    return QAdd(d, drop_unused_indices(d, indices))
end

# Collect every term of `r` into `d` (summing like terms); return `r`'s indices
# so the caller can accumulate the summation scope as the `+` chain would.
function merge_into!(d::QTermDict, r::QAdd)
    for (term, c) in r.arguments
        addto_key!(d, copy_key(term), c)
    end
    return r.indices
end

"""
    anticommutator(a, b) -> QAdd

Return the anticommutator ``\\{a, b\\} = a\\,b + b\\,a`` as a [`QAdd`](@ref).

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

"""
    substitute(expr::QField, rules::AbstractDict; replace_adjoint=true)

Substitute symbolic parameters and/or operators in `expr` using `rules`.

Scalar keys (anything other than a [`QSym`](@ref)) are substituted in
coefficients using SymbolicUtils. Operator keys are applied simultaneously to
the original operator leaves: replacement expressions are not searched again for
operator keys. This makes mode transformations such as `a => g*a + h*b`
well-defined even when the replacement contains the original operator.

When `replace_adjoint=true`, missing adjoint operator rules are added
automatically, e.g. `a => b` also implies `a' => b'` unless `a'` is already
specified explicitly.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x;

julia> substitute(x * a' * a, Dict(x => 2))
2 * a' * a
```

See also [`change_index`](@ref).
"""
function SymbolicUtils.substitute(op::QSym, rules::AbstractDict; replace_adjoint = true)
    return SymbolicUtils.substitute(
        single_qadd(CNUM_ONE, Op[op]), rules; replace_adjoint = replace_adjoint
    )
end

function SymbolicUtils.substitute(q::QAdd, rules::AbstractDict; replace_adjoint = true)
    iszero(q) && return q
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    return substitute_split(q, op_rules, scalar_rules)
end

# Never mutated: `substitute_cnum` early-returns on an empty dict and nothing else writes.
const EMPTY_SCALAR_RULES = Dict{Any, Any}()

# For a rule set already known to be all-operator. The public `substitute` has to copy an
# `AbstractDict` into two fresh `Dict{Any, Any}`, boxing every key and value, to find that out.
substitute_op_rules(q::QAdd, op_rules::Dict{Op, QAdd}) =
    iszero(q) ? q : substitute_split(q, op_rules, EMPTY_SCALAR_RULES)

function substitute_split(q::QAdd, op_rules::AbstractDict, scalar_rules::AbstractDict)
    out = QTermDict()
    indices = copy(q.indices)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    for (t, c) in q
        replacement_indices = stream_substitution_once!(
            out, t.ops, c, t.ne, op_rules, scalar_rules, phase_rules,
        )
        indices = merge_unique(indices, replacement_indices)
    end
    return QAdd(out, drop_unused_indices(out, indices))
end

function split_substitution_rules(rules::AbstractDict, replace_adjoint)
    op_rules = Dict{Any, Any}()
    scalar_rules = Dict{Any, Any}()
    for (k, v) in rules
        if k isa QSym
            op_rules[k] = v
        else
            scalar_rules[k] = v
        end
    end
    replace_adjoint && add_adjoint_rules!(op_rules)
    return op_rules, scalar_rules
end

function add_adjoint_rules!(op_rules::Dict{Any, Any})
    for (k, v) in collect(op_rules)
        k_adj = Base.adjoint(k)
        haskey(op_rules, k_adj) || (op_rules[k_adj] = adjoint_replacement(v))
    end
    return op_rules
end

function has_operator_rule(ops::Vector{Op}, dict::AbstractDict)
    for op in ops
        haskey(dict, op) && return true
    end
    return false
end

function stream_substitution_once!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
        op_rules::AbstractDict, scalar_rules::AbstractDict, phase_rules::Vector{Pair{Any, Any}},
    )
    if !has_operator_rule(ops, op_rules)
        new_c = substitute_cnum(c, scalar_rules, phase_rules)
        iszero_cnum(new_c) || stream!(out, copy(ops), new_c, ne)
        return EMPTY_INDICES
    end

    partials = Tuple{Vector{Op}, CNum, Vector{NonEqualPair}}[(Op[], c, EMPTY_NE)]
    indices = EMPTY_INDICES
    for op in ops
        replacement = replacement_qadd(op, op_rules)
        indices = merge_unique(indices, replacement.indices)
        next_partials = Tuple{Vector{Op}, CNum, Vector{NonEqualPair}}[]
        for (prefix_ops, prefix_c, prefix_ne) in partials
            for (rt, rc) in replacement
                new_ops = copy(prefix_ops)
                append!(new_ops, rt.ops)
                push!(
                    next_partials,
                    (new_ops, mul_cnum(prefix_c, rc), merge_ne(prefix_ne, rt.ne)),
                )
            end
        end
        partials = next_partials
    end

    for (new_ops, new_c0, new_ne) in partials
        new_c = substitute_cnum(new_c0, scalar_rules, phase_rules)
        iszero_cnum(new_c) && continue
        canonicalize!(out, new_ops, new_c, merge_ne(ne, new_ne))
    end
    return indices
end

function replacement_qadd(op::Op, dict::AbstractDict)
    haskey(dict, op) || return single_qadd(CNUM_ONE, Op[op])
    val = dict[op]
    val isa QAdd && return val
    val isa QSym && return single_qadd(CNUM_ONE, Op[val])
    (val isa Number || val isa SymbolicUtils.BasicSymbolic || val isa Coeff) &&
        return single_qadd(to_cnum(val), Op[])
    throw(ArgumentError("operator replacement for `$op` has unsupported value `$val`"))
end

adjoint_replacement(v::Coeff) = conj(v)
adjoint_replacement(v) = qadjoint(v)
