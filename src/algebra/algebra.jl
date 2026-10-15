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
    normal_order(expr::QField)

Bring `expr` into canonical, normal-ordered form.

Arithmetic already returns this form: `*`, `+`, `-`, `^`, `adjoint`,
[`commutator`](@ref), [`Σ`](@ref) and [`substitute`](@ref) all canonicalize as they
go, so `normal_order` leaves their results unchanged. It is there for code that
expects an explicit ordering step. Formal operator functions such as `cos` or
[`expim`](@ref) keep their form; only their polynomial arguments are ordered.

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
    is_native(x) && return x
    is_poly(x) && return reduce_trig(x)
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
            if depends_on_scope_index(c, term.ops, term.ne, idx)
                push!(used, idx)
                break
            end
        end
    end
    length(used) == length(indices) && return indices
    return used
end

"""
    simplify(expr::QField)

Simplify the scalar coefficients of `expr` symbolically.

Polynomial coefficients already cancel as an expression is built. `simplify` goes
further and applies symbolic identities, such as trigonometric ones, to every
coefficient, then drops the terms that vanish. Because this step is expensive, it
is best used once on a final result. Formal operator functions keep their form.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x;

julia> expr = (sin(x)^2 + cos(x)^2) * a' * a
(cos(x)^2 + sin(x)^2) * a' * a

julia> simplify(expr)
a' * a
```

See also [`normal_order`](@ref), [`expand`](@ref), [`expand_completeness`](@ref).
"""
SymbolicUtils.simplify(op::QSym; kwargs...) = single_qadd(CNUM_ONE, Op[op])

SymbolicUtils.simplify(q::QAdd; kwargs...) =
    map_coefficients(simplify_prefactor, normal_order(q))

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
    expand(expr::QField)

Multiply out the scalar coefficients of `expr`.

Polynomial coefficients are kept expanded already; `expand` also reaches products
inside functions such as `exp`. Formal operator functions are not expanded.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x y;

julia> expr = exp(x * (x + y)) * a;

julia> expand(expr)
exp(x^2 + x*y) * a
```

See also [`simplify`](@ref), [`taylor`](@ref).
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
    expand_completeness(expr)

Eliminate every ground-state projector in `expr` through completeness.

Each projector ``\\sigma^{gg}`` is replaced by ``1 - \\sum_{k \\neq g} \\sigma^{kk}``.
Arithmetic never does this by itself, which keeps expressions short. Expand when
the remaining projectors should be the independent variables, e.g. before
deriving equations for their expectation values.

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
    assume_distinct_index(expr, pairs::Vector{Tuple{Index, Index}})

Treat each pair of free indices in `pairs` as different sites.

Two free indices on the same subsystem may or may not refer to the same site, so their
operators are left as they are. Outside a [`Σ`](@ref), which settles this by itself, this
function supplies the missing information: the operators are then sorted into canonical
order, and ground-state projectors in the result are eliminated with
[`expand_completeness`](@ref).

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> σ(2, 1, j) * σ(1, 2, i)
σ_j₂₁ * σ_i₁₂

julia> assume_distinct_index(σ(2, 1, j) * σ(1, 2, i), [(i, j)])
σ_i₁₂ * σ_j₂₁
```

See also [`Σ`](@ref), [`expand_completeness`](@ref).
"""
function assume_distinct_index(q::QAdd, pairs::Vector{Tuple{Index, Index}})
    out = QTermDict()
    for (t, c) in q
        ne = t.ne
        for (a, b) in pairs
            ne = merge_ne_pair(ne, a, b)
        end
        temp = QTermDict()
        canonicalize!(temp, copy(t.ops), c, ne)
        settle_scope!(out, temp, bound_indices(t, c, q.indices))
    end
    return expand_completeness(QAdd(out, drop_unused_indices(out, copy(q.indices))))
end

const ZERO_QADD = QAdd(QTermDict(), Index[])
zero_qadd() = ZERO_QADD

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

function merge_into!(d::QTermDict, r::QAdd)
    for (term, c) in r.arguments
        addto_key!(d, copy_key(term), c)
    end
    return r.indices
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

"""
    substitute(expr::QField, rules::AbstractDict; replace_adjoint=true)

Replace parameters and operators in `expr` according to `rules`.

Keys may be parameters or operators, mixed freely. Operator rules act on the
original expression all at once, so a rule whose replacement contains its own
key, such as the shift `a => a + α`, is applied exactly once. A rule for `a` also
covers `a'`, unless `a'` has a rule of its own or `replace_adjoint = false`.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables Δ α::Number;

julia> substitute(Δ * a' * a, Dict(Δ => 2))
2 * a' * a

julia> substitute(a' * a, Dict(a => a + α))
conj(α)*α + conj(α) * a + α * a' + a' * a
```

See also [`change_index`](@ref), [`Displace`](@ref).
"""
function SymbolicUtils.substitute(op::QSym, rules::AbstractDict; replace_adjoint = true)
    return SymbolicUtils.substitute(
        single_qadd(CNUM_ONE, Op[op]), rules; replace_adjoint = replace_adjoint
    )
end

function SymbolicUtils.substitute(q::QAdd, rules::AbstractDict; replace_adjoint::Bool = true)
    iszero(q) && return q
    op_rules, scalar_rules = split_substitution_rules(rules, replace_adjoint)
    return substitute_split(q, op_rules, scalar_rules)
end

const OpRules = Dict{Op, QAdd}

const EMPTY_OP_RULES = OpRules()

substitute_op_rules(q::QAdd, op_rules::OpRules) =
    iszero(q) ? q : substitute_split(q, op_rules, EMPTY_SCALAR_RULES)

function substitute_split(q::QAdd, op_rules::OpRules, scalar_rules::ScalarRules)
    out = QTermDict()
    indices = copy(q.indices)
    phase_rules = nonreal_phase_substitutions(scalar_rules)
    for (t, c) in q
        temp = QTermDict()
        replacement_indices = stream_substitution_once!(
            temp, t.ops, c, t.ne, op_rules, scalar_rules, phase_rules,
        )
        settle_scope!(out, temp, bound_indices(t, c, q.indices))
        indices = merge_unique(indices, replacement_indices)
    end
    return QAdd(out, drop_unused_indices(out, indices))
end

function split_substitution_rules(rules::AbstractDict, replace_adjoint::Bool)
    op_rules = OpRules()
    scalar_rules = ScalarRules()
    for (k, v) in rules
        add_substitution_rule!(op_rules, scalar_rules, k, v)
    end
    replace_adjoint && add_adjoint_rules!(op_rules)
    return op_rules, scalar_rules
end

function add_substitution_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k::Op, v)
    op_rules[k] = replacement_qadd(k, v)
    return op_rules
end
function add_substitution_rule!(op_rules::OpRules, scalar_rules::ScalarRules, k, v)
    add_scalar_rule!(scalar_rules, k, v)
    return op_rules
end

function add_adjoint_rules!(op_rules::OpRules)
    for (k, v) in collect(op_rules)
        k_adj = Base.adjoint(k)
        haskey(op_rules, k_adj) || (op_rules[k_adj] = qadjoint(v))
    end
    return op_rules
end

function has_operator_rule(ops::Vector{Op}, dict::OpRules)
    for op in ops
        haskey(dict, op) && return true
    end
    return false
end

function stream_substitution_once!(
        out::QTermDict, ops::Vector{Op}, c::CNum, ne::Vector{NonEqualPair},
        op_rules::OpRules, scalar_rules::ScalarRules, phase_rules::PhaseRules,
    )
    if !has_operator_rule(ops, op_rules)
        new_c = substitute_cnum(c, scalar_rules, phase_rules)
        iszero_cnum(new_c) || stream!(out, copy(ops), new_c, ne)
        return EMPTY_INDICES
    end

    partials = Tuple{Vector{Op}, CNum, Vector{NonEqualPair}}[(Op[], c, EMPTY_NE)]
    indices = EMPTY_INDICES
    for op in ops
        replacement = get(() -> single_qadd(CNUM_ONE, Op[op]), op_rules, op)
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

replacement_qadd(::Op, v::QAdd) = v
replacement_qadd(::Op, v::Op) = single_qadd(CNUM_ONE, Op[v])
replacement_qadd(::Op, v::Union{Number, SymbolicUtils.BasicSymbolic, Coeff}) = single_qadd(to_cnum(v), Op[])
replacement_qadd(op::Op, v) =
    throw(ArgumentError("operator replacement for `$op` has unsupported value `$v`"))
