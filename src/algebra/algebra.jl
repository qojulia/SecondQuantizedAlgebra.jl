Base.:*(a::QSym, b::Coefficient) = single_qadd(to_cnum(b), Op[a])
Base.:*(b::Coefficient, a::QSym) = a * b

function Base.:*(a::QAdd, b::Coefficient)
    b isa Number && isone(b) && return a
    cb = to_cnum(b)
    d = QTermDict()
    for (term, c) in a.arguments
        new_c = c * cb
        iszero(new_c) && continue
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
        d[copy_key(term)] = -c
    end
    return QAdd(d, copy(a.indices))
end

Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Coefficient) = a + (-b)
Base.:-(a::Coefficient, b::QField) = a + (-b)

Base.:/(a::Union{QSym, QAdd}, b::Number) =
    b isa Integer && !(b isa Bool) && !iszero(b) ? a * (1 // b) : a * inv(b)
Base.:/(a::Union{QSym, QAdd}, b::Coefficient) = a * inv(b)

Base.://(a::Union{QSym, QAdd}, b::Integer) = a * (1 // b)
Base.://(a::Union{QSym, QAdd}, b::Coefficient) = a / b

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
    used = Index[idx for idx in indices if any_depends_on_index(d, idx)]
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

"""
    map_coefficients(f, x)

Apply `f` to every coefficient of `x` and rebuild it in canonical form. `x` is a
[`QAdd`](@ref) (or a single operator), a [`SuperOperator`](@ref) or a polynomial
[`Lindbladian`](@ref). `f` receives a [`Coeff`](@ref) and may return a `Number`, a
Symbolics `Num` or a `Coeff`; the result is converted exactly. Terms that map to zero are
dropped and like terms merge.

`f` acts per coefficient, so the operators are kept as they are. Summation scopes are kept:
a term whose coefficient loses its dependence on a summation index is multiplied by that
index's count, as with [`substitute`](@ref), so `Σᵢ gᵢ a'a` with `gᵢ ↦ g` gives `N g a'a`.
`f` must not introduce a summation index into a term it does not bind; that throws an
`ArgumentError`.

For a `Lindbladian`, `f` maps the coefficients of the Hamiltonian, of every jump and every
rate, and each channel is renormalized. This agrees with mapping the lowered
`SuperOperator` when `f` is a substitution of real parameters or preserves values, such as
`simplify`; a general `f` need not commute with lowering, because a lowered coefficient
`γ c_J conj(c_K)` is a product of the stored ones.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> @variables x;

julia> SecondQuantizedAlgebra.map_coefficients(c -> substitute(c, Dict(x => 2)), x * a' * a + a)
a + 2 * a' * a
```

See also [`simplify`](@ref), [`substitute`](@ref).
"""
function map_coefficients(f::F, q::QAdd) where {F}
    out = QTermDict()
    for (term, c) in q.arguments
        new_c = to_cnum(f(c))
        iszero(new_c) && continue
        if isempty(q.indices) || new_c === c
            addto_key!(out, copy_key(term), new_c)
            continue
        end
        mentions(x, idx) = depends_on_index_term(x, term.ops, idx)
        coef, ne = settle_mapped(mentions, c, new_c, term.ne, q.indices)
        key = ne === term.ne ? copy_key(term) : QTerm(copy(term.ops), ne)
        addto_key!(out, key, coef)
    end
    return QAdd(out, drop_unused_indices(out, q.indices))
end
map_coefficients(f::F, op::QSym) where {F} =
    map_coefficients(f, single_qadd(CNUM_ONE, Op[op]))

"""
    settle_mapped(mentions, c, new_c, ne, indices) -> (coef, ne)

Settle the scope of one term whose coefficient changed from `c` to `new_c`: a binder of the
term that `new_c` no longer needs is counted, and a scope index that `new_c` adds to a term
it did not bind is refused.
"""
function settle_mapped(
        mentions::M, c::CNum, new_c::CNum, ne::Vector{NonEqualPair}, indices::Vector{Index},
    ) where {M}
    bound = Index[]
    for idx in indices
        if mentions(c, idx) || any(p -> p[1] == idx || p[2] == idx, ne)
            push!(bound, idx)
        elseif mentions(new_c, idx)
            mapped_capture_error(idx)
        end
    end
    return eliminate_binders(mentions, new_c, ne, bound)
end

@noinline function mapped_capture_error(idx::Index)
    throw(
        ArgumentError(
            "the coefficient map introduced summation index $(index_name(idx)) into a term " *
                "that the sum does not bind; rename the index with `change_index`.",
        ),
    )
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
        iszero(new_c) && continue
        addto!(d, term.ops, new_c, term.ne)
    end
    return QAdd(d, copy(s.indices))
end
Symbolics.expand(op::QSym; kwargs...) = single_qadd(CNUM_ONE, Op[op])

expand_prefactor(x::CNum; kwargs...) = (is_native(x) || is_poly(x) || iszero(x)) ? x : cnum(Symbolics.expand(real(x)), Symbolics.expand(imag(x)))

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
