Base.:*(a::QSym, b::Number) = _single_qadd(_to_cnum(b), QSym[a])
Base.:*(b::Number, a::QSym) = a * b

function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict()
    for (term, c) in a.arguments
        new_c = _mul_cnum(c, cb)
        _iszero_cnum(new_c) && continue
        d[_copy_key(term)] = new_c
    end
    return QAdd(d, copy(a.indices))
end
Base.:*(a::Number, b::QAdd) = b * a

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, _EMPTY_INDICES)
end

function Base.:+(a::QAdd, b::QSym)
    d = _copy_args(a.arguments)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, copy(a.indices))
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = _copy_args(a.arguments)
    for (term, c) in b.arguments
        _addto_key!(d, _copy_key(term), c)
    end
    return QAdd(d, _merge_unique(a.indices, b.indices))
end

function Base.:+(a::QSym, b::Number)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, _EMPTY_OPS, _to_cnum(b))
    return QAdd(d, _EMPTY_INDICES)
end
Base.:+(a::Number, b::QSym) = b + a

function Base.:+(a::QAdd, b::Number)
    d = _copy_args(a.arguments)
    _addto!(d, _EMPTY_OPS, _to_cnum(b))
    return QAdd(d, copy(a.indices))
end
Base.:+(a::Number, b::QAdd) = b + a

Base.zero(::Type{QAdd}) = _zero_qadd()
Base.zero(::QAdd) = _zero_qadd()

Base.:-(a::QSym) = _single_qadd(_CNUM_NEG1, QSym[a])

function Base.:-(a::QAdd)
    d = QTermDict()
    for (term, c) in a.arguments
        d[_copy_key(term)] = _neg_cnum(c)
    end
    return QAdd(d, copy(a.indices))
end

Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

Base.:/(a::QSym, b::Number) = a * inv(b)
Base.:/(a::QAdd, b::Number) = a * inv(b)

Base.://(a::QSym, b::Integer) = a * (1 // b)
Base.://(a::QAdd, b::Integer) = a * (1 // b)

function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, _EMPTY_OPS)
    out = QTermDict()
    ops = QSym[a for _ in 1:n]
    _canonicalize!(out, ops, _CNUM_ONE, _EMPTY_NE)
    return QAdd(out, Index[])
end

function Base.:^(a::QAdd, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, _EMPTY_OPS)
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
normal_order(op::QSym) = _single_qadd(_CNUM_ONE, QSym[op])

function normal_order(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        _stream!(out, copy(t.ops), c, t.ne)
    end
    return QAdd(out, copy(q.indices))
end

function _simplify_prefactor(x::CNum)
    SymbolicUtils.unwrap(x.re) isa Number && SymbolicUtils.unwrap(x.im) isa Number && return x
    re = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(Symbolics.expand(real(x)))))
    im = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(Symbolics.expand(imag(x)))))
    return Complex(re, im)
end
function _drop_unused_indices(d::QTermDict, indices::Vector{Index})
    isempty(indices) && return indices
    used = Index[]
    for idx in indices
        for (term, c) in d
            if _depends_on_index_term(c, term.ops, idx)
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

julia> expr = (x^2 + 2x*y + y^2) * a' * a - (x + y)^2 * a' * a
(x^2 - ((x + y)^2) + 2x*y + y^2) * a' * a

julia> simplify(expr)
0
```

See also [`normal_order`](@ref), [`expand`](@ref), [`expand_completeness`](@ref).
"""
SymbolicUtils.simplify(op::QSym; kwargs...) = _single_qadd(_CNUM_ONE, QSym[op])

function SymbolicUtils.simplify(q::QAdd; kwargs...)
    nq = normal_order(q)
    out = QTermDict()
    for (term, c) in nq.arguments
        new_c = _simplify_prefactor(c)
        _iszero_cnum(new_c) && continue
        _addto_key!(out, _copy_key(term), new_c)
    end
    return QAdd(out, _drop_unused_indices(out, nq.indices))
end

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
        new_c = _expand_prefactor(c; kwargs...)
        _iszero_cnum(new_c) && continue
        _addto!(d, term.ops, new_c, term.ne)
    end
    return QAdd(d, copy(s.indices))
end
Symbolics.expand(op::QSym; kwargs...) = _single_qadd(_CNUM_ONE, QSym[op])

_expand_prefactor(x::CNum; kwargs...) = _iszero_cnum(x) ? x : Complex(Symbolics.expand(real(x)), Symbolics.expand(imag(x)))

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
1 + -σ₂₂
```

See also [`assume_distinct_index`](@ref), [`normal_order`](@ref).
"""
function expand_completeness(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        _expand_gs_ops(copy(t.ops), c) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    return QAdd(out, copy(q.indices))
end

expand_completeness(op::QSym) = expand_completeness(_single_qadd(_CNUM_ONE, QSym[op]))

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
            ne = _merge_ne_pair(ne, a, b)
        end
        _canonicalize!(out, copy(t.ops), c, ne)
    end
    return expand_completeness(QAdd(out, copy(q.indices)))
end

const _ZERO_QADD = QAdd(QTermDict(), Index[])
_zero_qadd() = _ZERO_QADD

# Different sites commute, so `[a, b] = 0` without routing through `*`.
@inline _same_site(a::QSym, b::QSym) =
    a.space_index == b.space_index && a.index == b.index

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

commutator(::Number, ::Number) = _zero_qadd()
commutator(::Number, ::QField) = _zero_qadd()
commutator(::QField, ::Number) = _zero_qadd()

function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _zero_qadd()
    isequal(a, b) && return _zero_qadd()
    # If exactly one direction needs a swap, [a, b] is the swap residual.
    forward = _can_commute(a, b)
    reverse = _can_commute(b, a)
    if !forward && reverse
        _, _, c, ops = _commute_pair(a, b)
        return _single_qadd(c, isempty(ops) ? _EMPTY_OPS : copy(ops))
    elseif forward && !reverse
        _, _, c, ops = _commute_pair(b, a)
        return _single_qadd(_neg_cnum(c), isempty(ops) ? _EMPTY_OPS : copy(ops))
    end
    return a * b - b * a
end

commutator(a::QAdd, b::QSym) = a * b - b * a
commutator(a::QSym, b::QAdd) = a * b - b * a

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return _zero_qadd()
    return a * b - b * a
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
    substitute(expr, d::Dict)

Substitute symbolic parameters and/or operators in `expr` using dictionary `d`.

Supports both symbolic coefficient replacement (for example `x => 2`) and
operator replacement. The result is re-canonicalized and returned as a
[`QAdd`](@ref).

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
function SymbolicUtils.substitute(op::QSym, d::Dict)
    return SymbolicUtils.substitute(_single_qadd(_CNUM_ONE, QSym[op]), d)
end

function SymbolicUtils.substitute(q::QAdd, d::Dict)
    out = QTermDict()
    for (t, c) in q
        _substitute_ops(copy(t.ops), c, d) do ops1, c1
            _stream!(out, ops1, c1, t.ne)
        end
    end
    return QAdd(out, copy(q.indices))
end
