# ============================================================================
#  algebra.jl — user-facing operator algebra
#
#  Defines the public API on top of the canonicalization passes (passes.jl) and
#  pipeline composition (pipelines.jl): scalar arithmetic, addition,
#  subtraction, powers, normal_order, simplify, commutator, expand,
#  expand_completeness, substitute.
# ============================================================================

# ---------------------------------------------------------------------------
#  Scalar multiplication
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
#  Addition
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
#  Subtraction, negation, division, power
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
#  normal_order
# ---------------------------------------------------------------------------

"""
    normal_order(expr::QField)

Apply normal ordering. Since `*` already produces canonical form, this is
typically an idempotent finalizer that re-routes each term through the
canonicalization pipeline.

For ground-state projector expansion (`σᵍᵍ = 1 - Σ σᵏᵏ`), use
[`expand_completeness`](@ref) instead.

See also [`simplify`](@ref), [`normal_to_symmetric`](@ref),
[`symmetric_to_normal`](@ref), [`expand_completeness`](@ref).
"""
normal_order(op::QSym) = _single_qadd(_CNUM_ONE, QSym[op])

function normal_order(q::QAdd)
    out = QTermDict()
    for (t, c) in q
        _stream!(out, copy(t.ops), c, t.ne)
    end
    return QAdd(out, copy(q.indices))
end

# ---------------------------------------------------------------------------
#  simplify
# ---------------------------------------------------------------------------

function _simplify_prefactor(x::CNum)
    _const_val(x.re) !== nothing && _const_val(x.im) !== nothing && return x
    re = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(Symbolics.expand(real(x)))))
    im = Num(SymbolicUtils.simplify(SymbolicUtils.unwrap(Symbolics.expand(imag(x)))))
    return Complex(re, im)
end
_simplify_prefactor(x::Number) = x

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
    simplify(expr::QField)

Apply normal ordering + algebraic identities + per-coefficient symbolic
simplification, then drop summation indices that no surviving term depends on.

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

# ---------------------------------------------------------------------------
#  expand (prefactor expansion)
# ---------------------------------------------------------------------------

"""
    expand(expr::QField)

Expand symbolic prefactors in each term, e.g. `(a+b)^2 → a^2 + 2ab + b^2`.
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

_expand_prefactor(x::CNum; kwargs...) = _iszero_cnum(x) ? x : Symbolics.expand(x; kwargs...)
_expand_prefactor(x::Number; kwargs...) = x

# ---------------------------------------------------------------------------
#  expand_completeness (opt-in σᵍᵍ rewriting)
# ---------------------------------------------------------------------------

"""
    expand_completeness(q) -> QAdd

Apply `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ` to every ground-state projector in `q` and
re-canonicalize each branch. Opt-in; `*` keeps σᵍᵍ atomic.
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

Re-canonicalize `q` under the given pairwise `≠` constraints on free indices,
then run [`expand_completeness`](@ref) so any ground-state projectors that
emerge after the constraint resolves same-site composition are expanded.

Use this when free indices `j` and `k` semantically range over distinct atoms
or modes but no `Σ` supplies the constraint. The pipeline cannot infer
"different symbol → different site" automatically, so this is the explicit
way to declare it.

```julia
# Declare j and k as distinct atomic sites:
assume_distinct_index([H, σ¹²_j · σ²¹_k], [(j, k)])
```
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

# ---------------------------------------------------------------------------
#  commutator / anticommutator
# ---------------------------------------------------------------------------

const _ZERO_QADD = QAdd(QTermDict(), Index[])
_zero_qadd() = _ZERO_QADD

# Two leaves are on the same site iff their (space_index, index) agree.
# Different sites commute → `[a, b] = 0` without routing through `*`.
@inline _same_site(a::QSym, b::QSym) =
    a.space_index == b.space_index && a.index == b.index

"""
    commutator(a, b) -> QAdd

Compute `[a, b] = a*b - b*a`.

For indexed expressions, the diagonal substitutions performed by `*` (via
`_accumulate_with_diag!`) produce the partial-collapse contributions when
summation indices share a Hilbert subspace with the other factor.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
commutator(a, a')    # 1
commutator(a', a)    # -1
```
"""
function commutator end

commutator(::Number, ::Number) = _zero_qadd()
commutator(::Number, ::QField) = _zero_qadd()
commutator(::QField, ::Number) = _zero_qadd()

function commutator(a::QSym, b::QSym)
    _same_site(a, b) || return _zero_qadd()
    isequal(a, b) && return _zero_qadd()
    # Fast path: if exactly one direction needs a swap, [a, b] is the swap residual.
    ab = _can_commute(a, b)
    ba = _can_commute(b, a)
    if !ab && ba
        _, _, c, ops = _commute_pair(a, b)
        return _single_qadd(c, isempty(ops) ? _EMPTY_OPS : copy(ops))
    elseif ab && !ba
        _, _, c, ops = _commute_pair(b, a)
        return _single_qadd(_neg_cnum(c), isempty(ops) ? _EMPTY_OPS : copy(ops))
    end
    # Both directions non-commuting (Pauli, Transition) — fall back to a*b - b*a.
    return a * b - b * a
end

commutator(a::QAdd, b::QSym) = a * b - b * a
commutator(a::QSym, b::QAdd) = a * b - b * a

function commutator(a::QAdd, b::QAdd)
    isequal(a, b) && return _zero_qadd()
    return a * b - b * a
end

"""
    anticommutator(a, b)

Compute `{a, b} = a*b + b*a`.
"""
anticommutator(a, b) = a * b + b * a

# ---------------------------------------------------------------------------
#  substitute (public API; pass lives in passes.jl)
# ---------------------------------------------------------------------------

"""
    substitute(expr, d::Dict)

Substitute symbolic parameters and/or operators in `expr` according to `d`.
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
