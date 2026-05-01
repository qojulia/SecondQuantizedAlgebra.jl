# --- Algebraic reductions (ordering-independent) ---

"""
    _apply_reductions(arg_c, ops) -> Vector{_OTerm}

Apply ordering-independent algebraic identities to a product of operators
using a worklist algorithm. Only applies:
- Transition composition: |i⟩⟨j| · |k⟩⟨l| → |i⟩⟨l| if j==k, else 0
- Pauli product rule: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ

Never applies commutation-based reordering (Fock, Spin, PhaseSpace).
"""
function _apply_reductions(arg_c::CNum, ops::Vector{QSym})
    isempty(ops) && return _OTerm[(arg_c, ops)]
    length(ops) == 1 && return _OTerm[(arg_c, ops)]

    worklist = _OTerm[(arg_c, copy(ops))]
    done = _OTerm[]

    while !isempty(worklist)
        c, ops_w = pop!(worklist)
        _reduce_product!(c, ops_w, worklist, done)
    end

    return done
end

"""
    _reduce_product!(c, ops, worklist, done)

Process one term. Scans adjacent operator pairs for algebraic identities only.
Pushes resulting terms to `worklist`, or to `done` if fully reduced.
"""
function _reduce_product!(
        c::CNum, ops::Vector{QSym},
        worklist::Vector{_OTerm}, done::Vector{_OTerm}
    )
    n = length(ops)

    if n <= 1
        push!(done, (c, ops))
        return
    end

    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]
        # `simplify(expr)` is the LazyOrder reductions-only pass: keep
        # σᵍᵍ atomic. Use `simplify(expr, h)` to opt into completeness.
        _try_algebraic_reduction!(a, b, i, c, ops, worklist, done, false) && return
    end

    return push!(done, (c, ops))
end

# --- Simplification (ordering-independent) ---

function _qsimplify(op::QSym)
    return QAdd(QTermDict(QSym[op] => _CNUM_ONE), Index[], Tuple{Index, Index}[])
end

function _qsimplify(s::QAdd)
    d = QTermDict()
    for (ops, c) in s.arguments
        _iszero_cnum(c) && continue
        for (oc, oops) in _apply_reductions(c, ops)
            _addto!(d, oops, oc)
        end
    end
    return QAdd(d, s.indices, s.non_equal)
end

# With Hilbert space (ground state rewriting)
function _qsimplify(op::QSym, h::HilbertSpace)
    return _apply_ground_state(_qsimplify(op), h)
end
function _qsimplify(s::QAdd, h::HilbertSpace)
    return _apply_ground_state(_qsimplify(s), h)
end

# Public API

"""
    simplify(expr::QField)
    simplify(expr::QField, h::HilbertSpace)

Simplify a quantum operator expression by applying ordering-independent algebraic identities.

Applied rules:
- **Transition composition**: ``|i\\rangle\\langle j| \\cdot |k\\rangle\\langle l| = \\delta_{jk} |i\\rangle\\langle l|``
- **Pauli product rule**: ``\\sigma_j \\sigma_k = \\delta_{jk} I + i\\epsilon_{jkl} \\sigma_l``
- Like-term collection and zero-term elimination

Does **not** apply commutation-based reordering (Fock `[a, a†]`, Spin `[Sⱼ, Sₖ]`,
PhaseSpace `[p, x]`). For that, use [`normal_order`](@ref).

The `h`-argument overload is the [`LazyOrder`](@ref) opt-in for ground-state
completeness: it additionally rewrites ground-state projectors
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|`` of every
[`NLevelSpace`](@ref) subspace. Without `h`, ``|g\\rangle\\langle g|`` is kept
atomic — useful for inspecting the un-expanded form.

Under [`NormalOrder`](@ref) (the default ordering convention) the eager arithmetic
already produces canonical form including completeness, so `simplify` is typically
idempotent and the `h`-overload is a no-op cleanup pass.

See also [`normal_order`](@ref), [`expand`](@ref).
"""
SymbolicUtils.simplify(expr::QField; kwargs...) = _qsimplify(expr)
SymbolicUtils.simplify(expr::QField, h::HilbertSpace; kwargs...) = _qsimplify(expr, h)

"""
    expand(expr::QField)

Expand symbolic prefactors in quantum operator expressions.

Applies `Symbolics.expand` to each prefactor in a [`QAdd`](@ref) sum,
e.g. expanding `(a+b)^2` into `a^2 + 2ab + b^2`. Zero terms are dropped after expansion.

See also [`simplify`](@ref).
"""
function Symbolics.expand(s::QAdd; kwargs...)
    d = QTermDict(ops => _expand_prefactor(c; kwargs...) for (ops, c) in s.arguments)
    _drop_zeros!(d)
    return QAdd(d, s.indices, s.non_equal)
end
function Symbolics.expand(op::QSym; kwargs...)
    return QAdd(QTermDict(QSym[op] => _CNUM_ONE), Index[], Tuple{Index, Index}[])
end

_expand_prefactor(x::CNum; kwargs...) = _iszero_cnum(x) ? x : Symbolics.expand(x; kwargs...)
_expand_prefactor(x::Number; kwargs...) = x
