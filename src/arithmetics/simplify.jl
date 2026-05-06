# --- Algebraic reductions (ordering-independent) ---

"""
    _apply_reductions(arg_c, ops) -> Vector{OrderedTerm}

Apply ordering-independent algebraic identities to a product of operators
using a worklist algorithm. Only applies:
- Transition composition: |i⟩⟨j| · |k⟩⟨l| → |i⟩⟨l| if j==k, else 0
- Pauli product rule: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ

Never applies commutation-based reordering (Fock, Spin, PhaseSpace).
"""
function _apply_reductions(arg_c::CNum, ops::Vector{QSym})
    isempty(ops) && return OrderedTerm[OrderedTerm(arg_c, ops)]
    length(ops) == 1 && return OrderedTerm[OrderedTerm(arg_c, ops)]

    worklist = OrderedTerm[OrderedTerm(arg_c, copy(ops))]
    done = OrderedTerm[]

    while !isempty(worklist)
        t = pop!(worklist)
        _reduce_product!(t.prefactor, t.ops, worklist, done)
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
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm}
    )
    n = length(ops)

    if n <= 1
        push!(done, OrderedTerm(c, ops))
        return
    end

    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]
        # `simplify(expr)` is the LazyOrder reductions-only pass: keep
        # σᵍᵍ atomic. Use `simplify(expr, h)` to opt into completeness.
        _try_algebraic_reduction!(a, b, i, c, ops, worklist, done, false) && return
    end

    return push!(done, OrderedTerm(c, ops))
end

# --- Simplification (ordering-independent) ---

function _qsimplify(op::QSym)
    return _single_qadd(_CNUM_ONE, QSym[op])
end

function _qsimplify(s::QAdd)
    d = QTermDict()
    for (term, c) in s.arguments
        _iszero_cnum(c) && continue
        for t in _apply_reductions(c, term.ops)
            _addto!(d, t.ops, t.prefactor, term.ne)
        end
    end
    # Drop summation indices that no surviving term depends on.
    indices = _drop_unused_indices(d, s.indices)
    return QAdd(d, indices)
end

"""
    _drop_unused_indices(d, indices) -> Vector{Index}

Filter `indices` to those some term in `d` actually depends on. Algebraic
reductions can collapse all index-bearing terms (e.g. `σ¹² · σ³¹ = 0`),
leaving the surrounding `Σ` over a now-vacuous index — this drops it so
the resulting `QAdd` doesn't claim a phantom summation.
"""
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
    d = QTermDict()
    for (term, c) in s.arguments
        new_c = _expand_prefactor(c; kwargs...)
        _iszero_cnum(new_c) && continue
        _addto!(d, term.ops, new_c, term.ne)
    end
    return QAdd(d, copy(s.indices))
end
function Symbolics.expand(op::QSym; kwargs...)
    return _single_qadd(_CNUM_ONE, QSym[op])
end

_expand_prefactor(x::CNum; kwargs...) = _iszero_cnum(x) ? x : Symbolics.expand(x; kwargs...)
_expand_prefactor(x::Number; kwargs...) = x
