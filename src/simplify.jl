"""
    simplify(expr, [ordering])

Simplify a quantum expression by:
1. Applying ordering-independent algebraic reductions (Transition composition/orthogonality, Pauli product rule)
2. Applying ordering-dependent commutation swaps (default: `NormalOrder()`)
3. Collecting like terms and summing prefactors

Returns `QAdd`. The `ordering` argument controls which commutation rules are applied.
Default is `NormalOrder()`.

The output element type is always widened to `Complex{T}` for real `T`, so that
commutation rules involving `im` (Pauli, Spin, PhaseSpace) are type-stable.
"""
simplify(expr::QField) = simplify(expr, NormalOrder())
simplify(expr::QField, h::HilbertSpace) = simplify(expr, NormalOrder(), h)

# Widen real types to Complex for type-stable simplification.
# Commutation rules for Pauli, Spin, and PhaseSpace multiply by `im`,
# so the output type must be able to hold complex values.
_complex_promote(::Type{T}) where {T <: Real} = Complex{T}
_complex_promote(::Type{T}) where {T} = T  # already Complex or symbolic

function simplify(op::QSym, ::OrderingConvention)
    CT = Complex{Int}
    return QAdd(QMul{CT}[QMul(one(CT), QSym[op])])
end

function simplify(m::QMul{T}, ord::OrderingConvention) where {T}
    CT = _complex_promote(T)
    terms = _simplify_qmul(convert(CT, m.arg_c), m.args_nc, ord)
    return _collect_like_terms(QAdd(terms))
end

function simplify(s::QAdd{T}, ord::OrderingConvention) where {T}
    CT = _complex_promote(T)
    all_terms = QMul{CT}[]
    for term in s.arguments
        append!(all_terms, _simplify_qmul(convert(CT, term.arg_c), term.args_nc, ord))
    end
    return _collect_like_terms(QAdd(all_terms))
end

# With Hilbert space (ground state rewriting)
function simplify(op::QSym, ord::OrderingConvention, h::HilbertSpace)
    return _collect_like_terms(_apply_ground_state(simplify(op, ord), h))
end
function simplify(m::QMul{T}, ord::OrderingConvention, h::HilbertSpace) where {T}
    CT = _complex_promote(T)
    terms = _simplify_qmul(convert(CT, m.arg_c), m.args_nc, ord)
    return _collect_like_terms(_apply_ground_state(QAdd(terms), h))
end
function simplify(s::QAdd{T}, ord::OrderingConvention, h::HilbertSpace) where {T}
    CT = _complex_promote(T)
    all_terms = QMul{CT}[]
    for term in s.arguments
        append!(all_terms, _simplify_qmul(convert(CT, term.arg_c), term.args_nc, ord))
    end
    return _collect_like_terms(_apply_ground_state(QAdd(all_terms), h))
end

# Levi-Civita lookup: _levi_civita[j][k] = ε_{jk(6-j-k)} for j,k ∈ {1,2,3}
# Zero-allocation, constant-time replacement for Combinatorics.levicivita.
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

"""
    _simplify_qmul(prefactor, ops, ordering) -> Vector{QMul{CT}}

Simplify a product of operators using a worklist algorithm. Returns a flat
`Vector{QMul{CT}}` of fully simplified terms (not wrapped in QAdd).

The worklist and done vectors are typed `Vector{QMul{CT}}` where `CT` is the
(pre-widened) prefactor type, ensuring full type stability throughout.
"""
function _simplify_qmul(arg_c::CT, ops::Vector{QSym}, ord::OrderingConvention) where {CT}
    isempty(ops) && return QMul{CT}[QMul(arg_c, QSym[])]
    length(ops) == 1 && return QMul{CT}[QMul(arg_c, copy(ops))]

    worklist = QMul{CT}[QMul(arg_c, copy(ops))]
    done = QMul{CT}[]

    while !isempty(worklist)
        term = pop!(worklist)
        _simplify_product!(term, ord, worklist, done)
    end

    return done
end

"""
    _simplify_product!(m::QMul{CT}, ord, worklist, done)

Process one QMul term. Scans adjacent operator pairs for the first applicable
rule, pushes resulting terms to `worklist`, or pushes to `done` if fully simplified.
Mutates `m.args_nc` in-place for efficiency.
"""
function _simplify_product!(
        m::QMul{CT}, ord::OrderingConvention,
        worklist::Vector{QMul{CT}}, done::Vector{QMul{CT}}
    ) where {CT}
    ops = m.args_nc
    c = m.arg_c
    n = length(ops)

    if n <= 1
        push!(done, m)
        return
    end

    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]

        ## Ordering-independent reductions

        # Transition: |i⟩⟨j| · |k⟩⟨l|
        if a isa Transition && b isa Transition && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
            if a.j == b.i
                ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.copy_index)
                deleteat!(ops, i + 1)
                push!(worklist, QMul(c, ops))
            else
                push!(done, QMul(zero(c), QSym[]))
            end
            return
        end

        # Pauli: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
        if a isa Pauli && b isa Pauli && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
            if a.axis == b.axis
                deleteat!(ops, i:(i + 1))
                push!(worklist, QMul(c, ops))
            else
                eps = _levi_civita[a.axis][b.axis]
                ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index)
                deleteat!(ops, i + 1)
                push!(worklist, QMul(c * im * eps, ops))
            end
            return
        end

        ## Ordering-dependent swaps
        if _apply_ordering_swap!(a, b, i, c, ops, ord, worklist)
            return
        end
    end

    # No rule applied — fully simplified
    return push!(done, m)
end

"""
    _apply_ordering_swap!(a, b, i, c::CT, ops, ord, worklist) -> Bool

Apply ordering-specific commutation rules in-place. Returns `true` if a rule
was applied (new terms pushed to `worklist`), `false` otherwise.

Same-space adjacent swaps preserve space_index ordering, so no re-sorting needed.
"""
function _apply_ordering_swap!(
        a, b, i, c::CT, ops,
        ::NormalOrder, worklist::Vector{QMul{CT}}
    ) where {CT}
    # Fock: a·a† → a†·a + 1
    if a isa Destroy && b isa Create && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, QMul(c, ops))
        return true
    end

    # Spin: [Sⱼ, Sₖ] = iϵⱼₖₗSₗ (swap out-of-order axes)
    if a isa Spin && b isa Spin && a.space_index == b.space_index && a.copy_index == b.copy_index && a.name == b.name && a.axis > b.axis
        eps = _levi_civita[a.axis][b.axis]
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index)
        deleteat!(ops, i + 1)
        push!(worklist, QMul(c * im * eps, ops))
        return true
    end

    # PhaseSpace: P·X → X·P - i
    # Position and Momentum have different names, so check space_index only
    if a isa Momentum && b isa Position && a.space_index == b.space_index && a.copy_index == b.copy_index
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, QMul(-im * c, ops))
        return true
    end

    return false
end

# TODO: SymmetricOrder for TWA (Truncated Wigner Approximation)
# Symmetric ordering symmetrizes operator products: (ab + ba)/2.
# For Fock: instead of swapping a·a† → a†·a + 1, symmetric ordering
# would use (a·a† + a†·a)/2 = a†·a + 1/2.
# For Spin: symmetric ordering of Sⱼ·Sₖ = (SⱼSₖ + SₖSⱼ)/2.
# Implement by adding:
#   function _apply_ordering_swap!(a, b, i, c::CT, ops,
#       ::SymmetricOrder, worklist::Vector{QMul{CT}}) where {CT}
#       ...
#   end

"""
    _collect_like_terms(expr::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
Remove zero-prefactor terms. Fully type-stable: input QAdd{T} → output QAdd{T}.
Uses a Dict for O(n) amortized collection instead of O(n²) linear scan.
"""
function _collect_like_terms(s::QAdd{T}) where {T}
    d = Dict{Vector{QSym}, T}()
    for term in s.arguments
        d[term.args_nc] = get(d, term.args_nc, zero(T)) + term.arg_c
    end

    result = QMul{T}[QMul(c, ops) for (ops, c) in d if !iszero(c)]
    isempty(result) && return QAdd(QMul{T}[QMul(zero(T), QSym[])])
    return QAdd(result)
end

# SymbolicUtils.simplify — also simplify each prefactor
function SymbolicUtils.simplify(s::QAdd; kwargs...)
    simplified = simplify(s)
    args = simplified.arguments
    TT = promote_type((typeof(_simplify_prefactor(t.arg_c; kwargs...)) for t in args)...)
    result = QMul{TT}[
        QMul(convert(TT, _simplify_prefactor(t.arg_c; kwargs...)), t.args_nc)
            for t in args
    ]
    return QAdd(result)
end
function SymbolicUtils.simplify(m::QMul; kwargs...)
    return SymbolicUtils.simplify(QAdd([m]); kwargs...)
end
function SymbolicUtils.simplify(op::QSym; kwargs...)
    return SymbolicUtils.simplify(simplify(op); kwargs...)
end

_simplify_prefactor(x::Number; kwargs...) = x
_simplify_prefactor(x; kwargs...) = SymbolicUtils.simplify(x; kwargs...)

# Symbolics.expand — distribute products, then expand prefactors
function Symbolics.expand(s::QAdd; kwargs...)
    args = s.arguments
    TT = promote_type((typeof(_expand_prefactor(t.arg_c; kwargs...)) for t in args)...)
    result = QMul{TT}[
        QMul(convert(TT, _expand_prefactor(t.arg_c; kwargs...)), t.args_nc)
            for t in args
    ]
    return QAdd(result)
end
function Symbolics.expand(m::QMul; kwargs...)
    return Symbolics.expand(QAdd([m]); kwargs...)
end
function Symbolics.expand(op::QSym; kwargs...)
    return QAdd(QMul{Int}[QMul(1, QSym[op])])
end

_expand_prefactor(x::Number; kwargs...) = x
_expand_prefactor(x; kwargs...) = Symbolics.expand(x; kwargs...)
