"""
    simplify(expr, [ordering])

Simplify a quantum expression by:
1. Applying ordering-independent algebraic reductions (Transition composition/orthogonality, Pauli product rule)
2. Applying ordering-dependent commutation swaps (default: `NormalOrder()`)
3. Collecting like terms and summing prefactors

Returns `QAdd`. The `ordering` argument controls which commutation rules are applied.
Default is `NormalOrder()`.
"""

"""
    _same_site(a::QSym, b::QSym) -> Bool

Two operators are on the same site iff they share space_index, copy_index, AND index.
Only operators on the same site can interact (apply commutation/composition rules).
"""
function _same_site(a::QSym, b::QSym)
    return a.space_index == b.space_index &&
           a.copy_index == b.copy_index &&
           a.index == b.index
end

# Internal quantum simplification — applies operator algebra rules
function _qsimplify(op::QSym, ::OrderingConvention)
    return QAdd(QMul[QMul(_to_cnum(1), QSym[op])])
end

function _qsimplify(m::QMul, ord::OrderingConvention)
    terms = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return _collect_like_terms(QAdd(terms))
end

function _qsimplify(s::QAdd, ord::OrderingConvention)
    all_terms = QMul[]
    for term in s.arguments
        append!(all_terms, _simplify_qmul(term.arg_c, term.args_nc, ord))
    end
    return _collect_like_terms(QAdd(all_terms, s.indices, s.non_equal))
end

# With Hilbert space (ground state rewriting)
function _qsimplify(op::QSym, ord::OrderingConvention, h::HilbertSpace)
    return _collect_like_terms(_apply_ground_state(_qsimplify(op, ord), h))
end
function _qsimplify(m::QMul, ord::OrderingConvention, h::HilbertSpace)
    terms = _simplify_qmul(m.arg_c, m.args_nc, ord)
    return _collect_like_terms(_apply_ground_state(QAdd(terms), h))
end
function _qsimplify(s::QAdd, ord::OrderingConvention, h::HilbertSpace)
    all_terms = QMul[]
    for term in s.arguments
        append!(all_terms, _simplify_qmul(term.arg_c, term.args_nc, ord))
    end
    return _collect_like_terms(_apply_ground_state(QAdd(all_terms, s.indices, s.non_equal), h))
end

# Public API: SymbolicUtils.simplify dispatches to _qsimplify
SymbolicUtils.simplify(expr::QField; kwargs...) = _qsimplify(expr, NormalOrder())
SymbolicUtils.simplify(expr::QField, h::HilbertSpace; kwargs...) = _qsimplify(expr, NormalOrder(), h)

# Levi-Civita lookup: _levi_civita[j][k] = ε_{jk(6-j-k)} for j,k ∈ {1,2,3}
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

"""
    _simplify_qmul(prefactor, ops, ordering) -> Vector{QMul}

Simplify a product of operators using a worklist algorithm.
"""
function _simplify_qmul(arg_c::CNum, ops::Vector{QSym}, ord::OrderingConvention)
    isempty(ops) && return QMul[QMul(arg_c, QSym[])]
    length(ops) == 1 && return QMul[QMul(arg_c, copy(ops))]

    # copy(ops) is critical: _simplify_product! mutates the ops vector in-place
    # (deleteat!, element swaps). Without this copy, the caller's vector is corrupted.
    worklist = QMul[QMul(arg_c, copy(ops))]
    done = QMul[]

    while !isempty(worklist)
        term = pop!(worklist)
        _simplify_product!(term, ord, worklist, done)
    end

    return done
end

"""
    _simplify_product!(m::QMul, ord, worklist, done)

Process one QMul term. Scans adjacent operator pairs for the first applicable
rule, pushes resulting terms to `worklist`, or pushes to `done` if fully simplified.
"""
function _simplify_product!(m::QMul, ord::OrderingConvention,
    worklist::Vector{QMul}, done::Vector{QMul})
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
        if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
            if a.j == b.i
                ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.copy_index, a.index)
                deleteat!(ops, i + 1)
                push!(worklist, QMul(c, ops))
            else
                push!(done, QMul(_to_cnum(0), QSym[]))
            end
            return
        end

        # Pauli: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
        if a isa Pauli && b isa Pauli && _same_site(a, b) && a.name == b.name
            if a.axis == b.axis
                deleteat!(ops, i:(i + 1))
                push!(worklist, QMul(c, ops))
            else
                eps = _levi_civita[a.axis][b.axis]
                ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index, a.index)
                deleteat!(ops, i + 1)
                push!(worklist, QMul(c * _to_cnum(im * eps), ops))
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
    _apply_ordering_swap!(a, b, i, c, ops, ord, worklist) -> Bool

Apply ordering-specific commutation rules in-place.
"""
function _apply_ordering_swap!(a, b, i, c::CNum, ops,
    ::NormalOrder, worklist::Vector{QMul})
    # Fock: a·a† → a†·a + 1
    if a isa Destroy && b isa Create && _same_site(a, b) && a.name == b.name
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, QMul(c, ops))
        return true
    end

    # Spin: [Sⱼ, Sₖ] = iϵⱼₖₗSₗ (swap out-of-order axes)
    if a isa Spin && b isa Spin && _same_site(a, b) && a.name == b.name && a.axis > b.axis
        eps = _levi_civita[a.axis][b.axis]
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index, a.index)
        deleteat!(ops, i + 1)
        push!(worklist, QMul(c * _to_cnum(im * eps), ops))
        return true
    end

    # PhaseSpace: P·X → X·P - i
    if a isa Momentum && b isa Position && _same_site(a, b)
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, QMul(c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, QMul(_to_cnum(-im) * c, ops))
        return true
    end

    return false
end

"""
    _collect_like_terms(expr::QAdd)

Group terms with identical `args_nc`, sum their `arg_c` prefactors.
"""
function _collect_like_terms(s::QAdd)
    d = Dict{Vector{QSym}, CNum}()
    for term in s.arguments
        d[term.args_nc] = get(d, term.args_nc, _to_cnum(0)) + term.arg_c
    end

    result = QMul[]
    for (ops, c) in d
        iszero(c) && continue
        cs = _simplify_prefactor(c)
        iszero(cs) && continue
        push!(result, QMul(cs, ops))
    end
    isempty(result) && return QAdd(QMul[QMul(_to_cnum(0), QSym[])], s.indices, s.non_equal)
    return QAdd(result, s.indices, s.non_equal)
end

# Symbolics.expand — expand symbolic prefactors (e.g. (a+b)² → a²+2ab+b²)
# Note: operator products are already distributed at the QAdd*QAdd level;
# this only calls Symbolics.expand on the c-number prefactors.
function Symbolics.expand(s::QAdd; kwargs...)
    result = QMul[
        QMul(_expand_prefactor(t.arg_c; kwargs...), t.args_nc)
            for t in s.arguments
    ]
    return QAdd(result, s.indices, s.non_equal)
end
function Symbolics.expand(m::QMul; kwargs...)
    return Symbolics.expand(QAdd(QMul[m]); kwargs...)
end
function Symbolics.expand(op::QSym; kwargs...)
    return QAdd(QMul[QMul(_to_cnum(1), QSym[op])])
end

function _simplify_prefactor(x::CNum)
    iszero(x) && return x
    r, i = Symbolics.unwrap(real(x)), Symbolics.unwrap(imag(x))
    # Fast path: skip Symbolics.simplify for purely numeric prefactors
    (SymbolicUtils.iscall(r) || SymbolicUtils.issym(r) ||
     SymbolicUtils.iscall(i) || SymbolicUtils.issym(i)) || return x
    return Symbolics.simplify(x)
end
_simplify_prefactor(x::Number) = x

_expand_prefactor(x::CNum; kwargs...) = iszero(x) ? x : Symbolics.expand(x; kwargs...)
_expand_prefactor(x::Number; kwargs...) = x
