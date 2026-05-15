"""
    _same_site(a::QSym, b::QSym) -> Bool

Two operators are on the same site iff they share space_index AND index.
Only operators on the same site can interact (apply commutation/composition rules).
"""
function _same_site(a::QSym, b::QSym)
    return a.space_index == b.space_index && a.index == b.index
end

# Levi-Civita lookup: _levi_civita[j][k] = ε_{jk(6-j-k)} for j,k ∈ {1,2,3}
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

# Worklist term: a (prefactor, operators) record used by the ordering and
# simplification worklists, and by the Weyl conversion.
struct OrderedTerm
    prefactor::CNum
    ops::Vector{QSym}
end

# Per-rule step return:
#   nothing  — no rule applied; caller scans the next pair
#   c::CNum  — single-output rule applied in place; caller resumes scan with c
#   HALT     — forking or terminal rule applied; outputs already pushed; caller stops
struct Halt end
const HALT = Halt()
const StepResult = Union{Nothing, Halt, CNum}

"""
    _try_reduction!(a, b, i, c, ops, worklist, done) -> StepResult

Apply an algebraic identity to the pair `(a, b)` at position `i`. A Transition
composition that would produce a ground-state projector ``|g\\rangle\\langle g|``
is rewritten via the completeness relation
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g} |k\\rangle\\langle k|`` directly,
keeping the result in canonical form.
"""
@inline function _try_reduction!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm}
    )::StepResult
    # Transition: |i⟩⟨j| · |k⟩⟨l|
    if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
        if a.j != b.i
            push!(done, OrderedTerm(_to_cnum(0), QSym[]))
            return HALT
        end
        if a.i == a.ground_state && b.j == a.ground_state
            _push_ground_state_expansion!(c, ops, i, a, worklist)
            return HALT
        end
        ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.index, a.ground_state, a.n_levels)
        deleteat!(ops, i + 1)
        return c
    end

    # Pauli: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
    if a isa Pauli && b isa Pauli && _same_site(a, b) && a.name == b.name
        if a.axis == b.axis
            deleteat!(ops, i:(i + 1))
            return c
        end
        eps = _levi_civita[a.axis][b.axis]
        ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
        deleteat!(ops, i + 1)
        return _mul_cnum(c, _to_cnum(im * eps))
    end

    return nothing
end

"""
    _push_ground_state_expansion!(c, ops, i, a, worklist)

Expand a ground-state projector that would be produced at position `i` (from
`ops[i] · ops[i+1]`) using completeness ``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g}|k\\rangle\\langle k|``.
Removes the pair `(ops[i], ops[i+1])` from `ops` and pushes one identity term plus
``n_{\\text{levels}} - 1`` expansion terms onto `worklist`. Uses `a.ground_state`,
`a.n_levels`, `a.name`, `a.space_index`, `a.index` to construct the replacement
projectors.
"""
function _push_ground_state_expansion!(
        c::CNum, ops::Vector{QSym}, i::Int, a::Transition, worklist::Vector{OrderedTerm}
    )
    n = length(ops)
    # Identity contribution: drop ops[i] and ops[i+1]
    id_ops = Vector{QSym}(undef, n - 2)
    @inbounds for k in 1:(i - 1)
        id_ops[k] = ops[k]
    end
    @inbounds for k in (i + 2):n
        id_ops[k - 2] = ops[k]
    end
    push!(worklist, OrderedTerm(c, id_ops))
    # -σᵏᵏ for each k ≠ ground_state, replacing the pair at position i
    neg_c = _neg_cnum(c)
    for k in 1:a.n_levels
        k == a.ground_state && continue
        new_ops = Vector{QSym}(undef, n - 1)
        @inbounds for kk in 1:(i - 1)
            new_ops[kk] = ops[kk]
        end
        new_ops[i] = Transition(a.name, k, k, a.space_index, a.index, a.ground_state, a.n_levels)
        @inbounds for kk in (i + 2):n
            new_ops[kk - 1] = ops[kk]
        end
        push!(worklist, OrderedTerm(neg_c, new_ops))
    end
    return nothing
end

# Push the commuted (swapped) branch onto `worklist` so the caller can continue
# in-place on the contracted branch with a single `deleteat!` / mutation.
@inline function _push_swapped!(
        ops::Vector{QSym}, i::Int, c::CNum, worklist::Vector{OrderedTerm}
    )
    swapped = copy(ops)
    swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
    push!(worklist, OrderedTerm(c, swapped))
    return nothing
end

"""
    _try_swap!(a, b, i, c, ops, worklist) -> StepResult

Apply a commutation rule to the pair `(a, b)` at position `i`. Forking rules
push the swapped branch onto `worklist` and return the new prefactor for
in-place continuation on the contracted branch; `_try_swap!` never returns
`HALT`.
"""
@inline function _try_swap!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        worklist::Vector{OrderedTerm}
    )::StepResult
    # Fock: a·a† = a†·a + 1
    if a isa Destroy && b isa Create && _same_site(a, b) && a.name == b.name
        _push_swapped!(ops, i, c, worklist)
        deleteat!(ops, i:(i + 1))
        return c
    end

    # Spin: [Sj, Sk] = i·ε_jkl·Sl (swap out-of-order axes)
    if a isa Spin && b isa Spin && _same_site(a, b) && a.name == b.name && a.axis > b.axis
        _push_swapped!(ops, i, c, worklist)
        eps = _levi_civita[a.axis][b.axis]
        ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
        deleteat!(ops, i + 1)
        return _mul_cnum(c, _to_cnum(im * eps))
    end

    # PhaseSpace: P·X = X·P − i
    if a isa Momentum && b isa Position && _same_site(a, b)
        _push_swapped!(ops, i, c, worklist)
        deleteat!(ops, i:(i + 1))
        return _mul_cnum(_CNUM_NEG_IM, c)
    end

    return nothing
end

"""
    _process_product!(c, ops, worklist, done)

Drive one term `(c, ops)` to completion. Single-output rules mutate `ops` in
place and the loop resumes with the new prefactor; forking rules push their
alternate branch onto `worklist` and the loop continues on the contracted
branch; terminal rules push directly to `done`. When no rule matches the term
is pushed to `done` as-is.
"""
function _process_product!(
        c::CNum, ops::Vector{QSym},
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm}
    )
    while true
        length(ops) <= 1 && (push!(done, OrderedTerm(c, ops)); return)

        res = _try_step!(c, ops, worklist, done)
        if res isa CNum
            c = res
        elseif res isa Halt
            return
        else
            push!(done, OrderedTerm(c, ops))
            return
        end
    end
    return
end

# Try the leftmost rule matching an adjacent pair in `ops`.
@inline function _try_step!(
        c::CNum, ops::Vector{QSym},
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm}
    )::StepResult
    n = length(ops)
    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]
        res = _try_reduction!(a, b, i, c, ops, worklist, done)
        res === nothing || return res
        res = _try_swap!(a, b, i, c, ops, worklist)
        res === nothing || return res
    end
    return nothing
end

"""
    _apply_ordering(arg_c, ops) -> Vector{OrderedTerm}

Apply the eager normal-ordering rules to a product: algebraic reductions,
commutation swaps, and ground-state completeness. Used both inside `*` to
maintain canonical form and by [`normal_order`](@ref) as an idempotent
finalizer.

The implementation seeds a worklist with `(arg_c, copy(ops))`, drains it via
[`_process_product!`](@ref), then runs a final completeness pass over any
surviving standalone `σᵍᵍ` operands that the in-worklist GS expansion didn't
see as adjacent.
"""
function _apply_ordering(arg_c::CNum, ops::Vector{QSym})
    length(ops) <= 1 && return _expand_gs_oterms(OrderedTerm[OrderedTerm(arg_c, ops)])

    worklist = OrderedTerm[OrderedTerm(arg_c, copy(ops))]
    done = OrderedTerm[]
    while !isempty(worklist)
        t = pop!(worklist)
        _process_product!(t.prefactor, t.ops, worklist, done)
    end
    return _expand_gs_oterms(done)
end
