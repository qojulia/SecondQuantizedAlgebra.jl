# --- Global ordering convention ---

const ORDERING = Ref{OrderingConvention}(NormalOrder())

"""
    set_ordering!(ord::OrderingConvention)

Set the global operator ordering convention. All subsequent `*` operations on
[`QSym`](@ref) operators will apply this convention eagerly.

# Examples
```julia
set_ordering!(LazyOrder())       # disable eager ordering
set_ordering!(NormalOrder())     # restore default
```

See also [`get_ordering`](@ref), [`NormalOrder`](@ref), [`LazyOrder`](@ref).
"""
set_ordering!(ord::OrderingConvention) = (ORDERING[] = ord)

"""
    get_ordering() -> OrderingConvention

Return the current global ordering convention (default: [`NormalOrder`](@ref)).

See also [`set_ordering!`](@ref).
"""
get_ordering() = ORDERING[]

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

# Worklist term: (prefactor, operators). No struct — just a tuple.
const _OTerm = Tuple{CNum, Vector{QSym}}

"""
    _try_algebraic_reduction!(a, b, i, c, ops, worklist, done, expand_gs) -> Bool

Apply ordering-independent algebraic identities for a single adjacent pair.
Returns `true` if a rule fired (and pushed results to worklist/done).
Shared by both the ordering worklist and the simplification worklist.

When `expand_gs == true`, a Transition composition that would produce a ground-state
projector ``|g\\rangle\\langle g|`` is rewritten via the completeness relation
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g} |k\\rangle\\langle k|`` directly,
keeping the [`NormalOrder`](@ref) result in canonical form. Pass `false` from
[`simplify`](@ref) (no Hilbert-space argument) so the user-constructed shape is
preserved; [`simplify`](@ref)`(expr, h)` invokes the explicit cleanup pass.
"""
function _try_algebraic_reduction!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        worklist::Vector{_OTerm}, done::Vector{_OTerm},
        expand_gs::Bool
    )
    # Transition: |i⟩⟨j| · |k⟩⟨l|
    if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
        if a.j == b.i
            if expand_gs && a.i == a.ground_state && b.j == a.ground_state
                _push_ground_state_expansion!(c, ops, i, a, worklist)
            else
                ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.index, a.ground_state, a.n_levels)
                deleteat!(ops, i + 1)
                push!(worklist, (c, ops))
            end
        else
            push!(done, (_to_cnum(0), QSym[]))
        end
        return true
    end

    # Pauli: σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
    if a isa Pauli && b isa Pauli && _same_site(a, b) && a.name == b.name
        if a.axis == b.axis
            deleteat!(ops, i:(i + 1))
            push!(worklist, (c, ops))
        else
            eps = _levi_civita[a.axis][b.axis]
            ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
            deleteat!(ops, i + 1)
            push!(worklist, (_mul_cnum(c, _to_cnum(im * eps)), ops))
        end
        return true
    end

    return false
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
        c::CNum, ops::Vector{QSym}, i::Int, a::Transition, worklist::Vector{_OTerm}
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
    push!(worklist, (c, id_ops))
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
        push!(worklist, (neg_c, new_ops))
    end
    return nothing
end

"""
    _apply_ordering(arg_c, ops, ordering) -> Vector{_OTerm}

Apply ordering rules to a product of operators using a worklist algorithm.
Returns a vector of fully-ordered (prefactor, ops) tuples.
"""
function _apply_ordering(arg_c::CNum, ops::Vector{QSym}, ord::OrderingConvention)
    isempty(ops) && return _OTerm[(arg_c, ops)]
    length(ops) == 1 && return _OTerm[(arg_c, ops)]

    worklist = _OTerm[(arg_c, copy(ops))]
    done = _OTerm[]

    while !isempty(worklist)
        c, ops_w = pop!(worklist)
        _order_product!(c, ops_w, ord, worklist, done)
    end

    return done
end

"""
    _order_product!(c, ops, ord, worklist, done)

Process one term. Scans adjacent operator pairs for the first applicable
rule, pushes resulting terms to `worklist`, or pushes to `done` if fully ordered.
"""
function _order_product!(
        c::CNum, ops::Vector{QSym}, ord::OrderingConvention,
        worklist::Vector{_OTerm}, done::Vector{_OTerm}
    )
    n = length(ops)

    if n <= 1
        push!(done, (c, ops))
        return
    end

    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]

        # Ordering-independent reductions (Transition, Pauli).
        # NormalOrder is canonical: eagerly expand σᵍᵍ via completeness when
        # produced by Transition composition. LazyOrder skips this entire
        # function via the identity overload of `_apply_ordering` below.
        _try_algebraic_reduction!(a, b, i, c, ops, worklist, done, true) && return

        # Ordering-dependent swaps
        _apply_ordering_swap!(a, b, i, c, ops, ord, worklist) && return
    end

    return push!(done, (c, ops))
end

# NormalOrder swaps
function _apply_ordering_swap!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        ::NormalOrder, worklist::Vector{_OTerm}
    )
    # Fock: a*a' -> a'*a + 1
    if a isa Destroy && b isa Create && _same_site(a, b) && a.name == b.name
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, (c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, (c, ops))
        return true
    end

    # Spin: [Sj, Sk] = i eps_jkl Sl (swap out-of-order axes)
    if a isa Spin && b isa Spin && _same_site(a, b) && a.name == b.name && a.axis > b.axis
        eps = _levi_civita[a.axis][b.axis]
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, (c, swapped))
        ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
        deleteat!(ops, i + 1)
        push!(worklist, (_mul_cnum(c, _to_cnum(im * eps)), ops))
        return true
    end

    # PhaseSpace: P*X -> X*P - i
    if a isa Momentum && b isa Position && _same_site(a, b)
        swapped = copy(ops)
        swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
        push!(worklist, (c, swapped))
        deleteat!(ops, i:(i + 1))
        push!(worklist, (_mul_cnum(_CNUM_NEG_IM, c), ops))
        return true
    end

    return false
end

# LazyOrder: skip all ordering (identity transform)
function _apply_ordering(arg_c::CNum, ops::Vector{QSym}, ::LazyOrder)
    return _OTerm[(arg_c, ops)]
end
