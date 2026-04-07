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

Two operators are on the same site iff they share space_index, copy_index, AND index.
Only operators on the same site can interact (apply commutation/composition rules).
"""
function _same_site(a::QSym, b::QSym)
    return a.space_index == b.space_index &&
        a.copy_index == b.copy_index &&
        a.index == b.index
end

# Levi-Civita lookup: _levi_civita[j][k] = ε_{jk(6-j-k)} for j,k ∈ {1,2,3}
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

# Worklist term: (prefactor, operators). No struct — just a tuple.
const _OTerm = Tuple{CNum, Vector{QSym}}

"""
    _try_algebraic_reduction!(a, b, i, c, ops, worklist, done) -> Bool

Apply ordering-independent algebraic identities for a single adjacent pair.
Returns `true` if a rule fired (and pushed results to worklist/done).
Shared by both the ordering worklist and the simplification worklist.
"""
function _try_algebraic_reduction!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        worklist::Vector{_OTerm}, done::Vector{_OTerm}
    )
    # Transition: |i⟩⟨j| · |k⟩⟨l|
    if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
        if a.j == b.i
            ops[i] = Transition(a.name, a.i, b.j, a.space_index, a.copy_index, a.index)
            deleteat!(ops, i + 1)
            push!(worklist, (c, ops))
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
            ops[i] = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index, a.index)
            deleteat!(ops, i + 1)
            push!(worklist, (_mul_cnum(c, _to_cnum(im * eps)), ops))
        end
        return true
    end

    return false
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

        # Ordering-independent reductions (Transition, Pauli)
        _try_algebraic_reduction!(a, b, i, c, ops, worklist, done) && return

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
        ops[i] = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.copy_index, a.index)
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
