# --- Global ordering convention ---

const _ORDERING_DEFAULT = Ref{OrderingConvention}(NormalOrder())
const ORDERING = ScopedValue{OrderingConvention}()

"""
    get_ordering() -> OrderingConvention

Return the active operator ordering convention. Inside a [`with_ordering`](@ref)
block the scoped value wins; otherwise returns the global default (initially
[`NormalOrder`](@ref), mutable via [`set_ordering!`](@ref)).
"""
function get_ordering()
    v = ScopedValues.get(ORDERING)
    return v === nothing ? _ORDERING_DEFAULT[] : something(v)
end

"""
    set_ordering!(ord::OrderingConvention)

Set the global default ordering convention. Affects subsequent `*` operations
unless overridden by an enclosing [`with_ordering`](@ref) scope.

Prefer [`with_ordering`](@ref) for transient changes — it is task-local and
cannot leak between tests.

See also [`get_ordering`](@ref), [`with_ordering`](@ref), [`NormalOrder`](@ref),
[`LazyOrder`](@ref).
"""
set_ordering!(ord::OrderingConvention) = (_ORDERING_DEFAULT[] = ord)

"""
    with_ordering(f, ord::OrderingConvention)

Run `f()` with `ord` as the active ordering convention. Restores on return,
even on exception. Task-local — concurrent tasks see their own value.

# Examples
```julia
with_ordering(LazyOrder()) do
    a * a'   # not reordered inside this block
end
```
"""
with_ordering(f, ord::OrderingConvention) = with(f, ORDERING => ord)

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
    _try_reduction!(a, b, i, c, ops, worklist, done, expand_gs) -> StepResult

Apply an ordering-independent algebraic identity to the pair `(a, b)` at
position `i`. Shared between the eager-ordering driver and the simplify pass.

When `expand_gs == true`, a Transition composition that would produce a ground-state
projector ``|g\\rangle\\langle g|`` is rewritten via the completeness relation
``|g\\rangle\\langle g| = 1 - \\sum_{k \\neq g} |k\\rangle\\langle k|`` directly,
keeping the [`NormalOrder`](@ref) result in canonical form. [`simplify`](@ref)
without a Hilbert space passes `false` to preserve the user-constructed shape;
[`simplify`](@ref)`(expr, h)` invokes the explicit cleanup pass.
"""
@inline function _try_reduction!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm},
        expand_gs::Bool
    )::StepResult
    # Transition: |i⟩⟨j| · |k⟩⟨l|
    if a isa Transition && b isa Transition && _same_site(a, b) && a.name == b.name
        if a.j != b.i
            push!(done, OrderedTerm(_to_cnum(0), QSym[]))
            return HALT
        end
        if expand_gs && a.i == a.ground_state && b.j == a.ground_state
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
    _try_swap!(a, b, i, c, ops, ord, worklist) -> StepResult

Apply an ordering-dependent commutation rule to the pair `(a, b)` at position
`i`. Forking rules push the swapped branch onto `worklist` and return the new
prefactor for in-place continuation on the contracted branch; `_try_swap!`
never returns `HALT`. The `LazyOrder` method is a no-op so the simplify pass
can share the same driver.
"""
@inline function _try_swap!(
        a::QSym, b::QSym, i::Int, c::CNum, ops::Vector{QSym},
        ::NormalOrder, worklist::Vector{OrderedTerm}
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

# LazyOrder skips swaps so the simplify pass can share `_process_product!`.
@inline _try_swap!(
    ::QSym, ::QSym, ::Int, ::CNum, ::Vector{QSym},
    ::LazyOrder, ::Vector{OrderedTerm}
)::StepResult = nothing

"""
    _process_product!(c, ops, ord, worklist, done, expand_gs)

Drive one term `(c, ops)` to completion. Single-output rules mutate `ops` in
place and the loop resumes with the new prefactor; forking rules push their
alternate branch onto `worklist` and the loop continues on the contracted
branch; terminal rules push directly to `done`. When no rule matches the term
is pushed to `done` as-is. Shared by [`_apply_ordering`](@ref) and
[`_apply_reductions`](@ref).
"""
function _process_product!(
        c::CNum, ops::Vector{QSym}, ord::OrderingConvention,
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm},
        expand_gs::Bool
    )
    while true
        length(ops) <= 1 && (push!(done, OrderedTerm(c, ops)); return)

        res = _try_step!(c, ops, ord, worklist, done, expand_gs)
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
        c::CNum, ops::Vector{QSym}, ord::OrderingConvention,
        worklist::Vector{OrderedTerm}, done::Vector{OrderedTerm},
        expand_gs::Bool
    )::StepResult
    n = length(ops)
    for i in 1:(n - 1)
        a, b = ops[i], ops[i + 1]
        res = _try_reduction!(a, b, i, c, ops, worklist, done, expand_gs)
        res === nothing || return res
        res = _try_swap!(a, b, i, c, ops, ord, worklist)
        res === nothing || return res
    end
    return nothing
end

"""
    _drive_worklist(arg_c, ops, ord, expand_gs) -> Vector{OrderedTerm}

Run the shared worklist algorithm: seed with `(arg_c, copy(ops))` and drain via
[`_process_product!`](@ref) under `ord`. Used by both [`_apply_ordering`](@ref)
(eager, `expand_gs = true`) and [`_apply_reductions`](@ref) (simplify pass,
`expand_gs = false` under `LazyOrder`).
"""
function _drive_worklist(
        arg_c::CNum, ops::Vector{QSym}, ord::OrderingConvention, expand_gs::Bool
    )
    length(ops) <= 1 && return OrderedTerm[OrderedTerm(arg_c, ops)]

    worklist = OrderedTerm[OrderedTerm(arg_c, copy(ops))]
    done = OrderedTerm[]
    while !isempty(worklist)
        t = pop!(worklist)
        _process_product!(t.prefactor, t.ops, ord, worklist, done, expand_gs)
    end
    return done
end

"""
    _apply_ordering(arg_c, ops, ord) -> Vector{OrderedTerm}

Apply the eager ordering rules of `ord` to a product. Under [`LazyOrder`](@ref)
the eager arithmetic is the identity — `simplify` and `normal_order` apply
rules explicitly via [`_apply_reductions`](@ref) and `_apply_ground_state`.
"""
_apply_ordering(arg_c::CNum, ops::Vector{QSym}, ord::OrderingConvention) =
    _drive_worklist(arg_c, ops, ord, true)

_apply_ordering(arg_c::CNum, ops::Vector{QSym}, ::LazyOrder) =
    OrderedTerm[OrderedTerm(arg_c, ops)]
