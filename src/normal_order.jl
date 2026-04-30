"""
    normal_order(expr::QField)
    normal_order(expr::QField, h::HilbertSpace)

Apply normal ordering to a quantum operator expression.

Places creation operators to the left of annihilation operators by applying
all commutation relations:
- **Fock**: ``[a, a^\\dagger] = 1``
- **Spin**: ``[S_j, S_k] = i\\epsilon_{jkl} S_l``
- **Phase space**: ``[p, x] = -i``

Also applies all algebraic identities (Transition composition, Pauli products),
simplifies prefactors, and collects like terms.

When a [`HilbertSpace`](@ref) `h` is provided, additionally rewrites ground-state
projectors of [`NLevelSpace`](@ref) subspaces via the completeness relation.

!!! note
    Under the default [`NormalOrder`](@ref) convention, operator products are already
    normal-ordered eagerly during `*`. This function is primarily useful after
    [`set_ordering!(LazyOrder())`](@ref set_ordering!) or for explicit re-ordering.

See also [`simplify`](@ref), [`normal_to_symmetric`](@ref), [`symmetric_to_normal`](@ref).
"""
function normal_order(op::QSym)
    return QAdd(QTermDict(QSym[op] => _CNUM_ONE), Index[], Tuple{Index, Index}[])
end

function normal_order(s::QAdd)
    d = QTermDict()
    for (ops, c) in s.arguments
        _iszero_cnum(c) && continue
        for (oc, oops) in _apply_ordering(c, ops, NormalOrder())
            _addto!(d, oops, oc)
        end
    end
    return QAdd(d, s.indices, s.non_equal)
end

function normal_order(op::QSym, h::HilbertSpace)
    return _apply_ground_state(normal_order(op), h)
end
function normal_order(s::QAdd, h::HilbertSpace)
    return _apply_ground_state(normal_order(s), h)
end

# ============================================================================
#  Normal ↔ Symmetric (Weyl) ordering conversion
# ============================================================================

"""
    normal_to_symmetric(expr) -> QAdd

Convert a normal-ordered expression to symmetric (Weyl) ordering.

The input **must** be in normal-ordered form (``a^\\dagger`` left of ``a`` per site).
Uses the combinatorial formula: for each normal-ordered term ``c \\cdot (a^\\dagger)^m a^n``,
the Weyl-ordered contribution is

``c \\sum_{k=0}^{\\min(m,n)} \\binom{m}{k} \\binom{n}{k} k! \\left(-\\tfrac{1}{2}\\right)^k (a^\\dagger)^{m-k} a^{n-k}``

For multi-site terms, the conversion applies independently per Fock site.
Non-Fock operators (Transition, Pauli, Spin, PhaseSpace) are left unchanged.

Exact inverse of [`symmetric_to_normal`](@ref).

See also [`symmetric_to_normal`](@ref), [`normal_order`](@ref).
"""
function normal_to_symmetric(s::QAdd)
    return _convert_ordering(s, -1 // 2)
end
function normal_to_symmetric(op::QSym)
    return normal_to_symmetric(QAdd(QTermDict(QSym[op] => _CNUM_ONE), Index[], Tuple{Index, Index}[]))
end

"""
    symmetric_to_normal(expr) -> QAdd

Convert a symmetric (Weyl) ordered expression to normal ordering.

The input **must** be in Weyl-ordered form. Uses the same combinatorial formula
as [`normal_to_symmetric`](@ref) but with ``+\\tfrac{1}{2}`` instead of ``-\\tfrac{1}{2}``.

Exact inverse: `symmetric_to_normal(normal_to_symmetric(x)) == x`.

See also [`normal_to_symmetric`](@ref), [`normal_order`](@ref).
"""
function symmetric_to_normal(s::QAdd)
    return _convert_ordering(s, 1 // 2)
end
function symmetric_to_normal(op::QSym)
    return symmetric_to_normal(QAdd(QTermDict(QSym[op] => _CNUM_ONE), Index[], Tuple{Index, Index}[]))
end

"""
    _convert_ordering(expr, α) -> QAdd

Apply the ordering conversion formula with parameter `α`:
- `α = -1/2`: normal → symmetric (Weyl)
- `α = +1/2`: symmetric (Weyl) → normal

For each term `c · ops`, groups Fock operators by site, counts `(m, n)` = (creates, destroys)
per site, and applies: `Σ_k C(m,k) C(n,k) k! α^k · (reduced ops)`.
"""
function _convert_ordering(s::QAdd, α::Rational)
    d = QTermDict()
    for (ops, c) in s.arguments
        _iszero_cnum(c) && continue
        _convert_term!(d, c, ops, α)
    end
    return QAdd(d, s.indices, s.non_equal)
end

function _convert_term!(d::QTermDict, c::CNum, ops::Vector{QSym}, α::Rational)
    # Find Fock operator groups by site key (space_index, index, name)
    sites = _fock_site_groups(ops)

    if isempty(sites)
        # No Fock operators — term unchanged
        _addto!(d, ops, c)
        return
    end

    # For each site, compute the conversion terms and combine via tensor product
    # Start with the single term [(c, ops)] and process each site
    current = _OTerm[(c, ops)]

    for (site_key, m, n) in sites
        next = _OTerm[]
        for (tc, tops) in current
            _apply_site_conversion!(next, tc, tops, site_key, m, n, α)
        end
        current = next
    end

    for (tc, tops) in current
        _addto!(d, tops, tc)
    end
    return
end

# Identify Fock (Create/Destroy) sites and count m (creates) and n (destroys) per site.
# Returns Vector of (site_key, m, n) tuples.
function _fock_site_groups(ops::Vector{QSym})
    sites = Tuple{Tuple{Int, Symbol, Symbol}, Int, Int}[]  # (key, m, n)
    counted = Dict{Tuple{Int, Symbol, Symbol}, Tuple{Int, Int}}()
    for op in ops
        if op isa Create || op isa Destroy
            key = (op.space_index, op.index.name, op.name)
            m_old, n_old = get(counted, key, (0, 0))
            if op isa Create
                counted[key] = (m_old + 1, n_old)
            else
                counted[key] = (m_old, n_old + 1)
            end
        end
    end
    for (key, (m, n)) in counted
        (m > 0 && n > 0) || continue  # only sites with both a† and a contribute
        push!(sites, (key, m, n))
    end
    return sites
end

# Apply the conversion formula for one Fock site to a set of terms.
# For site with m creates and n destroys:
# Σ_{k=1}^{min(m,n)} C(m,k) C(n,k) k! α^k · (remove k creates and k destroys from ops)
# plus the k=0 term (original, unchanged).
function _apply_site_conversion!(
        result::Vector{_OTerm}, c::CNum, ops::Vector{QSym},
        site_key::Tuple{Int, Symbol, Symbol}, m::Int, n::Int,
        α::Rational
    )
    # k=0: original term, unchanged
    push!(result, (c, copy(ops)))

    # k=1..min(m,n): remove k creates and k destroys, multiply by coefficient
    kmax = min(m, n)
    for k in 1:kmax
        coeff = _to_cnum(binomial(m, k) * binomial(n, k) * factorial(k) * α^k)
        new_ops = _remove_fock_ops(ops, site_key, k)
        push!(result, (_mul_cnum(c, coeff), new_ops))
    end
    return
end

# Remove k Create and k Destroy operators matching the given site key.
function _remove_fock_ops(
        ops::Vector{QSym}, site_key::Tuple{Int, Symbol, Symbol}, k::Int
    )
    new_ops = copy(ops)
    creates_removed = 0
    destroys_removed = 0

    # Remove Creates (scanning left to right)
    i = 1
    while i <= length(new_ops) && creates_removed < k
        op = new_ops[i]
        if op isa Create &&
                (op.space_index, op.index.name, op.name) == site_key
            deleteat!(new_ops, i)
            creates_removed += 1
        else
            i += 1
        end
    end

    # Remove Destroys (scanning left to right)
    i = 1
    while i <= length(new_ops) && destroys_removed < k
        op = new_ops[i]
        if op isa Destroy &&
                (op.space_index, op.index.name, op.name) == site_key
            deleteat!(new_ops, i)
            destroys_removed += 1
        else
            i += 1
        end
    end

    return new_ops
end

# Ground state rewriting: |g⟩⟨g| = 1 - Σ_{k≠g} |k⟩⟨k|

function _apply_ground_state(expr::QAdd, h::NLevelSpace)
    d = QTermDict()
    for (ops, c) in expr.arguments
        expanded = _expand_ground_state(c, ops, h)
        for (eops, ec) in expanded.arguments
            _addto!(d, eops, ec)
        end
    end
    return QAdd(d, expr.indices, expr.non_equal)
end

function _apply_ground_state(expr::QAdd, h::ProductSpace)
    result = expr
    for space in h.spaces
        result = _apply_ground_state(result, space)
    end
    return result
end

_apply_ground_state(expr::QAdd, ::FockSpace) = expr
_apply_ground_state(expr::QAdd, ::PauliSpace) = expr
_apply_ground_state(expr::QAdd, ::SpinSpace) = expr
_apply_ground_state(expr::QAdd, ::PhaseSpace) = expr

function _expand_ground_state(c::CNum, ops::Vector{QSym}, h::NLevelSpace)
    g = h.ground_state
    n = h.n
    for (idx, op) in enumerate(ops)
        if op isa Transition && op.i == g && op.j == g
            expanded = _OTerm[]
            # Identity term (remove the transition)
            id_ops = QSym[ops[j] for j in eachindex(ops) if j != idx]
            push!(expanded, (c, id_ops))
            # Subtraction terms
            for k in 1:n
                k == g && continue
                new_op = Transition(op.name, k, k, op.space_index, op.index)
                sub_ops = QSym[j == idx ? new_op : ops[j] for j in eachindex(ops)]
                push!(expanded, (-c, sub_ops))
            end
            # Recurse: expanded terms may still contain ground state projections
            d = QTermDict()
            for (tc, tops) in expanded
                sub = _expand_ground_state(tc, tops, h)
                for (sops, sc) in sub.arguments
                    _addto!(d, sops, sc)
                end
            end
            return QAdd(d, Index[], Tuple{Index, Index}[])
        end
    end
    d = QTermDict(ops => c)
    return QAdd(d, Index[], Tuple{Index, Index}[])
end
