# ============================================================================
#  Normal <-> Symmetric (Weyl) ordering conversion
# ============================================================================

# Local data structure for the Weyl-conversion helpers.
struct OrderedTerm
    prefactor::CNum
    ops::Vector{QSym}
end

"""
    normal_to_symmetric(expr) -> QAdd

Convert a normal-ordered expression to symmetric (Weyl) ordering.

Exact inverse of [`symmetric_to_normal`](@ref).
"""
function normal_to_symmetric(s::QAdd)
    return _convert_ordering(s, -1 // 2)
end
normal_to_symmetric(op::QSym) = normal_to_symmetric(_single_qadd(_CNUM_ONE, QSym[op]))

"""
    symmetric_to_normal(expr) -> QAdd

Convert a symmetric (Weyl) ordered expression to normal ordering.

Exact inverse: `symmetric_to_normal(normal_to_symmetric(x)) == x`.
"""
function symmetric_to_normal(s::QAdd)
    return _convert_ordering(s, 1 // 2)
end
symmetric_to_normal(op::QSym) = symmetric_to_normal(_single_qadd(_CNUM_ONE, QSym[op]))

function _convert_ordering(s::QAdd, α::Rational)
    d = QTermDict()
    for (term, c) in s.arguments
        _iszero_cnum(c) && continue
        _convert_term!(d, c, term.ops, term.ne, α)
    end
    return QAdd(d, copy(s.indices))
end

function _convert_term!(
        d::QTermDict, c::CNum, ops::Vector{QSym},
        ne::Vector{NonEqualPair}, α::Rational
    )
    sites = _fock_site_groups(ops)
    if isempty(sites)
        _addto!(d, ops, c, ne)
        return
    end

    current = OrderedTerm[OrderedTerm(c, ops)]
    for (site_key, m, n) in sites
        next = OrderedTerm[]
        for t in current
            _apply_site_conversion!(next, t.prefactor, t.ops, site_key, m, n, α)
        end
        current = next
    end

    for t in current
        _addto!(d, t.ops, t.prefactor, ne)
    end
    return
end

function _fock_site_groups(ops::Vector{QSym})
    sites = Tuple{Tuple{Int, Symbol, Symbol}, Int, Int}[]
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
        (m > 0 && n > 0) || continue
        push!(sites, (key, m, n))
    end
    return sites
end

function _apply_site_conversion!(
        result::Vector{OrderedTerm}, c::CNum, ops::Vector{QSym},
        site_key::Tuple{Int, Symbol, Symbol}, m::Int, n::Int,
        α::Rational
    )
    push!(result, OrderedTerm(c, copy(ops)))
    kmax = min(m, n)
    for k in 1:kmax
        coeff = _to_cnum(binomial(m, k) * binomial(n, k) * factorial(k) * α^k)
        new_ops = _remove_fock_ops(ops, site_key, k)
        push!(result, OrderedTerm(_mul_cnum(c, coeff), new_ops))
    end
    return
end

function _remove_fock_ops(
        ops::Vector{QSym}, site_key::Tuple{Int, Symbol, Symbol}, k::Int
    )
    new_ops = copy(ops)
    creates_removed = 0
    destroys_removed = 0
    i = 1
    while i <= length(new_ops) && creates_removed < k
        op = new_ops[i]
        if op isa Create && (op.space_index, op.index.name, op.name) == site_key
            deleteat!(new_ops, i)
            creates_removed += 1
        else
            i += 1
        end
    end
    i = 1
    while i <= length(new_ops) && destroys_removed < k
        op = new_ops[i]
        if op isa Destroy && (op.space_index, op.index.name, op.name) == site_key
            deleteat!(new_ops, i)
            destroys_removed += 1
        else
            i += 1
        end
    end
    return new_ops
end
