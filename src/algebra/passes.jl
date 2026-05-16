# ============================================================================
#  Passes: the low-level building blocks of the canonicalization pipeline.
#
#  Each pass takes a `sink::F` first (so `do` blocks work), an `ops::Vector{QSym}`
#  and a `c::CNum`. Forking passes may invoke the sink multiple times; folding
#  passes invoke it once. Aliasing contracts are documented per pass.
# ============================================================================

# ---------------------------------------------------------------------------
#  _partial_sort!: distinct-site reorder
# ---------------------------------------------------------------------------

"""
    _partial_sort!(ops::Vector{QSym}, ne::Vector{NonEqualPair})

In-place stable partial-sort using `_site_compare`. Distinct-site adjacents are
placed in canonical order; `Equal` and `Undetermined` pairs are left in
physical order so the sibling passes can interpret them.
"""
function _partial_sort!(ops::Vector{QSym}, ne::Vector{NonEqualPair})
    n = length(ops)
    n < 2 && return ops
    # Insertion sort: stable, O(n²) worst case, n is small in practice.
    for i in 2:n
        j = i
        while j > 1
            cmp = _site_compare(ops[j - 1], ops[j], ne)
            cmp === Greater || break
            ops[j - 1], ops[j] = ops[j], ops[j - 1]
            j -= 1
        end
    end
    return ops
end

# ---------------------------------------------------------------------------
#  _canonicalize_to_dict!: terminal sink
# ---------------------------------------------------------------------------

"""
    _canonicalize_to_dict!(out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair})

Terminal sink for the pipeline. Constructs `QTerm(ops, ne)`, looks up in `out`,
sums `c`, drops zero-coefficient entries. Takes ownership of `ops`; the caller
must not mutate after this call. Pre-condition: `ops` is in canonical form.
"""
@inline function _canonicalize_to_dict!(
        out::QTermDict, ops::Vector{QSym}, c::CNum, ne::Vector{NonEqualPair},
    )
    _iszero_cnum(c) && return out
    return _addto_key!(out, QTerm(ops, _canonical_ne(ne)), c)
end

# ---------------------------------------------------------------------------
#  _reduce_ops: fold same-site composition
# ---------------------------------------------------------------------------

"""
    _reduce_ops(sink, ops, c)

Fold adjacent provably-same-site pairs via `_reduce_pair`. Single output per
input. Mutates `ops` in place; sink receives the same Vector.

Pre: `ops` is in canonical partial-sort order. Post: no adjacent provably-
same-site pair `(a, b)` for which `_reduce_pair(a, b)` returns non-`nothing`.
"""
@inline function _reduce_ops(sink::F, ops::Vector{QSym}, c::CNum) where {F}
    n = length(ops)
    n < 2 && (sink(ops, c); return)
    i = 1
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        res = _reduce_pair(a, b)
        if res === nothing
            i += 1
        elseif res isa CNum
            # Pair contracts to a scalar (e.g. σⁱʲσʲᵏ when j≠k yielding 0, or Pauli δ).
            _iszero_cnum(res) && return
            c = _mul_cnum(c, res)
            deleteat!(ops, i:(i + 1))
            i = max(i - 1, 1)
        else
            new_op, factor = res
            c = _mul_cnum(c, factor)
            ops[i] = new_op
            deleteat!(ops, i + 1)
            i = max(i - 1, 1)
        end
    end
    sink(ops, c)
    return
end

# ---------------------------------------------------------------------------
#  _commute_ops: same-site commutation, forks per swap
# ---------------------------------------------------------------------------

"""
    _commute_ops(sink, ops, c)

Apply commutation swaps for adjacent same-site pairs whose `_can_commute` is
false. May fork: emits the swapped branch through the sink, plus a residual
branch where the pair is replaced by `residual_ops` (empty for identity-residual
operators like Fock/PhaseSpace; one operator for Spin commutators).

Pre-condition: pairs that reduce locally (Transition, Pauli) should have been
folded by `_reduce_ops` first. The residual branch mutates `ops` in place; the
swap branch receives a `copy(ops)`.
"""
@inline function _commute_ops(sink::F, ops::Vector{QSym}, c::CNum) where {F}
    _commute_ops_at(sink, ops, c, 1)
    return
end

function _commute_ops_at(sink::F, ops::Vector{QSym}, c::CNum, start::Int) where {F}
    i = start
    while i < length(ops)
        a, b = ops[i], ops[i + 1]
        cmp = _site_compare(a, b, _EMPTY_NE)
        if cmp === Equal && !_can_commute(a, b)
            sw_b, sw_a, residual_coeff, residual_ops = _commute_pair(a, b)
            # Swap branch: stepping back to i-1 because ops[i-1] may now form
            # a new non-canonical pair with sw_b.
            swapped = copy(ops)
            swapped[i] = sw_b
            swapped[i + 1] = sw_a
            _commute_ops_at(sink, swapped, c, max(i - 1, 1))
            # Residual branch: replace pair with residual_ops (empty or one op).
            new_c = _mul_cnum(c, residual_coeff)
            _iszero_cnum(new_c) && return
            deleteat!(ops, i:(i + 1))
            for (k, op) in enumerate(residual_ops)
                insert!(ops, i + k - 1, op)
            end
            c = new_c
            i = max(i - 1, 1)
        else
            i += 1
        end
    end
    sink(ops, c)
    return
end

# ---------------------------------------------------------------------------
#  _expand_gs_ops: ground-state completeness rewriting (opt-in)
# ---------------------------------------------------------------------------

"""
    _expand_gs_ops(sink, ops, c)

Apply `σᵍᵍ → 1 - Σ_{k≠g} σᵏᵏ` to every `σᵍᵍ` in `ops`. May fork by `n_levels`
per `σᵍᵍ`; recurses to handle multiple `σᵍᵍ` atoms in one term. Each branch
receives its own `copy(ops)`; the original is not mutated.
"""
function _expand_gs_ops(sink::F, ops::Vector{QSym}, c::CNum) where {F}
    idx = 0
    for k in eachindex(ops)
        _ground_state_expand(ops[k]) === nothing || (idx = k; break)
    end
    idx == 0 && (sink(ops, c); return)

    op = ops[idx]
    (g, n_levels, _site) = _ground_state_expand(op)
    id_ops = QSym[ops[k] for k in eachindex(ops) if k != idx]
    _expand_gs_ops(sink, id_ops, c)
    neg_c = _neg_cnum(c)
    for k in 1:n_levels
        k == g && continue
        new_ops = copy(ops)
        new_ops[idx] = Transition(op.name, k, k, op.space_index, op.index, g, n_levels)
        _expand_gs_ops(sink, new_ops, neg_c)
    end
    return
end

# ---------------------------------------------------------------------------
#  _substitute_ops: per-term substitution, may fork on QAdd values
# ---------------------------------------------------------------------------

"""
    _substitute_ops(sink, ops, c, d)

Walk `ops` applying substitutions from `d`. Supported value types:

- `QSym`            — replaced in place
- `Number` / `CNum` — folded into the coefficient (operator removed)
- `QAdd`            — spliced; forks once per term of the `QAdd`
- symbolic key      — passes through to `_substitute_cnum` on the coefficient
"""
function _substitute_ops(sink::F, ops::Vector{QSym}, c::CNum, d) where {F}
    for i in eachindex(ops)
        op = ops[i]
        haskey(d, op) || continue
        val = d[op]
        if val isa QSym
            new_ops = copy(ops)
            new_ops[i] = val
            return _substitute_ops(sink, new_ops, c, d)
        elseif val isa Number
            new_c = _mul_cnum(c, _to_cnum(val))
            _iszero_cnum(new_c) && return
            new_ops = QSym[ops[k] for k in eachindex(ops) if k != i]
            return _substitute_ops(sink, new_ops, new_c, d)
        elseif val isa QAdd
            for (vt, vc) in val
                spliced_ops = copy(ops)
                deleteat!(spliced_ops, i)
                for (k, vop) in enumerate(vt.ops)
                    insert!(spliced_ops, i + k - 1, vop)
                end
                new_c = _mul_cnum(c, vc)
                _substitute_ops(sink, spliced_ops, new_c, d)
            end
            return
        end
    end
    new_c = _substitute_cnum(c, d)
    _iszero_cnum(new_c) && return
    sink(ops, new_c)
    return
end

function _substitute_cnum(c::CNum, d)
    rep, imp = real(c), imag(c)
    new_rep = SymbolicUtils.substitute(Symbolics.value(rep), d)
    new_imp = SymbolicUtils.substitute(Symbolics.value(imp), d)
    return Num(new_rep) + im * Num(new_imp)
end
