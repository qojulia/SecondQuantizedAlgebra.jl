const QTermDict = Dict{Vector{QSym}, CNum}

function _addto!(d::QTermDict, key::Vector{QSym}, val::CNum)
    existing = get(d, key, nothing)
    if existing === nothing
        _iszero_cnum(val) || (d[key] = val)
    else
        new_val = _add_cnum(existing, val)
        if _iszero_cnum(new_val)
            delete!(d, key)
        else
            d[key] = new_val
        end
    end
    return d
end

function _drop_zeros!(d::QTermDict)
    for (ops, c) in d
        _iszero_cnum(c) && delete!(d, ops)
    end
    return d
end

function _merge_unique(a::Vector{T}, b::Vector{T}) where {T}
    isempty(a) && return b
    isempty(b) && return a
    result = copy(a)
    for x in b
        x ∉ result && push!(result, x)
    end
    return result
end

"""
    QAdd <: QTerm

The sole compound expression type — a sum of eagerly-ordered operator products.

All arithmetic on [`QSym`](@ref) operators returns `QAdd`. Internally stores a
`Dict{Vector{QSym}, CNum}` mapping operator sequences to their `Complex{Num}` prefactors.
Like terms are automatically collected on construction and zero-prefactor terms are dropped.

# Fields
- `arguments::QTermDict` — operator sequence → prefactor mapping
- `indices::Vector{Index}` — summation indices (empty for a regular sum)
- `non_equal::Vector{Tuple{Index,Index}}` — pairwise inequality constraints on indices

# Iteration
Iterating over a `QAdd` yields `Pair{Vector{QSym}, CNum}` entries from the internal dict.
For deterministic ordering (printing, comparison), use [`sorted_arguments`](@ref).

See also [`prefactor`](@ref), [`operators`](@ref), [`Σ`](@ref).
"""
struct QAdd <: QTerm
    arguments::QTermDict
    indices::Vector{Index}
    non_equal::Vector{Tuple{Index, Index}}
end

# Convenience: single-term QAdd
function _single_qadd(c::CNum, ops::Vector{QSym})
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
    return QAdd(QTermDict(ops => c), Index[], Tuple{Index, Index}[])
end

# Internal: build QAdd from ordering worklist output, applying eager ordering
function _qadd_from_oterms(
        terms::Vector{_OTerm}, indices::Vector{Index},
        non_equal::Vector{Tuple{Index, Index}}
    )
    d = QTermDict()
    for (c, ops) in terms
        for (oc, oops) in _apply_ordering(c, ops, ORDERING[])
            _addto!(d, oops, oc)
        end
    end
    return QAdd(d, indices, non_equal)
end

Base.length(a::QAdd) = length(a.arguments)
Base.iszero(a::QAdd) = isempty(a.arguments)

function Base.isequal(a::QAdd, b::QAdd)
    isequal(a.arguments, b.arguments) || return false
    a.indices == b.indices || return false
    a.non_equal == b.non_equal || return false
    return true
end
Base.:(==)(a::QAdd, b::QAdd) = isequal(a, b)
function Base.hash(q::QAdd, h::UInt)
    return hash(:QAdd, hash(q.arguments, hash(q.indices, hash(q.non_equal, h))))
end

function Base.adjoint(q::QAdd)
    d = QTermDict()
    for (ops, c) in q.arguments
        adj_ops = QSym[adjoint(op) for op in reverse(ops)]
        _site_sort!(adj_ops)
        _addto!(d, adj_ops, conj(c))
    end
    return QAdd(d, q.indices, q.non_equal)
end

# --- Iteration: yields Pair{Vector{QSym}, CNum} from Dict directly ---

Base.iterate(q::QAdd) = iterate(q.arguments)
Base.iterate(q::QAdd, state) = iterate(q.arguments, state)
Base.eltype(::Type{QAdd}) = Pair{Vector{QSym}, CNum}

# --- Sorted term access for printing and TermInterface ---

"""
    sorted_arguments(q::QAdd) -> Vector{QAdd}

Return each term of `q` as a single-entry [`QAdd`](@ref), in deterministic sort order.

Sort key: `(term length, operator keys...)`. Used internally by `show` and
TermInterface for reproducible output.
"""
function sorted_arguments(q::QAdd)
    isempty(q.arguments) && return QAdd[]
    pairs = sort!(collect(q.arguments); by = p -> _term_sort_key(p.first))
    return QAdd[_single_qadd(c, ops) for (ops, c) in pairs]
end

_term_sort_key(ops::Vector{QSym}) = (length(ops), map(_full_op_key, ops)...)
_full_op_key(op::QSym) = (_sort_key(op)..., _type_order(op), op.name)
_type_order(::Destroy) = 0
_type_order(::Create) = 1
_type_order(::Transition) = 2
_type_order(::Pauli) = 3
_type_order(::Spin) = 4
_type_order(::Position) = 5
_type_order(::Momentum) = 6

"""
    Base.getindex(q::QAdd, key::Vector{QSym}) -> CNum

Look up the prefactor for a given operator sequence. Returns zero if absent.
"""
Base.getindex(q::QAdd, key::Vector{QSym}) = get(q.arguments, key, _CNUM_ZERO)

# --- QAdd accessor helpers ---

"""
    prefactor(s::QAdd) -> CNum

Return the `Complex{Num}` prefactor of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term expressions,
iterate over the `QAdd` directly.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
prefactor(2 * a' * a)   # 2 + 0im
```
"""
function prefactor(s::QAdd)
    length(s.arguments) == 1 || throw(ArgumentError("prefactor requires a single-term expression, got $(length(s.arguments)) terms"))
    return first(values(s.arguments))
end

"""
    operators(s::QAdd) -> Vector{QSym}

Return the ordered operator sequence of a single-term [`QAdd`](@ref).

Throws `ArgumentError` if `s` contains more than one term. For multi-term expressions,
iterate over the `QAdd` directly.

# Examples
```julia
h = FockSpace(:f)
@qnumbers a::Destroy(h)
operators(a' * a)   # [Create(:a, 1, 1, NO_INDEX), Destroy(:a, 1, 1, NO_INDEX)]
```
"""
function operators(s::QAdd)
    length(s.arguments) == 1 || throw(ArgumentError("operators requires a single-term expression, got $(length(s.arguments)) terms"))
    return first(keys(s.arguments))
end

# --- Eagerly ordered product helper ---

function _ordered_qadd(c::CNum, ops::Vector{QSym})
    _iszero_cnum(c) && return QAdd(QTermDict(), Index[], Tuple{Index, Index}[])
    ordered = _apply_ordering(c, ops, ORDERING[])
    d = QTermDict()
    for (oc, oops) in ordered
        _addto!(d, oops, oc)
    end
    return QAdd(d, Index[], Tuple{Index, Index}[])
end

# ============================================================================
#  Multiplication — always returns QAdd, eagerly ordered
# ============================================================================

function Base.:*(a::QSym, b::QSym)
    ops = QSym[a, b]
    _site_sort!(ops)
    return _ordered_qadd(_CNUM_ONE, ops)
end

Base.:*(a::QSym, b::Number) = _single_qadd(_to_cnum(b), QSym[a])
Base.:*(b::Number, a::QSym) = a * b

function Base.:*(a::QAdd, b::Number)
    cb = _to_cnum(b)
    d = QTermDict(ops => _mul_cnum(c, cb) for (ops, c) in a.arguments)
    _drop_zeros!(d)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:*(a::Number, b::QAdd) = b * a

function Base.:*(a::QAdd, b::QSym)
    d = QTermDict()
    for (ops, c) in a.arguments
        n = length(ops)
        new_ops = Vector{QSym}(undef, n + 1)
        copyto!(new_ops, 1, ops, 1, n)
        new_ops[n + 1] = b
        _site_sort!(new_ops)
        for (oc, oops) in _apply_ordering(c, new_ops, ORDERING[])
            _addto!(d, oops, oc)
        end
    end
    product = QAdd(d, copy(a.indices), copy(a.non_equal))
    if has_index(b.index) && !isempty(a.indices)
        return _apply_diagonal_split(product, b.index)
    end
    return product
end
function Base.:*(a::QSym, b::QAdd)
    d = QTermDict()
    for (ops, c) in b.arguments
        n = length(ops)
        new_ops = Vector{QSym}(undef, n + 1)
        new_ops[1] = a
        copyto!(new_ops, 2, ops, 1, n)
        _site_sort!(new_ops)
        for (oc, oops) in _apply_ordering(c, new_ops, ORDERING[])
            _addto!(d, oops, oc)
        end
    end
    product = QAdd(d, copy(b.indices), copy(b.non_equal))
    if has_index(a.index) && !isempty(b.indices)
        return _apply_diagonal_split(product, a.index)
    end
    return product
end

function Base.:*(a::QAdd, b::QAdd)
    if !isempty(a.indices) && !isempty(b.indices)
        for idx in a.indices
            if idx in b.indices
                throw(
                    ArgumentError(
                        "Summation index $(idx.name) appears in both factors. " *
                            "Use `change_index` to re-index one side, or use different indices."
                    )
                )
            end
        end
    end
    d = QTermDict()
    for (ops_a, c_a) in a.arguments, (ops_b, c_b) in b.arguments
        new_ops = vcat(ops_a, ops_b)
        _site_sort!(new_ops)
        for (oc, oops) in _apply_ordering(_mul_cnum(c_a, c_b), new_ops, ORDERING[])
            _addto!(d, oops, oc)
        end
    end
    indices = _merge_unique(a.indices, b.indices)
    non_equal = _merge_unique(a.non_equal, b.non_equal)
    product = QAdd(d, indices, non_equal)
    if !isempty(a.indices) || !isempty(b.indices)
        if !isempty(a.indices)
            for ext_idx in _get_indices_from_terms(b)
                has_index(ext_idx) && (product = _apply_diagonal_split(product, ext_idx))
            end
        end
        if !isempty(b.indices)
            for ext_idx in _get_indices_from_terms(a)
                has_index(ext_idx) && (product = _apply_diagonal_split(product, ext_idx))
            end
        end
        if !isempty(a.indices) && !isempty(b.indices)
            for idx_a in a.indices, idx_b in b.indices
                if idx_a.space_index == idx_b.space_index
                    product = _apply_diagonal_split(product, idx_a)
                    product = _apply_diagonal_split(product, idx_b)
                end
            end
        end
    end
    return product
end

# ============================================================================
#  Addition — always returns QAdd
# ============================================================================

function Base.:+(a::QSym, b::QSym)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, Index[], Tuple{Index, Index}[])
end

function Base.:+(a::QAdd, b::QSym)
    d = copy(a.arguments)
    _addto!(d, QSym[b], _CNUM_ONE)
    return QAdd(d, a.indices, a.non_equal)
end
Base.:+(a::QSym, b::QAdd) = b + a

function Base.:+(a::QAdd, b::QAdd)
    d = copy(a.arguments)
    for (ops, c) in b.arguments
        _addto!(d, ops, c)
    end
    indices = _merge_unique(a.indices, b.indices)
    non_equal = _merge_unique(a.non_equal, b.non_equal)
    return QAdd(d, indices, non_equal)
end

function Base.:+(a::QSym, b::Number)
    d = QTermDict()
    _addto!(d, QSym[a], _CNUM_ONE)
    _addto!(d, QSym[], _to_cnum(b))
    return QAdd(d, Index[], Tuple{Index, Index}[])
end
Base.:+(a::Number, b::QSym) = b + a
function Base.:+(a::QAdd, b::Number)
    d = copy(a.arguments)
    _addto!(d, QSym[], _to_cnum(b))
    return QAdd(d, a.indices, a.non_equal)
end
Base.:+(a::Number, b::QAdd) = b + a

Base.zero(::Type{QAdd}) = _zero_qadd()
Base.zero(::QAdd) = _zero_qadd()

# ============================================================================
#  Subtraction, negation, division, power
# ============================================================================

Base.:-(a::QSym) = _single_qadd(_CNUM_NEG1, QSym[a])
Base.:-(a::QAdd) = QAdd(QTermDict(ops => _neg_cnum(c) for (ops, c) in a.arguments), a.indices, a.non_equal)
Base.:-(a::QField, b::QField) = a + (-b)
Base.:-(a::QField, b::Number) = a + (-b)
Base.:-(a::Number, b::QField) = a + (-b)

Base.:/(a::QSym, b::Number) = a * inv(b)
Base.:/(a::QAdd, b::Number) = a * inv(b)

function Base.:^(a::QSym, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, QSym[])
    ops = QSym[a for _ in 1:n]
    return _ordered_qadd(_CNUM_ONE, ops)
end

function Base.:^(a::QAdd, n::Integer)
    n >= 0 || throw(ArgumentError("Negative powers not supported"))
    n == 0 && return _single_qadd(_CNUM_ONE, QSym[])
    result = a
    for _ in 2:n
        result = result * a
    end
    return result
end

# ============================================================================
#  Index helpers
# ============================================================================

function _depends_on_index_term(c::CNum, ops::Vector{QSym}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    isym = SymbolicUtils.unwrap(idx.sym)
    for part in (real(c), imag(c))
        vars = Symbolics.get_variables(part)
        any(v -> isequal(v, isym), vars) && return true
    end
    return false
end

function _any_depends_on_index(s::QAdd, idx::Index)
    for (ops, c) in s.arguments
        _depends_on_index_term(c, ops, idx) && return true
    end
    return false
end

"""
    _diagonal_split(qadd, sum_idx, non_equal)

Extract diagonal terms for summation index `sum_idx`. For each term containing a free
index `j` on the same space as `sum_idx` (and not already constrained `sum_idx != j`),
create the diagonal term via `change_index(term, sum_idx, j)` and record the constraint.

Returns `(diag_terms::Vector{_OTerm}, new_non_equal::Vector{Tuple{Index,Index}})`.
"""
function _diagonal_split(
        qadd::QAdd, sum_idx::Index,
        non_equal::Vector{Tuple{Index, Index}}
    )
    diag_terms = _OTerm[]
    new_ne = copy(non_equal)
    seen = Index[]
    for (ops, c) in qadd.arguments
        empty!(seen)
        for op in ops
            idx = op.index
            has_index(idx) || continue
            idx ∈ seen && continue
            push!(seen, idx)
            if idx != sum_idx &&
                    idx.space_index == sum_idx.space_index &&
                    !((sum_idx, idx) in new_ne) &&
                    !((idx, sum_idx) in new_ne)
                new_ops = QSym[change_index(op, sum_idx, idx) for op in ops]
                new_c = change_index(c, sum_idx, idx)
                push!(diag_terms, (new_c, new_ops))
                push!(new_ne, (sum_idx, idx))
            end
        end
    end
    return diag_terms, new_ne
end

function _get_indices_from_terms(s::QAdd)
    inds = Index[]
    for (ops, _) in s.arguments
        for op in ops
            idx = op.index
            has_index(idx) && idx ∉ inds && push!(inds, idx)
        end
    end
    return inds
end

"""
    _apply_diagonal_split(s::QAdd, ext_idx::Index) -> QAdd

Apply diagonal splitting for each summation index in `s` that shares a space with `ext_idx`.
Called after multiplication introduces `ext_idx` into the terms.
"""
function _apply_diagonal_split(s::QAdd, ext_idx::Index)
    result = s
    for sum_idx in s.indices
        if sum_idx != ext_idx &&
                sum_idx.space_index == ext_idx.space_index &&
                !((sum_idx, ext_idx) in result.non_equal) &&
                !((ext_idx, sum_idx) in result.non_equal)
            diag_terms, new_ne = _diagonal_split(result, sum_idx, result.non_equal)
            result = QAdd(result.arguments, result.indices, new_ne)
            if !isempty(diag_terms)
                result = result + _qadd_from_oterms(diag_terms, Index[], Tuple{Index, Index}[])
            end
        end
    end
    return result
end

# ============================================================================
#  Σ — symbolic sums
# ============================================================================

"""
    Σ(expr, i::Index; non_equal::Vector{Index}=Index[])
    Σ(expr, i::Index, j::Index, rest::Index...)
    ∑(expr, i::Index, ...)

Create a symbolic sum ``\\sum_{i=1}^{N}`` of `expr` over index `i`.

Returns a [`QAdd`](@ref) with summation metadata. If `expr` does not depend on `i`,
the sum is evaluated eagerly as `range * expr`.

Diagonal splitting is performed automatically: when `expr` contains a free index `j`
on the same subspace as `i`, the sum is split into off-diagonal (`i ≠ j`) and diagonal
(`i = j`) contributions.

# Arguments
- `expr` — the expression to sum (a [`QAdd`](@ref), [`QSym`](@ref), or `Number`)
- `i::Index` — the summation index
- `non_equal::Vector{Index}` — indices that `i` must not equal (pairwise constraints)

Multiple indices can be passed to create nested sums: `Σ(expr, i, j)` is equivalent
to `Σ(Σ(expr, i), j)`.

Unicode input: `\\sum<tab>`. ASCII alias: `∑ = Σ`.

# Examples
```julia
h = FockSpace(:site) ⊗ FockSpace(:cavity)
i = Index(h, :i, N, 1)
@qnumbers a::Destroy(h, 1)
a_i = IndexedOperator(a, i)
H = Σ(a_i' * a_i, i)           # Σ(i=1:N) a†_i a_i
```

See also [`expand_sums`](@ref), [`Index`](@ref), [`IndexedOperator`](@ref).
"""
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    if !_any_depends_on_index(expr, i)
        return expr * i.range
    end
    ne_pairs = Tuple{Index, Index}[(i, j) for j in non_equal]
    all_ne = vcat(expr.non_equal, ne_pairs)
    diag_terms, all_ne = _diagonal_split(expr, i, all_ne)
    all_indices = vcat(expr.indices, [i])
    off_diag = QAdd(expr.arguments, all_indices, all_ne)
    isempty(diag_terms) && return off_diag
    return off_diag + _qadd_from_oterms(diag_terms, Index[], Tuple{Index, Index}[])
end

function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_CNUM_ONE, QSym[expr]), i, non_equal)
end
function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_to_cnum(expr), QSym[]), i, non_equal)
end

function Σ(expr, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ

"""
    expand_sums(expr)

Identity function — returns `expr` unchanged.

Diagonal splitting of symbolic sums is now performed eagerly at construction time
by [`Σ`](@ref) and `QAdd` multiplication. This function is retained for backward
compatibility.
"""
expand_sums(s::QAdd) = s
expand_sums(op::QSym) = op
