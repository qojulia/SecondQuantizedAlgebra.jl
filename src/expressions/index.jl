"""
    IndexedVariable(name::Symbol, i::Index) -> Num

Create a symbolic single-indexed parameter ``\\mathrm{name}(i)``.

Use this for site-dependent quantities in Hamiltonians, for example detunings
or local couplings.

# Examples
```jldoctest
julia> h = FockSpace(:site) ⊗ FockSpace(:cavity);

julia> i = Index(h, :i, 10, 1);

julia> ω = IndexedVariable(:ω, i)
ω(i)
```

See also [`DoubleIndexedVariable`](@ref).
"""
function IndexedVariable(name::Symbol, i::Index)
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        name;
        type = SymbolicUtils.FnType{Tuple{Int}, Real, Nothing},
        shape = UnitRange{Int}[]
    )
    return Num(f(SymbolicUtils.unwrap(i.sym)))
end

"""Metadata key: marks a `DoubleIndexedVariable` node where equal indices must give zero."""
struct NotIdentical end

"""
    DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical::Bool=true) -> Num

Create a symbolic two-index parameter ``\\mathrm{name}(i, j)``.

Use this for pairwise interactions such as hopping amplitudes or coupling
matrices.

# Keyword arguments
- `identical::Bool = true` — if `false`, the variable evaluates to zero when `i == j`,
  enforcing that the parameter is only defined for distinct sites.

# Examples
```jldoctest
julia> h = FockSpace(:site) ⊗ FockSpace(:cavity);

julia> i = Index(h, :i, 10, 1);

julia> j = Index(h, :j, 10, 1);

julia> J = DoubleIndexedVariable(:J, i, j; identical = false)
J(i, j)
```

See also [`IndexedVariable`](@ref).
"""
function DoubleIndexedVariable(
        name::Symbol, i::Index, j::Index;
        identical::Bool = true
    )
    if !identical && i == j
        return Num(0)
    end
    f = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
        name;
        type = SymbolicUtils.FnType{Tuple{Int, Int}, Real, Nothing},
        shape = UnitRange{Int}[]
    )
    node = f(SymbolicUtils.unwrap(i.sym), SymbolicUtils.unwrap(j.sym))
    if !identical
        node = SymbolicUtils.setmetadata(node, NotIdentical, true)
    end
    return Num(node)
end

"""
    change_index(expr, from::Index, to::Index)

Rename an index in an expression by replacing `from` with `to` everywhere.

Operator indices are swapped directly. Symbolic prefactors (e.g. [`IndexedVariable`](@ref),
[`DoubleIndexedVariable`](@ref)) are substituted via `Symbolics.substitute` using
the `sym` fields of the indices. `DoubleIndexedVariable` nodes with `identical=false`
automatically evaluate to zero if the substitution makes both arguments equal.

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 5, h); j = Index(h, :j, 5, h);

julia> expr = IndexedVariable(:ω, i) * IndexedOperator(a, i);

julia> change_index(expr, i, j)
ω(j) * a_j
```

See also [`get_indices`](@ref).
"""
change_index(x::Number, ::Index, ::Index) = x
function change_index(x::Num, from::Index, to::Index)
    raw = SymbolicUtils.unwrap(x)
    isym = SymbolicUtils.unwrap(from.sym)
    vars = Symbolics.get_variables(raw)
    any(v -> isequal(v, isym), vars) || return x
    result = Symbolics.substitute(raw, Dict(isym => SymbolicUtils.unwrap(to.sym)))
    return _check_not_identical(result)
end
function change_index(x::CNum, from::Index, to::Index)
    return Complex(change_index(real(x), from, to), change_index(imag(x), from, to))
end

function change_index(op::Destroy, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Destroy(op.name, op.space_index, idx)
end
function change_index(op::Create, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Create(op.name, op.space_index, idx)
end
function change_index(op::Transition, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Transition(op.name, op.i, op.j, op.space_index, idx, op.ground_state, op.n_levels)
end
function change_index(op::Pauli, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Pauli(op.name, op.axis, op.space_index, idx)
end
function change_index(op::Spin, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Spin(op.name, op.axis, op.space_index, idx)
end
function change_index(op::Position, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Position(op.name, op.space_index, idx)
end
function change_index(op::Momentum, from::Index, to::Index)
    idx = op.index == from ? to : op.index
    return Momentum(op.name, op.space_index, idx)
end

function change_index(s::QAdd, from::Index, to::Index)
    out = QTermDict()
    needs = !isempty(s.indices)
    for (term, c) in s.arguments
        new_c = change_index(c, from, to)
        new_ops = QSym[change_index(op, from, to) for op in term.ops]
        new_ne = _substitute_ne(term.ne, from, to)
        if needs
            _accumulate_with_diag!(out, new_ops, new_c, s.indices, new_ne)
        else
            _canonicalize!(out, new_ops, new_c, new_ne)
        end
    end
    new_indices = Index[idx == from ? to : idx for idx in s.indices]
    return QAdd(out, new_indices)
end

"""
    get_indices(expr) -> Vector{Index}

Return all symbolic indices that appear in `expr`.

Returns a `Vector{Index}` of unique indices found in operator fields and
summation metadata. Excludes the sentinel `NO_INDEX`.

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 5, h);

julia> get_indices(IndexedOperator(a, i))
1-element Vector{Index}:
 i
```
"""
get_indices(::Number) = Index[]
get_indices(x::Num) = get_indices(SymbolicUtils.unwrap(x))
function get_indices(op::QSym)
    return has_index(op.index) ? Index[op.index] : Index[]
end
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for term in keys(s.arguments)
        for op in term.ops
            has_index(op.index) && op.index ∉ inds && push!(inds, op.index)
        end
    end
    return inds
end

"""
    create_index_arrays(indices::Vector{Index}, ranges::Vector{<:AbstractRange}) -> Vector

Cartesian product of `ranges`, returned as a flat vector.
Single-index case returns `collect(only(ranges))`.
"""
function create_index_arrays(indices::Vector{Index}, ranges::Vector{<:AbstractRange})
    length(indices) == 1 && return collect(only(ranges))
    return vec(collect(Iterators.product(ranges...)))
end

"""Check if a substituted BasicSymbolic node with NotIdentical metadata has equal args → 0."""
function _check_not_identical(result::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(result) &&
            SymbolicUtils.hasmetadata(result, NotIdentical) &&
            length(SymbolicUtils.arguments(result)) == 2
        a1, a2 = SymbolicUtils.arguments(result)
        isequal(a1, a2) && return Num(0)
    end
    return Num(result)
end
_check_not_identical(result::Number) = Num(result)

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
    for (term, c) in s.arguments
        _depends_on_index_term(c, term.ops, idx) && return true
    end
    return false
end

"""
    _diagonal_split!(off_diag, diag, sum_idx) -> nothing

For every term in `off_diag` that depends on `sum_idx`, locate each free index
`j` on the same Hilbert subspace (and not already constrained `sum_idx ≠ j`),
emit the diagonal substitution `sum_idx → j` into `diag`, and re-key the
off-diagonal entry under the augmented constraint `ne ∪ {(sum_idx, j)}`.
Mutates both dicts in place; the caller composes them via `+`.
"""
function _diagonal_split!(off_diag::QTermDict, diag::QTermDict, sum_idx::Index)
    seen = Index[]
    for (term, c) in collect(off_diag)
        _depends_on_index_term(c, term.ops, sum_idx) || continue
        empty!(seen)
        current_term = term
        current_c = c
        for op in term.ops
            idx = op.index
            has_index(idx) || continue
            idx ∈ seen && continue
            push!(seen, idx)
            idx == sum_idx && continue
            idx.space_index == sum_idx.space_index || continue
            _ne_contains(current_term.ne, sum_idx, idx) && continue

            new_c = change_index(current_c, sum_idx, idx)
            new_ops = QSym[change_index(o, sum_idx, idx) for o in current_term.ops]
            _canonicalize!(diag, new_ops, new_c, _drop_ne_with(current_term.ne, sum_idx))

            delete!(off_diag, current_term)
            current_term = _term_key(
                current_term.ops,
                _merge_ne_pair(current_term.ne, sum_idx, idx),
            )
            _addto_key!(off_diag, current_term, current_c)
        end
    end
    return nothing
end

"""
    Σ(expr, i::Index, non_equal::Vector{Index} = Index[])
    Σ(expr, i::Index, j::Index, rest::Index...)
    ∑(expr, i::Index, ...)

Build symbolic sums over indices, for example ``\\sum_i expr``.

Returns a [`QAdd`](@ref) carrying `i` in its `indices` field. If `expr` does not
depend on `i`, the sum is evaluated eagerly as `i.range * expr`.

The optional `non_equal` records pairwise inequality constraints `i ≠ j` per
term. Diagonal splitting is performed automatically: if `expr` carries another
free index `j` on the same Hilbert subspace as `i`, the contribution at `i = j`
is emitted as a separate diagonal term and the off-diagonal term gains the
constraint `(i, j)`.

Multiple positional indices create nested sums: `Σ(expr, i, j)` is equivalent
to `Σ(Σ(expr, i), j)`. The Unicode alias `∑` is also exported.

# Examples

```jldoctest
julia> h = FockSpace(:site);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 3, h);

julia> Σ(IndexedOperator(a', i) * IndexedOperator(a, i), i)
Σ(i=1:3) a_i' * a_i
```

See also [`Index`](@ref), [`IndexedOperator`](@ref), [`constraint_pairs`](@ref).
"""
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    if !_any_depends_on_index(expr, i)
        return expr * i.range
    end

    off_diag = QTermDict()
    diag = QTermDict()

    for (term, c) in expr.arguments
        # User-supplied non_equal constraints apply when the term depends on i.
        ne_aug = term.ne
        if _depends_on_index_term(c, term.ops, i)
            for j in non_equal
                ne_aug = _merge_ne_pair(ne_aug, i, j)
            end
        end

        # Identify (i, ext_idx) pairs that need diagonal contributions.
        distinct = _distinct_op_indices(term.ops)
        diag_pairs = Tuple{Index, Index}[]
        if _depends_on_index_term(c, term.ops, i)
            for ext_idx in distinct
                ext_idx == i && continue
                ext_idx.space_index == i.space_index || continue
                _ne_contains(ne_aug, i, ext_idx) && continue
                push!(diag_pairs, (i, ext_idx))
            end
        end

        # Off-diagonal: canonicalize under ne_aug augmented with every diag pair,
        # so partial_sort can use those inequalities.
        full_ne = ne_aug
        for (a, b) in diag_pairs
            full_ne = _merge_ne_pair(full_ne, a, b)
        end
        _canonicalize!(off_diag, copy(term.ops), c, full_ne)

        for (_, ext_idx) in diag_pairs
            sub_ops = QSym[change_index(o, i, ext_idx) for o in term.ops]
            sub_c = change_index(c, i, ext_idx)
            sub_ne = _drop_ne_with(term.ne, i)
            _canonicalize!(diag, sub_ops, sub_c, sub_ne)
        end
    end

    result = QAdd(off_diag, vcat(expr.indices, [i]))
    isempty(diag) && return result
    return result + QAdd(diag, Index[])
end

function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_CNUM_ONE, QSym[expr]), i, non_equal)
end

function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(_single_qadd(_to_cnum(expr), QSym[]), i, non_equal)
end

function Σ(expr::Union{QField, Number}, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ
