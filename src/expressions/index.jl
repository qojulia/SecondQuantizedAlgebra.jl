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
    return Num(f(SymbolicUtils.unwrap(index_sym(i))))
end
IndexedVariable(name::AbstractString, args...) = name_must_be_symbol(name)

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
    node = f(SymbolicUtils.unwrap(index_sym(i)), SymbolicUtils.unwrap(index_sym(j)))
    if !identical
        node = SymbolicUtils.setmetadata(node, NotIdentical, true)
    end
    return Num(node)
end
DoubleIndexedVariable(name::AbstractString, args...; kwargs...) = name_must_be_symbol(name)

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
    isym = SymbolicUtils.unwrap(index_sym(from))
    vars = Symbolics.get_variables(raw)
    any(v -> isequal(v, isym), vars) || return x
    result = Symbolics.substitute(raw, Dict(isym => SymbolicUtils.unwrap(index_sym(to))))
    return check_not_identical(result)
end
function change_index(x::CNum, from::Index, to::Index)
    is_native(x) && return x
    tail = x.tail
    if tail isa RawSymbolicCoeff
        raw = change_index(Num(tail.expr), from, to)
        isequal(raw, Num(tail.expr)) && return x
        result = SymbolicUtils.unwrap(raw)
        is_numeric_const(result) && return const_coeff(result)
        return from_raw_arithmetic(result, tail.real_slot)
    end
    raw_re, raw_im = realimag(x)
    re = change_index(raw_re, from, to)
    im = change_index(raw_im, from, to)
    return cnum(re, im)
end

with_index(op::Op, i::Index) = IndexedOperator(op, i)

function change_index(op::QSym, from::Index, to::Index)
    return with_index(op, op.index == from ? to : op.index)
end

function change_index(s::QAdd, from::Index, to::Index)
    out = QTermDict()
    needs = !isempty(s.indices)
    for (term, c) in s.arguments
        ne_becomes_contradictory(term.ne, from, to) && continue
        new_c = change_index(c, from, to)
        new_ops = Op[change_index(op, from, to) for op in term.ops]
        new_ne = substitute_ne(term.ne, from, to)
        if needs
            accumulate_with_diag!(out, new_ops, new_c, s.indices, new_ne)
        else
            canonicalize!(out, new_ops, new_c, new_ne)
        end
    end
    new_indices = Index[idx == from ? to : idx for idx in s.indices]
    return QAdd(out, new_indices)
end

"""
    change_index(expr, pairs::AbstractDict{Index, Index})

Simultaneous (batched) variant of [`change_index`](@ref): substitute every
key in `pairs` by its mapped value in one pass.

Unlike a sequence of single-pair `change_index` calls, this never produces
an intermediate state where a renamed-into index temporarily collides with
another encountered index. Use it when the rename describes a permutation
or any map whose image overlaps its domain, e.g. swapping two indices:

```julia
julia> change_index(expr, Dict(i => j, j => i))   # swap, not fuse
```

A pair `(idx => idx)` is a no-op; `change_index(expr, Dict())` returns
`expr` unchanged.
"""
change_index(x::Number, ::AbstractDict{Index, Index}) = x
function change_index(x::Num, pairs::AbstractDict{Index, Index})
    isempty(pairs) && return x
    raw = SymbolicUtils.unwrap(x)
    sub = Dict{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic}()
    for (from, to) in pairs
        sub[SymbolicUtils.unwrap(index_sym(from))] = SymbolicUtils.unwrap(index_sym(to))
    end
    vars = Symbolics.get_variables(raw)
    any(v -> haskey(sub, v), vars) || return x
    result = Symbolics.substitute(raw, sub)
    return check_not_identical(result)
end
function change_index(x::CNum, pairs::AbstractDict{Index, Index})
    is_native(x) && return x
    tail = x.tail
    if tail isa RawSymbolicCoeff
        raw = change_index(Num(tail.expr), pairs)
        isequal(raw, Num(tail.expr)) && return x
        result = SymbolicUtils.unwrap(raw)
        is_numeric_const(result) && return const_coeff(result)
        return from_raw_arithmetic(result, tail.real_slot)
    end
    raw_re, raw_im = realimag(x)
    re = change_index(raw_re, pairs)
    im = change_index(raw_im, pairs)
    return cnum(re, im)
end

function change_index(op::QSym, pairs::AbstractDict{Index, Index})
    new_idx = get(pairs, op.index, op.index)
    return new_idx === op.index ? op : with_index(op, new_idx)
end

function change_index(s::QAdd, pairs::AbstractDict{Index, Index})
    isempty(pairs) && return s
    out = QTermDict()
    needs = !isempty(s.indices)
    for (term, c) in s.arguments
        ne_becomes_contradictory(term.ne, pairs) && continue
        new_c = change_index(c, pairs)
        new_ops = Op[change_index(op, pairs) for op in term.ops]
        new_ne = substitute_ne(term.ne, pairs)
        if needs
            accumulate_with_diag!(out, new_ops, new_c, s.indices, new_ne)
        else
            canonicalize!(out, new_ops, new_c, new_ne)
        end
    end
    new_indices = Index[get(pairs, idx, idx) for idx in s.indices]
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

"""
Zero every `NotIdentical`-tagged node anywhere in a substituted expression whose
two indices became equal. Nested occurrences matter: a `DoubleIndexedVariable`
with `identical=false` is almost always a factor in a larger product or sum, and
the surrounding `*`/`+` only collapses once that factor is replaced by zero.
"""
function check_not_identical(result::SymbolicUtils.BasicSymbolic)
    zeros = collect_identical_zeros!(Dict{SymbolicUtils.BasicSymbolic, Int}(), result)
    isempty(zeros) && return Num(result)
    return Num(Symbolics.substitute(result, zeros))
end
check_not_identical(result::Number) = Num(result)

"""Accumulate into `acc` every `NotIdentical` node in `x` with two equal args, mapped to 0."""
function collect_identical_zeros!(acc, x)
    (x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x)) || return acc
    if SymbolicUtils.hasmetadata(x, NotIdentical) &&
            length(SymbolicUtils.arguments(x)) == 2
        a1, a2 = SymbolicUtils.arguments(x)
        isequal(a1, a2) && (acc[x] = 0; return acc)
    end
    for a in SymbolicUtils.arguments(x)
        collect_identical_zeros!(acc, a)
    end
    return acc
end

function contains_leaf(x::RawExpression, v::RawExpression)
    isequal(x, v) && return true
    SymbolicUtils.iscall(x) || return false
    for a in SymbolicUtils.arguments(x)
        a isa RawExpression && contains_leaf(a, v) && return true
    end
    return false
end

function depends_on_index_term(c::CNum, ops::Vector{Op}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    is_native(c) && return false
    isym = SymbolicUtils.unwrap(index_sym(idx))
    tail = c.tail
    if tail isa Poly
        for m in tail.terms, s in m.syms
            contains_leaf(s, isym) && return true
        end
        return false
    end
    tail isa RawSymbolicCoeff || return false
    return contains_leaf(tail.expr, isym)
end

# A summation index is also observable through a counting constraint: Σ_{i≠j} T has N-1
# terms even when T itself does not mention i.
function depends_on_scope_index(
        c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair}, idx::Index,
    )
    depends_on_index_term(c, ops, idx) && return true
    for (α, β) in ne
        (α == idx || β == idx) && return true
    end
    return false
end

"""
    eliminate_binder(ne, idx) -> (eliminated, factor, reduced_ne)

Count the values of the summation index `idx` in a term that no longer mentions it through
its operators or coefficient. With `X` the indices on the same space that `ne` keeps distinct
from `idx`, the count is `index_range(idx) - |X|`. It is only known when the members of `X`
are themselves provably distinct; otherwise `eliminated` is `false` and the term must keep
the binder together with its constraints.
"""
function eliminate_binder(ne::Vector{NonEqualPair}, idx::Index)
    partners = Index[]
    mentions = false
    for (α, β) in ne
        other = α == idx ? β : β == idx ? α : idx
        other == idx && continue
        mentions = true
        other.space_index == idx.space_index && push!(partners, other)
    end
    for p in 1:length(partners), q in (p + 1):length(partners)
        ne_contains(ne, partners[p], partners[q]) ||
            return (false, CNUM_ONE, ne)
    end
    factor = to_cnum(index_range(idx) - length(partners))
    mentions || return (true, factor, ne)
    reduced = NonEqualPair[p for p in ne if p[1] != idx && p[2] != idx]
    return (true, factor, isempty(reduced) ? EMPTY_NE : reduced)
end

"""
    eliminate_binders(mentions, coef, ne, scope) -> (coef, ne)

Eliminate every binder of `scope` for which `mentions(coef, idx)` is `false`, multiplying
`coef` by its constrained count. Elimination runs to a fixpoint: removing one binder drops
its constraints, which can make the count of a binder kept earlier known.
"""
function eliminate_binders(
        mentions::F, coef::CNum, ne::Vector{NonEqualPair}, scope::Vector{Index},
    ) where {F}
    pending = EMPTY_INDICES
    for idx in scope
        mentions(coef, idx) && continue
        eliminated, factor, ne = eliminate_binder(ne, idx)
        if eliminated
            coef = mul_cnum(coef, factor)
        else
            pending === EMPTY_INDICES && (pending = Index[])
            push!(pending, idx)
        end
    end
    progress = true
    while progress && !isempty(pending)
        progress = false
        k = 1
        while k <= length(pending)
            eliminated, factor, ne = eliminate_binder(ne, pending[k])
            if eliminated
                coef = mul_cnum(coef, factor)
                deleteat!(pending, k)
                progress = true
            else
                k += 1
            end
        end
    end
    return (coef, ne)
end

function any_depends_on_index(s::QAdd, idx::Index)
    for (term, c) in s.arguments
        depends_on_scope_index(c, term.ops, term.ne, idx) && return true
    end
    return false
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
    i.range_id == 0 && throw(ArgumentError("Σ: index $(index_name(i)) has no summation range"))
    if isempty(non_equal) && !any_depends_on_index(expr, i)
        return expr * index_range(i)
    end

    off_diag = QTermDict()
    diag = QTermDict()

    for (term, c) in expr.arguments
        sigma_term!(off_diag, diag, term, c, i, non_equal)
    end

    binds = any(((t, c),) -> depends_on_scope_index(c, t.ops, t.ne, i), off_diag)
    result = QAdd(off_diag, binds ? vcat(expr.indices, [i]) : expr.indices)
    isempty(diag) && return result
    return result + QAdd(diag, Index[])
end

# Split one term of `Σ_i` into its off-diagonal part (`i` distinct from every other index
# on its space) and one diagonal part per index `i` may coincide with.
function sigma_term!(
        off_diag::QTermDict, diag::QTermDict, term::QTerm, c::CNum,
        i::Index, non_equal::Vector{Index},
    )
    ne_aug = term.ne
    for j in non_equal
        ne_aug = merge_ne_pair(ne_aug, i, j)
    end

    if !depends_on_index_term(c, term.ops, i)
        emit_scaled_by_scope!(off_diag, copy(term.ops), c, ne_aug, Index[i])
        return nothing
    end

    # Identify (i, ext_idx) pairs that need diagonal contributions.
    diag_pairs = Tuple{Index, Index}[]
    for ext_idx in distinct_op_indices(term.ops)
        may_coincide(i, ext_idx, ne_aug) && push!(diag_pairs, (i, ext_idx))
    end

    # Off-diagonal: canonicalize under ne_aug augmented with every diag pair,
    # so partial_sort can use those inequalities.
    full_ne = ne_aug
    for (a, b) in diag_pairs
        full_ne = merge_ne_pair(full_ne, a, b)
    end
    # Σ_i: bake `i.range` into emitted entries whose final ops lose their
    # i-dependence (e.g. the `+1` residual from `a_i*a_i' = a_i'*a_i + 1`),
    # so per-term scope matches the explicit `Σ_i` written by the user.
    emit_scaled_by_scope!(off_diag, copy(term.ops), c, full_ne, Index[i])

    for (_, ext_idx) in diag_pairs
        sub_ops = Op[change_index(o, i, ext_idx) for o in term.ops]
        sub_c = change_index(c, i, ext_idx)
        # Propagate NE onto the diagonal (i≠β becomes ext_idx≠β); dropping it
        # lets a later sum re-admit ext_idx=β and double-count the diagonal.
        sub_ne = diagonal_ne(ne_aug, diag_pairs, i, ext_idx)
        canonicalize!(diag, sub_ops, sub_c, sub_ne)
    end
    return nothing
end

function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(single_qadd(CNUM_ONE, Op[expr]), i, non_equal)
end

function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(single_qadd(to_cnum(expr), Op[]), i, non_equal)
end

function Σ(expr::Union{QField, Number}, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ
