"""
    IndexedVariable(name::Symbol, i::Index)

Site-dependent parameter `name(i)`, such as a local detuning or coupling.

The parameter is real. Under a sum it follows the index, so it is evaluated site by site
once the sum is unrolled.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> g = IndexedVariable(:g, i)
g(i)

julia> Σ(g * IndexedOperator(Transition(h, :σ, 2, 1), i), i)
Σ(i=1:10) g(i) * σ_i₂₁
```

See also [`DoubleIndexedVariable`](@ref), [`Index`](@ref), [`index_sym`](@ref).
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
    DoubleIndexedVariable(name::Symbol, i::Index, j::Index; identical = true)

Pair parameter `name(i, j)`, such as a hopping amplitude or a dipole-dipole coupling.

With `identical = false` the parameter vanishes whenever both indices refer to the same
site, so the diagonal terms of a double sum drop out.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> J = DoubleIndexedVariable(:J, i, j; identical = false);

julia> Σ(J * σ(2, 1, i) * σ(1, 2, j), i, j)
Σ(i=1:10)Σ(j=1:10)(i≠j) J(i, j) * σ_i₂₁ * σ_j₁₂
```

See also [`IndexedVariable`](@ref), [`Σ`](@ref).
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
    change_index(expr, pairs::AbstractDict{Index, Index})

Replace the index `from` by `to` in operators and coefficients of `expr`.

The result is brought back into canonical form, so operators that now share a site
compose, and a [`DoubleIndexedVariable`](@ref) declared with `identical = false` vanishes
when its two indices become equal. The second form renames several indices at once, which
makes swaps such as `Dict(i => j, j => i)` possible.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> change_index(σ(2, 1, i) * σ(1, 2, j), i, j)
σ_j₂₂

julia> change_index(σ(2, 1, i) * σ(1, 2, j), Dict(i => j, j => i))
σ_j₂₁ * σ_i₁₂
```

See also [`get_indices`](@ref), [`substitute`](@ref).
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
        return raw_coeff(result, tail.known_real)
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
        return raw_coeff(result, tail.known_real)
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

Indices occurring in `expr`, both summed and free.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h); j = Index(h, :j, 10, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> get_indices(Σ(σ(2, 1, i) * σ(1, 2, j), i))
2-element Vector{Index}:
 i
 j
```

See also [`change_index`](@ref), [`has_index`](@ref).
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


function depends_on_index_term(c::CNum, ops::Vector{Op}, idx::Index)
    for op in ops
        op.index == idx && return true
    end
    return coefficient_depends_on(c, SymbolicUtils.unwrap(index_sym(idx)))
end

function depends_on_scope_index(
        c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair}, idx::Index,
    )
    depends_on_index_term(c, ops, idx) && return true
    for (α, β) in ne
        (α == idx || β == idx) && return true
    end
    return false
end

bound_indices(c::CNum, ops::Vector{Op}, ne::Vector{NonEqualPair}, indices::Vector{Index}) =
    Index[idx for idx in indices if depends_on_scope_index(c, ops, ne, idx)]
bound_indices(t::QTerm, c::CNum, indices::Vector{Index}) = bound_indices(c, t.ops, t.ne, indices)

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
        (α == idx || β == idx) || continue
        mentions = true
        other = α == idx ? β : α
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
            coef = coef * factor
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
                coef = coef * factor
                deleteat!(pending, k)
                progress = true
            else
                k += 1
            end
        end
    end
    return (coef, ne)
end

function any_depends_on_index(d::QTermDict, idx::Index)
    for (term, c) in d
        depends_on_scope_index(c, term.ops, term.ne, idx) && return true
    end
    return false
end
any_depends_on_index(s::QAdd, idx::Index) = any_depends_on_index(s.arguments, idx)

"""
    Σ(expr, i::Index, [non_equal::Vector{Index}])
    Σ(expr, i::Index, j::Index, more::Index...)
    ∑(expr, i, ...)

Symbolic sum of `expr` over all sites of the index `i`.

Several indices give a nested sum. The sum stays symbolic and follows the algebra like any
other expression. Where its index meets another index on the same subsystem, the terms in
which the two sites coincide are split off, so the remaining sum carries the constraint
`i ≠ j`. The optional `non_equal` adds such constraints by hand, excluding those sites from
the sum. A summand that does not depend on `i` is multiplied by the number of sites.

`∑` (typed `\\sum<tab>`) is the same function.

# Examples
```jldoctest
julia> @variables N;

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, N, h); j = Index(h, :j, N, h);

julia> σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m);

julia> Σ(σ(2, 1, i) * σ(1, 2, j), i, j)
Σ(j=1:N) σ_j₂₂ + Σ(i=1:N)Σ(j=1:N)(i≠j) σ_i₂₁ * σ_j₁₂

julia> Σ(σ(2, 2, j), i, [j])
(-1 + N) * σ_j₂₂
```

See also [`Index`](@ref), [`IndexedOperator`](@ref), [`assume_distinct_index`](@ref).
"""
function Σ(expr::QAdd, i::Index, non_equal::Vector{Index} = Index[])
    i.range_id == 0 && throw(ArgumentError("Σ: index $(index_name(i)) has no summation range"))
    if isempty(non_equal) && !any_depends_on_index(expr, i)
        return expr * index_range(i)
    end

    off_diag = QTermDict()
    diag = QTermDict()

    for (term, c) in expr.arguments
        scope = vcat(bound_indices(term, c, expr.indices), [i])
        sigma_term!(off_diag, diag, term, c, i, non_equal, scope)
    end

    result = QAdd(off_diag, drop_unused_indices(off_diag, vcat(expr.indices, [i])))
    isempty(diag) && return result
    return result + QAdd(diag, drop_unused_indices(diag, copy(expr.indices)))
end

function sigma_term!(
        off_diag::QTermDict, diag::QTermDict, term::QTerm, c::CNum,
        i::Index, non_equal::Vector{Index}, scope::Vector{Index},
    )
    ne_aug = term.ne
    for j in non_equal
        ne_aug = merge_ne_pair(ne_aug, i, j)
    end

    if !depends_on_index_term(c, term.ops, i)
        emit_scaled_by_scope!(off_diag, copy(term.ops), c, ne_aug, scope)
        return nothing
    end

    diag_pairs = Tuple{Index, Index}[]
    for ext_idx in distinct_op_indices(term.ops)
        may_coincide(i, ext_idx, ne_aug) && push!(diag_pairs, (i, ext_idx))
    end

    full_ne = ne_aug
    for (a, b) in diag_pairs
        full_ne = merge_ne_pair(full_ne, a, b)
    end
    emit_scaled_by_scope!(off_diag, copy(term.ops), c, full_ne, scope)

    for (_, ext_idx) in diag_pairs
        sub_ops = Op[change_index(o, i, ext_idx) for o in term.ops]
        sub_c = change_index(c, i, ext_idx)
        sub_ne = diagonal_ne(ne_aug, diag_pairs, i, ext_idx)
        canonicalize!(diag, sub_ops, sub_c, sub_ne)
    end
    return nothing
end

function Σ(expr::QSym, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(to_qadd(expr), i, non_equal)
end

function Σ(expr::Number, i::Index, non_equal::Vector{Index} = Index[])
    return Σ(single_qadd(to_cnum(expr), Op[]), i, non_equal)
end

function Σ(expr::Union{QField, Number}, i::Index, j::Index, rest::Index...)
    inner = Σ(expr, i)
    return Σ(inner, j, rest...)
end

const ∑ = Σ
