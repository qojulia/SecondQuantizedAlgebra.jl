"""
    IndexedVariable(name::Symbol, i::Index) -> Num

Create a symbolic single-indexed variable ``\\mathrm{name}(i)`` as a `Num` from Symbolics.jl.

Useful for site-dependent parameters in indexed Hamiltonians (e.g. detunings, couplings).

# Examples
```julia
h = FockSpace(:site) ⊗ FockSpace(:cavity)
i = Index(h, :i, 10, 1)
ω = IndexedVariable(:ω, i)    # ω(i) — symbolic parameter depending on site i
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

Create a symbolic double-indexed variable ``\\mathrm{name}(i, j)`` as a `Num` from Symbolics.jl.

Useful for pairwise interaction parameters (e.g. coupling matrices, hopping amplitudes).

# Keyword arguments
- `identical::Bool = true` — if `false`, the variable evaluates to zero when `i == j`,
  enforcing that the parameter is only defined for distinct sites.

# Examples
```julia
h = FockSpace(:site) ⊗ FockSpace(:cavity)
i = Index(h, :i, 10, 1)
j = Index(h, :j, 10, 1)
J = DoubleIndexedVariable(:J, i, j; identical=false)  # J(i,j), zero when i==j
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

Replace every occurrence of index `from` with `to` throughout an expression tree.

Operator indices are swapped directly. Symbolic prefactors (e.g. [`IndexedVariable`](@ref),
[`DoubleIndexedVariable`](@ref)) are substituted via `Symbolics.substitute` using
the `sym` fields of the indices. `DoubleIndexedVariable` nodes with `identical=false`
automatically evaluate to zero if the substitution makes both arguments equal.

See also [`get_indices`](@ref).
"""
change_index(x::Number, ::Index, ::Index) = x
function change_index(x::Num, from::Index, to::Index)
    raw = SymbolicUtils.unwrap(x)
    result = Symbolics.substitute(
        raw,
        Dict(SymbolicUtils.unwrap(from.sym) => SymbolicUtils.unwrap(to.sym))
    )
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
    new_d = QTermDict()
    for (ops, c) in s.arguments
        new_c = change_index(c, from, to)
        new_ops = QSym[change_index(op, from, to) for op in ops]
        _addto!(new_d, new_ops, new_c)
    end
    new_indices = [idx == from ? to : idx for idx in s.indices]
    new_ne = Tuple{Index, Index}[
        (a == from ? to : a, b == from ? to : b)
            for (a, b) in s.non_equal
    ]
    return QAdd(new_d, new_indices, new_ne)
end

"""
    get_indices(expr) -> Vector{Index}

Collect all symbolic summation indices appearing in `expr`.

Returns a `Vector{Index}` of unique indices found in operator fields and
summation metadata. Excludes the sentinel `NO_INDEX`.
"""
get_indices(::Number) = Index[]
get_indices(x::Num) = get_indices(SymbolicUtils.unwrap(x))
function get_indices(op::QSym)
    return has_index(op.index) ? Index[op.index] : Index[]
end
function get_indices(s::QAdd)
    inds = copy(s.indices)
    for (ops, _) in s.arguments
        for op in ops
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

# --- Shared helper for NotIdentical check ---

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
