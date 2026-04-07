"""
    Index(h::HilbertSpace, name::Symbol, range, space)
    Index(h::HilbertSpace, name::Symbol, range, space_index::Int)

Symbolic summation index for many-body systems with identical subsystems.

An `Index` labels a summation variable that runs from `1` to `range` over
operators on a particular subspace. It is used with [`Σ`](@ref) to build
symbolic sums and with [`IndexedOperator`](@ref) to attach indices to operators.

# Fields
- `name::Symbol` — display name (e.g. `:i`, `:j`)
- `range::Num` — upper bound of the summation (integer or symbolic)
- `space_index::Int` — which subspace in a [`ProductSpace`](@ref)
- `sym::Num` — Symbolics symbolic variable for algebraic substitution

# Examples
```julia
h = FockSpace(:site) ⊗ FockSpace(:cavity)
i = Index(h, :i, 10, FockSpace(:site))   # index i = 1:10 on the site space
j = Index(h, :j, 10, 1)                  # same, using integer space index
```

See also [`has_index`](@ref), [`IndexedOperator`](@ref), [`Σ`](@ref),
[`change_index`](@ref), [`insert_index`](@ref).
"""
struct Index
    name::Symbol
    range::Num
    space_index::Int
    sym::Num
end

const NO_INDEX = Index(:_, Num(0), 0, Num(0))

"""
    has_index(idx::Index) -> Bool

Return `true` if `idx` is a real summation index (not the sentinel `NO_INDEX`).
"""
has_index(idx::Index) = idx.space_index != 0

function Base.:(==)(a::Index, b::Index)
    a === b && return true
    return a.name == b.name && a.space_index == b.space_index && isequal(a.range, b.range)
end
function Base.hash(a::Index, h::UInt)
    return hash(:Index, hash(a.name, hash(a.space_index, h)))
end

# Construction from HilbertSpace
function Index(h::HilbertSpace, name::Symbol, range, space::HilbertSpace)
    si = _find_space_index(h, space)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, Num(range), si, Num(sym_var))
end
function Index(h::HilbertSpace, name::Symbol, range, si::Int)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, Num(range), si, Num(sym_var))
end

# Find space index in ProductSpace
_find_space_index(::HilbertSpace, ::HilbertSpace) = 1
function _find_space_index(h::ProductSpace, space::HilbertSpace)
    for (i, s) in enumerate(h.spaces)
        actual = _unwrap_space(s)
        actual == space && return i
    end
    throw(ArgumentError("Space $space not found in $h"))
end
