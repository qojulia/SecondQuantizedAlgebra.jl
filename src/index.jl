"""
    Index

Symbolic summation index for many-body systems.

Fields:
- `name::Symbol` — display name
- `range::Num` — upper bound (Int or symbolic)
- `space_index::Int` — which space in ProductSpace
- `sym::Num` — Symbolics symbolic integer for variable substitution
"""
struct Index
    name::Symbol
    range::Num
    space_index::Int
    sym::Num
end

const NO_INDEX = Index(:_, Num(0), 0, Num(0))
has_index(idx::Index) = idx.space_index != 0

function Base.:(==)(a::Index, b::Index)
    return a.name == b.name && isequal(a.range, b.range) && a.space_index == b.space_index
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
