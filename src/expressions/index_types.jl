"""
    Index(h::HilbertSpace, name::Symbol, range, space::HilbertSpace)
    Index(h::HilbertSpace, name::Symbol, range, subsystem::Int)
    i(k::Integer)

Symbolic label for one of many identical sites.

An index runs from `1` to `range`, an integer or a symbolic parameter such as `N`, over the
factor of `h` given either by the space itself or by its position. Attached to operators with
[`IndexedOperator`](@ref) and to parameters with [`IndexedVariable`](@ref), it is summed over
with [`Σ`](@ref). Two indices are equal only if they have the same name, range and subsystem.

Calling an index with an integer, `i(3)`, gives the index of a concrete site. It is named
`i_3`, keeps the range and subsystem of `i`, and addresses that site after a sum has been
written out site by site, for instance in numerical conversion or in cumulant equations.

# Examples
```jldoctest
julia> @variables N;

julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> i = Index(h, :i, N, NLevelSpace(:atom, 2));

julia> i == Index(h, :i, N, 2)
true

julia> σ = Transition(h, :σ, 2, 1);

julia> IndexedOperator(σ, i), IndexedOperator(σ, i(3))
(σ_i₂₁, σ_i_3₂₁)
```

See also [`IndexedOperator`](@ref), [`Σ`](@ref), [`index_name`](@ref), [`index_range`](@ref).
"""
struct Index
    name_id::Int32
    range_id::Int32
    space_index::Int32
    slot::Int32
end

const NO_INDEX = Index(Int32(0), Int32(0), Int32(0), Int32(0))

"""Metadata key for the site number of a concrete index such as `i(3)` on its symbol."""
struct IndexSlot end

"""
    has_index(i::Index) -> Bool

Whether `i` is an index rather than the placeholder carried by unindexed operators.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> σ = Transition(h, :σ, 2, 1);

julia> has_index(operator_index(σ)), has_index(operator_index(IndexedOperator(σ, i)))
(false, true)
```

See also [`operator_index`](@ref), [`get_indices`](@ref).
"""
has_index(idx::Index) = idx.space_index != 0

const NonEqualPair = Tuple{Index, Index}

const EMPTY_NE = NonEqualPair[]
const EMPTY_INDICES = Index[]

"""
    IndexMap

A renaming of indices: one `from => to` pair, or a dictionary for several at once. Indices
it does not mention map to themselves.
"""
const IndexMap = Union{Pair{Index, Index}, AbstractDict{Index, Index}}

reindex(m::Pair{Index, Index}, idx::Index) = idx == m.first ? m.second : idx
reindex(m::AbstractDict{Index, Index}, idx::Index) = get(m, idx, idx)

renames_nothing(::Pair{Index, Index}) = false
renames_nothing(m::AbstractDict{Index, Index}) = isempty(m)

function Base.:(==)(a::Index, b::Index)
    a === b && return true
    return a.name_id == b.name_id && a.space_index == b.space_index && a.range_id == b.range_id
end
function Base.hash(a::Index, h::UInt)
    return hash(:Index, hash(a.name_id, hash(a.space_index, h)))
end
index_key(idx::Index) = (idx.space_index, name_rank(idx.name_id), range_rank(idx.range_id), idx.slot)
Base.isless(a::Index, b::Index) = index_key(a) < index_key(b)

Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, space::HilbertSpace) =
    Index(h, name, range, find_space_index(h, space))
function Index(h::HilbertSpace, name::Symbol, range::Union{Int, Num}, si::Int)
    return Index(name_id(name), intern_range(range), Int32(si), Int32(0))
end
Index(h::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)

function (i::Index)(k::Integer)
    name = Symbol(name_from_id(i.name_id), "_", k)
    return Index(name_id(name), i.range_id, i.space_index, Int32(k))
end

"""
    index_slot(i::Index) -> Int
    index_slot(x) -> Int

Site number of a concrete index such as `i(3)`, or `0` for a symbolic index.

The second form reads the site number from the symbolic variable of an index, as found in
the coefficients of an expression.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> index_slot(i), index_slot(i(3))
(0, 3)
```

See also [`Index`](@ref), [`index_sym`](@ref).
"""
index_slot(idx::Index) = Int(idx.slot)
index_slot(x) = SymbolicUtils.getmetadata(SymbolicUtils.unwrap(x), IndexSlot, 0)

"""
    rename(i::Index, name::Symbol)
    rename(op::Op, name::Symbol)

Copy of an index or operator under a new name, with everything else kept.

For an index this gives a second index over the same sites, for instance the partner index
of a pair interaction. The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: rename

julia> @variables N;

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, N, h);

julia> j = rename(i, :j)
j

julia> j == Index(h, :j, N, h)
true
```

See also [`Index`](@ref), [`change_index`](@ref).
"""
rename(idx::Index, name::Symbol) = Index(name_id(name), idx.range_id, idx.space_index, idx.slot)
rename(idx::Index, name::AbstractString) = name_must_be_symbol(name)

"""
    index_range(i::Index)

Number of sites the index `i` runs over, as given at construction.

# Examples
```jldoctest
julia> @variables N;

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, N, h);

julia> index_range(i)
N
```

See also [`index_name`](@ref), [`Index`](@ref).
"""
index_range(idx::Index) = range_from_id(idx.range_id)

"""
    index_name(i::Index) -> Symbol

Name of the index `i`, as given at construction.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> index_name(i), index_name(i(3))
(:i, :i_3)
```

See also [`index_range`](@ref), [`index_sym`](@ref), [`Index`](@ref).
"""
index_name(idx::Index) = name_from_id(idx.name_id)

"""
    index_sym(i::Index) -> Num

Symbolic variable that stands for the index `i` in coefficients.

Indexed parameters such as `g(i)` are functions of this variable, so it is what to
substitute when a coefficient is evaluated at a given site.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> g = IndexedVariable(:g, i);

julia> substitute(g, Dict(index_sym(i) => 2))
g(2)
```

See also [`IndexedVariable`](@ref), [`index_name`](@ref).
"""
function index_sym(idx::Index)
    idx.name_id == 0 && return Num(Int(idx.slot))
    base = base_sym_from_id(idx.name_id)
    idx.slot == 0 && return base
    return Num(SymbolicUtils.setmetadata(SymbolicUtils.unwrap(base), IndexSlot, Int(idx.slot)))
end

find_space_index(::HilbertSpace, ::HilbertSpace) = 1
function find_space_index(h::ProductSpace, space::HilbertSpace)
    for (i, s) in enumerate(h.spaces)
        s == space && return i
    end
    throw(ArgumentError("Space $space not found in $h"))
end
