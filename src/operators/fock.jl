"""
    Destroy <: QSym

Bosonic annihilation operator ``a`` on a [`FockSpace`](@ref).

Satisfies the canonical commutation relation ``[a, a^\\dagger] = 1``.
The adjoint `a'` returns the corresponding [`Create`](@ref) operator.

# Construction
```julia
h = FockSpace(:cavity)
a = Destroy(h, :a)              # single-space
h2 = FockSpace(:a) ⊗ FockSpace(:b)
a = Destroy(h2, :a, 1)          # on first subspace of ProductSpace
```
Or via the [`@qnumbers`](@ref) macro:
```julia
@qnumbers a::Destroy(h)
```
"""
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Destroy(name::Symbol, si::Int) = Destroy(name, si, NO_INDEX)

"""
    Create <: QSym

Bosonic creation operator ``a^\\dagger`` on a [`FockSpace`](@ref).

The adjoint of [`Destroy`](@ref). Satisfies `[a, a'] = 1` (applied eagerly by `*`).
Constructed implicitly via `adjoint(::Destroy)`, or directly with the same signatures
as [`Destroy`](@ref).
"""
struct Create <: QSym
    name::Symbol
    space_index::Int
    index::Index
end
Create(name::Symbol, si::Int) = Create(name, si, NO_INDEX)

Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Create(name, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one FockSpace.
Destroy(h::ProductSpace, name::Symbol) = Destroy(h, name, _unique_subspace_index(h, FockSpace))
Create(h::ProductSpace, name::Symbol) = Create(h, name, _unique_subspace_index(h, FockSpace))

"""
    IndexedOperator(op::QSym, i::Index) -> QSym

Attach a symbolic summation index to an operator.

Returns a new operator of the same type with `index = i`. To clear the index, use
`IndexedOperator(op, NO_INDEX)`.

# Examples
```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
i = Index(h, :i, 5, h)
a_i = IndexedOperator(a, i)   # symbolic: a_i
```
"""
function IndexedOperator end
IndexedOperator(op::Destroy, i::Index) = Destroy(op.name, op.space_index, i)
IndexedOperator(op::Create, i::Index) = Create(op.name, op.space_index, i)

Base.adjoint(op::Destroy) = Create(op.name, op.space_index, op.index)
Base.adjoint(op::Create) = Destroy(op.name, op.space_index, op.index)

Base.isequal(a::Destroy, b::Destroy) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.isequal(a::Create, b::Create) = a.name == b.name && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::Destroy, b::Destroy) = isequal(a, b)
Base.:(==)(a::Create, b::Create) = isequal(a, b)

Base.hash(a::Destroy, h::UInt) = hash(:Destroy, hash(a.name, hash(a.space_index, hash(a.index, h))))
Base.hash(a::Create, h::UInt) = hash(:Create, hash(a.name, hash(a.space_index, hash(a.index, h))))

"""
    ladder(op::QSym)

Returns 0 for creation operators, 1 for annihilation operators.
Used for canonical ordering within operator product sequences.
"""
ladder(::Create) = 0
ladder(::Destroy) = 1

# --- Operator hooks ---

# Same-type same-site site relationship is always Equal (same name/space/index).
# Different name or different space -> distinct (sorted by name lex order, then space).
function _site_compare(a::Destroy, b::Destroy, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end
function _site_compare(a::Create, b::Create, ne::Vector{NonEqualPair})::SiteCmp
    return _site_compare(
        Destroy(a.name, a.space_index, a.index),
        Destroy(b.name, b.space_index, b.index), ne
    )
end

# Cross-type same-site returns Equal; canonical direction lives in _can_commute.
# Distinct-site uses the Destroy/Destroy comparison (lex name, then space).
function _site_compare(a::Create, b::Destroy, ne::Vector{NonEqualPair})::SiteCmp
    return _site_compare(
        Destroy(a.name, a.space_index, a.index),
        Destroy(b.name, b.space_index, b.index), ne
    )
end
function _site_compare(a::Destroy, b::Create, ne::Vector{NonEqualPair})::SiteCmp
    return _site_compare(
        Destroy(a.name, a.space_index, a.index),
        Destroy(b.name, b.space_index, b.index), ne
    )
end

# Same-site commutation: a·a† carries a residual (the identity term);
# a†·a is already in canonical order and commutes freely.
_can_commute(a::Destroy, b::Create) = false
_can_commute(a::Create, b::Destroy) = true
_can_commute(a::Destroy, b::Destroy) = true
_can_commute(a::Create, b::Create) = true

# _commute_pair returns (swap_b, swap_a, residual_coeff, residual_ops):
# aa† = a†a + 1; residual op vector is empty (identity branch).
_commute_pair(a::Destroy, b::Create) = (b, a, _CNUM_ONE, _EMPTY_OPS)

# Reductions: ladder operators on the same site don't reduce locally.
# (a·a is not a closed-form simplification; only a·a† triggers the commutation residual.)
_reduce_pair(a::Destroy, b::Create) = nothing
_reduce_pair(a::Create, b::Destroy) = nothing
_reduce_pair(a::Destroy, b::Destroy) = nothing
_reduce_pair(a::Create, b::Create) = nothing
