"""
    PauliSpace(name::Symbol)

Hilbert space for a two-level system described by Pauli operators.

Supports [`Pauli`](@ref) operators ``\\sigma_x, \\sigma_y, \\sigma_z`` with the algebra
``\\sigma_j \\sigma_k = \\delta_{jk} I + i\\epsilon_{jkl} \\sigma_l``.

# Examples
```julia
h = PauliSpace(:spin)
@qnumbers σx::Pauli(h, :σ, 1) σy::Pauli(h, :σ, 2) σz::Pauli(h, :σ, 3)
```
"""
struct PauliSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

"""
    Pauli <: QSym

Pauli operator ``\\sigma_x, \\sigma_y, \\sigma_z`` on a [`PauliSpace`](@ref).

The `axis` field selects the component: `1` = x, `2` = y, `3` = z.
Pauli operators are Hermitian (`σ' == σ`) and satisfy the product rule
``\\sigma_j \\sigma_k = \\delta_{jk} I + i\\epsilon_{jkl} \\sigma_l``.

# Construction
```julia
h = PauliSpace(:s)
σx = Pauli(h, :σ, 1)           # σx on single space
σy = Pauli(h, :σ, 2)           # σy

hp = FockSpace(:c) ⊗ PauliSpace(:s)
σz = Pauli(hp, :σ, 3, 2)       # σz on 2nd subspace of ProductSpace
```
"""
struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    index::Index
    function Pauli(name::Symbol, axis::Int, si::Int, idx::Index)
        1 <= axis <= 3 || throw(ArgumentError("Pauli axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, si, idx)
    end
end
Pauli(name::Symbol, axis::Int, si::Int) = Pauli(name, axis, si, NO_INDEX)

Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    h.spaces[idx] isa PauliSpace || throw(ArgumentError("Space at index $idx is not a PauliSpace"))
    return Pauli(name, axis, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one PauliSpace.
Pauli(h::ProductSpace, name::Symbol, axis::Int) =
    Pauli(h, name, axis, _unique_subspace_index(h, PauliSpace))

IndexedOperator(op::Pauli, i::Index) = Pauli(op.name, op.axis, op.space_index, i)

Base.adjoint(op::Pauli) = op

Base.isequal(a::Pauli, b::Pauli) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::Pauli, b::Pauli) = isequal(a, b)

Base.hash(a::Pauli, h::UInt) = hash(:Pauli, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.index, h)))))

ladder(::Pauli) = 0

# --- Operator hooks ---

# Per-site Levi-Civita lookup (3×3 antisymmetric). _levi_civita[j][k] = ε_{jk(6-j-k)}.
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))

function _site_compare(a::Pauli, b::Pauli, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    a.index == b.index && return Equal
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

_can_commute(a::Pauli, b::Pauli) = false   # always compose on same site

# σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ
function _reduce_pair(a::Pauli, b::Pauli)
    a.name == b.name || return nothing
    a.space_index == b.space_index || return nothing
    a.index == b.index || return nothing
    if a.axis == b.axis
        return _CNUM_ONE     # σⱼ² = 1
    else
        eps = _levi_civita[a.axis][b.axis]
        new = Pauli(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
        return (new, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE))
    end
end

_commute_pair(a::Pauli, b::Pauli) = error("unreachable: Pauli uses _reduce_pair")
