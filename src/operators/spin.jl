"""
    SpinSpace(name::Symbol)

Hilbert space for collective spin angular momentum operators.

Supports [`Spin`](@ref) operators ``S_x, S_y, S_z`` with the commutation relation
``[S_j, S_k] = i\\epsilon_{jkl} S_l``. The algebraic rules are independent of the
spin size ``S`` — that only enters via the `QuantumOpticsBase.SpinBasis` chosen
for numeric evaluation.

# Examples
```julia
SpinSpace(:S)
```
"""
struct SpinSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::SpinSpace, b::SpinSpace) = a.name == b.name
Base.hash(a::SpinSpace, h::UInt) = hash(:SpinSpace, hash(a.name, h))

"""
    Spin <: QSym

Angular momentum operator ``S_x, S_y, S_z`` on a [`SpinSpace`](@ref).

The `axis` field selects the component: `1` = x, `2` = y, `3` = z.
Spin operators are Hermitian (`S' == S`) and satisfy ``[S_j, S_k] = i\\epsilon_{jkl} S_l``
(applied eagerly during `*`).

# Construction
```julia
h = SpinSpace(:S)
Sx = Spin(h, :S, 1)             # Sx on single space
Sy = Spin(h, :S, 2)             # Sy

hp = FockSpace(:c) ⊗ SpinSpace(:S)
Sz = Spin(hp, :S, 3, 2)         # Sz on 2nd subspace
```
"""
struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    index::Index
    function Spin(name::Symbol, axis::Int, si::Int, idx::Index)
        1 <= axis <= 3 || throw(ArgumentError("Spin axis must be 1, 2, or 3, got $axis"))
        return new(name, axis, si, idx)
    end
end
Spin(name::Symbol, axis::Int, si::Int) = Spin(name, axis, si, NO_INDEX)

# Construction from Hilbert spaces
Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    h.spaces[idx] isa SpinSpace || throw(ArgumentError("Space at index $idx is not a SpinSpace"))
    return Spin(name, axis, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one SpinSpace.
Spin(h::ProductSpace, name::Symbol, axis::Int) =
    Spin(h, name, axis, _unique_subspace_index(h, SpinSpace))

# IndexedOperator convenience
IndexedOperator(op::Spin, i::Index) = Spin(op.name, op.axis, op.space_index, i)

# Adjoint — Hermitian
Base.adjoint(op::Spin) = op

# Equality
Base.isequal(a::Spin, b::Spin) = a.name == b.name && a.axis == b.axis && a.space_index == b.space_index && a.index == b.index
Base.:(==)(a::Spin, b::Spin) = isequal(a, b)

# Hashing
Base.hash(a::Spin, h::UInt) = hash(:Spin, hash(a.name, hash(a.axis, hash(a.space_index, hash(a.index, h)))))

# Ladder (not applicable)
ladder(::Spin) = 0

# --- Operator hooks ---

function _site_compare(a::Spin, b::Spin, ne::Vector{NonEqualPair})::SiteCmp
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    a.name == b.name || return a.name < b.name ? Less : Greater
    if a.index != b.index
        _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
        return Undetermined
    end
    # Same site: return Equal; axis canonical order lives in _can_commute / _commute_pair.
    return Equal
end

# Same-site spins commute only when in canonical axis order (ascending).
_can_commute(a::Spin, b::Spin) = a.axis <= b.axis

# [Sj, Sk] = iϵⱼₖₗSl. The pair (Sj·Sk) with j>k commutes to (Sk·Sj + iε·Sl).
# Returns (swap_b, swap_a, residual_coeff, residual_ops) where residual_ops is
# the single contracted spin Sl on the third axis.
function _commute_pair(a::Spin, b::Spin)
    a.name == b.name || error("unreachable")
    a.space_index == b.space_index || error("unreachable")
    a.index == b.index || error("unreachable")
    a.axis > b.axis || error("unreachable: _commute_pair called on in-order pair")
    eps = _levi_civita[a.axis][b.axis]
    contracted = Spin(a.name, 6 - a.axis - b.axis, a.space_index, a.index)
    return (b, a, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE), QSym[contracted])
end

_reduce_pair(::Spin, ::Spin) = nothing
