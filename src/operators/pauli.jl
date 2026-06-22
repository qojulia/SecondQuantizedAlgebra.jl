"""
    PauliSpace(name::Symbol) <: HilbertSpace

Hilbert space for a two-level system. Hosts [`Pauli`](@ref) operators
satisfying ``\\sigma_j \\sigma_k = \\delta_{jk} I + i \\epsilon_{jkl} \\sigma_l``.

# Examples

```jldoctest
julia> PauliSpace(:spin)
ℋ(spin)
```

See also [`Pauli`](@ref), [`SpinSpace`](@ref).
"""
struct PauliSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

"""
    Pauli(h::PauliSpace, name::Symbol, axis::Int) -> Op

Pauli operator ``\\sigma_x, \\sigma_y, \\sigma_z`` on a [`PauliSpace`](@ref).
The `axis` argument selects the component: `1 = x`, `2 = y`, `3 = z`. Hermitian
(`σ' == σ`) and satisfies
``\\sigma_j \\sigma_k = \\delta_{jk} I + i \\epsilon_{jkl} \\sigma_l``. Returns an
[`Op`](@ref) tagged `OP_PAULI` with the axis stored in `l1`.

# Examples

```jldoctest
julia> h = PauliSpace(:s);

julia> σx = Pauli(h, :σ, 1); σy = Pauli(h, :σ, 2);

julia> σx * σy
im * σz

julia> σx * σx
1
```

See also [`PauliSpace`](@ref), [`Spin`](@ref).
"""
function Pauli(name::Symbol, axis::Int, si::Int, idx::Index)
    1 <= axis <= 3 || throw(ArgumentError("Pauli axis must be 1, 2, or 3, got $axis"))
    return Op(OP_PAULI, name, si, idx, axis, 0, 0, 0)
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

# Per-site Levi-Civita lookup (3×3 antisymmetric); shared by Pauli and Spin
# hooks in `operators.jl`. _levi_civita[j][k] = ε_{jk(6-j-k)}.
const _levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))
