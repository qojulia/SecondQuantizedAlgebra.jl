"""
    SpinSpace(name::Symbol) <: HilbertSpace

Hilbert space for collective spin angular momentum. Hosts [`Spin`](@ref)
operators satisfying ``[S_j, S_k] = i \\epsilon_{jkl} S_l``. The algebra is
independent of the spin size ``S``, which enters only via the
`QuantumOpticsBase.SpinBasis` chosen for numeric evaluation.

# Examples

```jldoctest
julia> SpinSpace(:S)
ℋ(S)
```

See also [`Spin`](@ref), [`PauliSpace`](@ref).
"""
struct SpinSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::SpinSpace, b::SpinSpace) = a.name == b.name
Base.hash(a::SpinSpace, h::UInt) = hash(:SpinSpace, hash(a.name, h))

"""
    Spin(h::SpinSpace, name::Symbol, axis::Int) -> Op

Angular momentum operator ``S_x, S_y, S_z`` on a [`SpinSpace`](@ref). The
`axis` argument selects the component: `1 = x`, `2 = y`, `3 = z`. Hermitian
(`S' == S`) and satisfies ``[S_j, S_k] = i \\epsilon_{jkl} S_l`` (applied
eagerly by `*`). Returns an [`Op`](@ref) tagged `OP_SPIN` with the axis stored
in `l1`.

# Examples

```jldoctest
julia> h = SpinSpace(:S);

julia> Sx = Spin(h, :S, 1); Sy = Spin(h, :S, 2);

julia> Sy * Sx
-im * Sz + Sx * Sy
```

See also [`SpinSpace`](@ref), [`Pauli`](@ref).
"""
function Spin(name::Union{Symbol, Int32}, axis::Int, si::Int, idx::Index)
    1 <= axis <= 3 || throw(ArgumentError("Spin axis must be 1, 2, or 3, got $axis"))
    return Op(OP_SPIN, _name_id(name), si, idx, axis, 0, 0, 0)
end
Spin(name::Union{Symbol, Int32}, axis::Int, si::Int) = Spin(name, axis, si, NO_INDEX)

Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
function Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    h.spaces[idx] isa SpinSpace || throw(ArgumentError("Space at index $idx is not a SpinSpace"))
    return Spin(name, axis, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one SpinSpace.
Spin(h::ProductSpace, name::Symbol, axis::Int) =
    Spin(h, name, axis, _unique_subspace_index(h, SpinSpace))
