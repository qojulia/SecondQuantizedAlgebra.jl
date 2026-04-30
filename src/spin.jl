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
(applied eagerly under [`NormalOrder`](@ref)).

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
