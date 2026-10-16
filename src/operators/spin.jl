"""
    SpinSpace(name::Symbol) <: HilbertSpace

Hilbert space of an angular momentum, home of the [`Spin`](@ref) operators.

The spin length ``S`` does not enter the algebra. It is fixed only by the basis
chosen for numerics, so the same expressions hold for every ``S``.

# Examples

```jldoctest
julia> h = SpinSpace(:S)
ℋ(S)

julia> Sz = Spin(h, :S, 3)
Sz
```

See also [`Spin`](@ref), [`PauliSpace`](@ref).
"""
struct SpinSpace <: HilbertSpace
    name::Symbol
end
SpinSpace(name::AbstractString) = name_must_be_symbol(name)

"""
    Spin(h::SpinSpace, name::Symbol, axis::Int)
    Spin(h::ProductSpace, name::Symbol, axis::Int, [subsystem::Int])

Angular momentum component along `axis`, where `1`, `2`, `3` stand for x, y, z.

Spin operators are Hermitian and obey ``[S_j, S_k] = i \\epsilon_{jkl} S_l``.
Products are ordered by axis, ``x`` before ``y`` before ``z``. Squares such as
``S_x^2`` stay as they are, since the spin length is not fixed. On a
[`ProductSpace`](@ref), `subsystem` selects the [`SpinSpace`](@ref) factor and may
be omitted when there is only one.

# Examples

```jldoctest
julia> h = SpinSpace(:S);

julia> Sx = Spin(h, :S, 1); Sy = Spin(h, :S, 2);

julia> Sy * Sx
-im * Sz + Sx * Sy
```

See also [`SpinSpace`](@ref), [`Pauli`](@ref).
"""
Spin(name::Symbol, axis::Int, si::Int) = Op(OP_SPIN, name_id(name), si, NO_INDEX, check_axis(:Spin, axis), 0, 0, 0)
Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int = unique_subspace_index(h, SpinSpace)) =
    on_factor(Spin, h, SpinSpace, idx, name, axis)

Spin(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
