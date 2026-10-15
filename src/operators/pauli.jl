"""
    PauliSpace(name::Symbol) <: HilbertSpace

Hilbert space of a two-level system written with Pauli matrices, home of the
[`Pauli`](@ref) operators.

# Examples

```jldoctest
julia> h = PauliSpace(:qubit)
ℋ(qubit)

julia> σz = Pauli(h, :σ, 3);

julia> σz * σz
1
```

See also [`Pauli`](@ref), [`NLevelSpace`](@ref), [`SpinSpace`](@ref).
"""
struct PauliSpace <: HilbertSpace
    name::Symbol
end
PauliSpace(name::AbstractString) = name_must_be_symbol(name)
Base.:(==)(a::PauliSpace, b::PauliSpace) = a.name == b.name
Base.hash(a::PauliSpace, h::UInt) = hash(:PauliSpace, hash(a.name, h))

"""
    Pauli(h::PauliSpace, name::Symbol, axis::Int)
    Pauli(h::ProductSpace, name::Symbol, axis::Int, [subsystem::Int])

Pauli matrix along `axis`, where `1`, `2`, `3` stand for x, y, z.

Pauli operators are Hermitian and multiply as
``\\sigma_j \\sigma_k = \\delta_{jk} + i \\epsilon_{jkl} \\sigma_l``, so every product
on one qubit reduces to a single Pauli matrix or a number. On a
[`ProductSpace`](@ref), `subsystem` selects the [`PauliSpace`](@ref) factor and may
be omitted when there is only one.

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
    return Op(OP_PAULI, name_id(name), si, idx, axis, 0, 0, 0)
end
Pauli(name::Symbol, axis::Int, si::Int) = Pauli(name, axis, si, NO_INDEX)

Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
function Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range"))
    h.spaces[idx] isa PauliSpace || throw(ArgumentError("Space at index $idx is not a PauliSpace"))
    return Pauli(name, axis, idx)
end

Pauli(h::ProductSpace, name::Symbol, axis::Int) =
    Pauli(h, name, axis, unique_subspace_index(h, PauliSpace))

Pauli(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)

const levi_civita = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))
