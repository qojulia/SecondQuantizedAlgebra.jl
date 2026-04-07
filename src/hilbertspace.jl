"""
    HilbertSpace

Abstract supertype for Hilbert spaces.

Concrete subtypes: [`FockSpace`](@ref), [`NLevelSpace`](@ref), [`PauliSpace`](@ref),
[`SpinSpace`](@ref), [`PhaseSpace`](@ref), [`ProductSpace`](@ref), [`ClusterSpace`](@ref).

Compose spaces with [`⊗`](@ref) (or [`tensor`](@ref)):
```julia
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
```
"""
abstract type HilbertSpace end

"""
    FockSpace <: HilbertSpace
    FockSpace(name::Symbol)
    FockSpace(; name::Symbol)

Hilbert space for a bosonic mode.

Supports [`Destroy`](@ref) and [`Create`](@ref) operators with the canonical
commutation relation ``[a, a^\\dagger] = 1``.

# Examples
```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
a' * a  # number operator in normal order
```
"""
struct FockSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::FockSpace, b::FockSpace) = a.name == b.name
Base.hash(a::FockSpace, h::UInt) = hash(:FockSpace, hash(a.name, h))

"""
    ProductSpace{T} <: HilbertSpace

Composite Hilbert space formed by the tensor product of multiple subspaces.
The type parameter `T` is a concrete `Tuple` type for type-stable storage.

Constructed via [`⊗`](@ref) or [`tensor`](@ref), not directly:
```julia
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
```

Operators on a `ProductSpace` require a positional index specifying which subspace
they act on:
```julia
@qnumbers a::Destroy(h, 1) σ::Transition(h, :σ, 1, 2, 2)
```
"""
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
Base.:(==)(a::ProductSpace, b::ProductSpace) = a.spaces == b.spaces
Base.hash(a::ProductSpace, h::UInt) = hash(:ProductSpace, hash(a.spaces, h))

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) from multiple Hilbert spaces. Flattens nested
`ProductSpace` arguments so that `(A ⊗ B) ⊗ C == A ⊗ B ⊗ C`.

Unicode input: `\\otimes<tab>`. ASCII alias: [`tensor`](@ref).

# Examples
```julia
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
```
"""
⊗(a::HilbertSpace, b::HilbertSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::HilbertSpace) = ProductSpace((a.spaces..., b))
⊗(a::HilbertSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a

"""
    tensor(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) from multiple Hilbert spaces.
ASCII alias for [`⊗`](@ref).

# Examples
```julia
h = tensor(FockSpace(:a), FockSpace(:b))  # equivalent to FockSpace(:a) ⊗ FockSpace(:b)
```
"""
tensor(args::Vararg{HilbertSpace}) = ⊗(args...)

Base.isless(h1::HilbertSpace, h2::HilbertSpace) = isless(h1.name, h2.name)
Base.isless(h1::ProductSpace, h2::ProductSpace) = isless(h1.spaces, h2.spaces)

# Unwrap ClusterSpace to get the original space for validation.
# Base method returns the space as-is; ClusterSpace overrides in cluster.jl.
_unwrap_space(h::HilbertSpace) = h
