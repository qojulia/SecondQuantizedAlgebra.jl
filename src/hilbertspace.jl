"""
    HilbertSpace

Abstract type for representing Hilbert spaces. Used at construction time
for operator validation — not stored on operators at runtime.
"""
abstract type HilbertSpace end

"""
    FockSpace <: HilbertSpace

Hilbert space for bosonic operators (quantum harmonic oscillator).
"""
struct FockSpace <: HilbertSpace
    name::Symbol
end
Base.:(==)(a::FockSpace, b::FockSpace) = a.name == b.name
Base.hash(a::FockSpace, h::UInt) = hash(:FockSpace, hash(a.name, h))

"""
    ProductSpace{T} <: HilbertSpace

Composite Hilbert space consisting of multiple subspaces.
Uses a Tuple for fully concrete storage.
"""
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
Base.:(==)(a::ProductSpace, b::ProductSpace) = a.spaces == b.spaces
Base.hash(a::ProductSpace, h::UInt) = hash(:ProductSpace, hash(a.spaces, h))

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) consisting of multiple subspaces.
Unicode `\\otimes<tab>`.
"""
⊗(a::HilbertSpace, b::HilbertSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::HilbertSpace) = ProductSpace((a.spaces..., b))
⊗(a::HilbertSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a

"""
    tensor(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref). Alias for [`⊗`](@ref).
"""
tensor(args::Vararg{HilbertSpace}) = ⊗(args...)

Base.isless(h1::HilbertSpace, h2::HilbertSpace) = isless(h1.name, h2.name)
Base.isless(h1::ProductSpace, h2::ProductSpace) = isless(h1.spaces, h2.spaces)

# Unwrap ClusterSpace to get the original space for validation.
# Base method returns the space as-is; ClusterSpace overrides in cluster.jl.
_unwrap_space(h::HilbertSpace) = h
