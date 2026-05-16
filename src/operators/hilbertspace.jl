"""
    HilbertSpace

Abstract supertype for Hilbert spaces. Concrete subtypes:
[`FockSpace`](@ref), [`NLevelSpace`](@ref), [`PauliSpace`](@ref),
[`SpinSpace`](@ref), [`PhaseSpace`](@ref), and [`ProductSpace`](@ref) for
tensor products. Compose with [`⊗`](@ref) or [`tensor`](@ref).
"""
abstract type HilbertSpace end

"""
    FockSpace(name::Symbol) <: HilbertSpace

Hilbert space for a bosonic mode. Hosts [`Destroy`](@ref) and [`Create`](@ref)
operators satisfying ``[a, a^\\dagger] = 1``.

# Examples

```jldoctest
julia> FockSpace(:cavity)
ℋ(cavity)
```

See also [`Destroy`](@ref), [`Create`](@ref), [`⊗`](@ref).
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
Constructed via [`⊗`](@ref) or [`tensor`](@ref), not directly.

Operators on a `ProductSpace` take a positional index specifying which
subspace they act on, e.g. `Destroy(h, :a, 1)`.

# Examples

```jldoctest
julia> FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
ℋ(cavity) ⊗ ℋ(atom)
```

See also [`⊗`](@ref), [`tensor`](@ref).
"""
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
Base.:(==)(a::ProductSpace, b::ProductSpace) = a.spaces == b.spaces
Base.hash(a::ProductSpace, h::UInt) = hash(:ProductSpace, hash(a.spaces, h))

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) from multiple Hilbert spaces. Flattens nested
`ProductSpace` arguments so that `(A ⊗ B) ⊗ C == A ⊗ B ⊗ C`. Unicode input
`\\otimes<tab>`; ASCII alias [`tensor`](@ref).

# Examples

```jldoctest
julia> FockSpace(:a) ⊗ FockSpace(:b) ⊗ NLevelSpace(:atom, 2)
ℋ(a) ⊗ ℋ(b) ⊗ ℋ(atom)
```

See also [`ProductSpace`](@ref), [`tensor`](@ref).
"""
⊗(a::HilbertSpace, b::HilbertSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::HilbertSpace) = ProductSpace((a.spaces..., b))
⊗(a::HilbertSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a

"""
    tensor(spaces::HilbertSpace...)

ASCII alias for [`⊗`](@ref): create a [`ProductSpace`](@ref) from multiple
Hilbert spaces.

# Examples

```jldoctest
julia> tensor(FockSpace(:a), FockSpace(:b))
ℋ(a) ⊗ ℋ(b)
```

See also [`⊗`](@ref), [`ProductSpace`](@ref).
"""
tensor(args::Vararg{HilbertSpace}) = ⊗(args...)

Base.isless(h1::HilbertSpace, h2::HilbertSpace) = isless(h1.name, h2.name)
Base.isless(h1::ProductSpace, h2::ProductSpace) = isless(h1.spaces, h2.spaces)

"""
    length(h::HilbertSpace) -> Int

Number of subspaces in `h`. Single Hilbert spaces ([`FockSpace`](@ref),
[`NLevelSpace`](@ref), etc.) report `1`; a [`ProductSpace`](@ref) reports the
number of factor spaces it contains.
"""
Base.length(::HilbertSpace) = 1
Base.length(h::ProductSpace) = length(h.spaces)

"""
    _unique_subspace_index(h::ProductSpace, ::Type{T}) -> Int

Return the index of the unique subspace of type `T` inside `h`. Throws
`ArgumentError` if there are zero or more than one such subspaces — the caller
must then specify the subspace index explicitly.
"""
function _unique_subspace_index(h::ProductSpace, ::Type{T}) where {T <: HilbertSpace}
    matches = findall(s -> s isa T, collect(h.spaces))
    if isempty(matches)
        throw(ArgumentError("No $T found in ProductSpace; specify the subspace index explicitly"))
    elseif length(matches) > 1
        throw(ArgumentError("$(length(matches)) $T subspaces found in ProductSpace; specify which one with the subspace index argument"))
    end
    return only(matches)
end
