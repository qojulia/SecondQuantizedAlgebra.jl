"""
    HilbertSpace

Supertype of all Hilbert spaces.

The concrete spaces are [`FockSpace`](@ref), [`NLevelSpace`](@ref),
[`CollectiveNLevelSpace`](@ref), [`PauliSpace`](@ref), [`SpinSpace`](@ref),
[`PhaseSpace`](@ref) and their tensor product [`ProductSpace`](@ref). The type is public but
not exported.

# Examples
```jldoctest
julia> FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2) isa SecondQuantizedAlgebra.HilbertSpace
true
```

See also [`ProductSpace`](@ref), [`⊗`](@ref).
"""
abstract type HilbertSpace end

@noinline name_must_be_symbol(name::AbstractString) =
    throw(ArgumentError("name must be a `Symbol`, not a `String`; use `:$name` instead of `\"$name\"`"))

"""
    FockSpace(name::Symbol) <: HilbertSpace

Hilbert space of a bosonic mode, home of the ladder operators [`Destroy`](@ref) and
[`Create`](@ref). The `name` is a label used for printing.

# Examples

```jldoctest
julia> h = FockSpace(:cavity)
ℋ(cavity)

julia> a = Destroy(h, :a);

julia> a * a'
1 + a' * a
```

See also [`Destroy`](@ref), [`PhaseSpace`](@ref), [`⊗`](@ref).
"""
struct FockSpace <: HilbertSpace
    name::Symbol
end
FockSpace(name::AbstractString) = name_must_be_symbol(name)

"""
    ProductSpace{T} <: HilbertSpace

Tensor product of Hilbert spaces, built with [`⊗`](@ref) or [`tensor`](@ref).

The factors are numbered `1, 2, …` from left to right. An operator on a product
space names the factor it acts on by this number, e.g. `Destroy(h, :a, 1)`. The
number may be left out when only one factor has the right type; otherwise the
constructor throws an `ArgumentError`.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
ℋ(cavity) ⊗ ℋ(atom)

julia> a = Destroy(h, :a, 1); σ = Transition(h, :σ, 1, 2);

julia> a * σ
a * σ₁₂
```

See also [`⊗`](@ref), [`acts_on`](@ref).
"""
struct ProductSpace{T <: Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
Base.:(==)(a::ProductSpace, b::ProductSpace) = a.spaces == b.spaces
Base.hash(a::ProductSpace, h::UInt) = hash(:ProductSpace, hash(a.spaces, h))

"""
    ⊗(space::HilbertSpace, spaces::HilbertSpace...)

Tensor product of Hilbert spaces, returning a [`ProductSpace`](@ref).

Nested products are flattened, so `(A ⊗ B) ⊗ C == A ⊗ (B ⊗ C)` and the factors
are always numbered left to right. Type it as `\\otimes<tab>`, or use the ASCII
alias [`tensor`](@ref).

# Examples

```jldoctest
julia> A = FockSpace(:a); B = FockSpace(:b); C = NLevelSpace(:atom, 2);

julia> A ⊗ B ⊗ C
ℋ(a) ⊗ ℋ(b) ⊗ ℋ(atom)

julia> (A ⊗ B) ⊗ C == A ⊗ (B ⊗ C)
true
```

See also [`ProductSpace`](@ref), [`tensor`](@ref).
"""
⊗(a::HilbertSpace, b::HilbertSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::HilbertSpace) = ProductSpace((a.spaces..., b))
⊗(a::HilbertSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
⊗(a::HilbertSpace) = a

@doc """
    tensor(space::HilbertSpace, spaces::HilbertSpace...)

ASCII alias for [`⊗`](@ref), the tensor product of Hilbert spaces.

# Examples

```jldoctest
julia> tensor(FockSpace(:a), FockSpace(:b))
ℋ(a) ⊗ ℋ(b)
```

See also [`⊗`](@ref), [`ProductSpace`](@ref).
""" tensor

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
    unique_subspace_index(h::ProductSpace, ::Type{T}) -> Int

Index of the only factor of `h` that is a `T`. Throws an `ArgumentError` unless exactly one
factor is a `T`.
"""
function unique_subspace_index(h::ProductSpace, ::Type{T}) where {T <: HilbertSpace}
    n = count(s -> s isa T, h.spaces)
    n == 0 && throw(ArgumentError("No $T found in ProductSpace; specify the subspace index explicitly"))
    n > 1 && throw(ArgumentError("$n $T subspaces found in ProductSpace; specify which one with the subspace index argument"))
    return position_of(s -> s isa T, h.spaces)
end

"""
    subspace(h::ProductSpace, idx::Int, ::Type{T}) -> T

Factor `idx` of `h`, checked to exist and to be a `T`. Throws an `ArgumentError` otherwise.
"""
function subspace(h::ProductSpace, idx::Int, ::Type{T}) where {T <: HilbertSpace}
    n = length(h.spaces)
    1 <= idx <= n || throw(ArgumentError("Index $idx out of range for ProductSpace with $n spaces"))
    space = h.spaces[idx]
    space isa T || throw(ArgumentError("Space at index $idx is a $(nameof(typeof(space))), expected $(nameof(T))"))
    return space
end

"""
    on_factor(constructor, h::ProductSpace, ::Type{T}, idx::Int, name::Symbol, args...)

Operator built by `constructor` on factor `idx` of `h`, which must be a `T`.

The operator is built and validated on the single space, then moved to factor `idx`, so
every role constructor gets its `ProductSpace` methods from its single-space one.
"""
on_factor(constructor::F, h::ProductSpace, ::Type{T}, idx::Int, name::Symbol, args...) where {F, T <: HilbertSpace} =
    set_acts_on(constructor(subspace(h, idx, T), name, args...), idx)
