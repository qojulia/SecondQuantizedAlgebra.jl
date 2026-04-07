"""
    ClusterSpace(original_space::HilbertSpace, N, order::Int)

Hilbert space representing `N` identical copies of `original_space`, with
correlations tracked up to the specified `order`.

Used for mean-field and cluster expansions: operators on a `ClusterSpace` can be
expanded into `order` distinct copies via [`cluster_expand`](@ref).

# Arguments
- `original_space` — the single-site Hilbert space being replicated
- `N` — total number of copies (`Int` or symbolic `Num`, used for scaling prefactors)
- `order` — number of distinct copies to track (must be ≥ 1)

# Examples
```julia
h_single = FockSpace(:site)
h_cluster = ClusterSpace(h_single, 100, 2)  # 100 copies, track 2-body correlations
h = h_cluster ⊗ FockSpace(:cavity)
```

See also [`cluster_expand`](@ref), [`has_cluster`](@ref).
"""
struct ClusterSpace{H <: HilbertSpace, T} <: HilbertSpace
    original_space::H
    N::T
    order::Int
    function ClusterSpace(original_space::H, N::T, order::Int) where {H <: HilbertSpace, T}
        order >= 1 || throw(ArgumentError("Order must be >= 1, got $order"))
        return new{H, T}(original_space, N, order)
    end
end

Base.:(==)(a::ClusterSpace, b::ClusterSpace) = a.original_space == b.original_space && a.N == b.N && a.order == b.order
Base.hash(a::ClusterSpace, h::UInt) = hash(:ClusterSpace, hash(a.original_space, hash(a.N, hash(a.order, h))))

# ClusterSpace specialization (base method in hilbertspace.jl)
_unwrap_space(h::ClusterSpace) = h.original_space

"""
    has_cluster(h::HilbertSpace) -> Bool

Return `true` if `h` is a [`ClusterSpace`](@ref) or a [`ProductSpace`](@ref) containing one.
"""
has_cluster(::HilbertSpace) = false
has_cluster(::ClusterSpace) = true
function has_cluster(h::ProductSpace)
    for space in h.spaces
        space isa ClusterSpace && return true
    end
    return false
end

"""
    cluster_expand(op::QSym, order::Int) -> Vector{QSym}
    cluster_expand(op::QSym, h::ProductSpace) -> Vector{QSym}

Create `order` copies of operator `op`, each with a distinct `copy_index`
(`1, 2, ..., order`) and name suffixed `_1`, `_2`, etc.

When called with a [`ProductSpace`](@ref), the order is read from the
[`ClusterSpace`](@ref) at `op.space_index`.

# Examples
```julia
h = FockSpace(:site)
a = Destroy(h, :a)
copies = cluster_expand(a, 3)  # [a_1, a_2, a_3]
```
"""
function cluster_expand(op::QSym, order::Int)
    return [_with_copy(op, Symbol(op.name, :_, i), i) for i in 1:order]
end

function cluster_expand(op::QSym, h::ProductSpace)
    space = h.spaces[op.space_index]
    space isa ClusterSpace || throw(ArgumentError("Space at index $(op.space_index) is not a ClusterSpace"))
    return cluster_expand(op, space.order)
end

# Per-type copy constructors
_with_copy(op::Destroy, name::Symbol, ci::Int) = Destroy(name, op.space_index, ci)
_with_copy(op::Create, name::Symbol, ci::Int) = Create(name, op.space_index, ci)
_with_copy(op::Transition, name::Symbol, ci::Int) = Transition(name, op.i, op.j, op.space_index, ci)
_with_copy(op::Pauli, name::Symbol, ci::Int) = Pauli(name, op.axis, op.space_index, ci)
_with_copy(op::Spin, name::Symbol, ci::Int) = Spin(name, op.axis, op.space_index, ci)
_with_copy(op::Position, name::Symbol, ci::Int) = Position(name, op.space_index, ci)
_with_copy(op::Momentum, name::Symbol, ci::Int) = Momentum(name, op.space_index, ci)
