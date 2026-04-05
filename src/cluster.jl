"""
    ClusterSpace{H,T} <: HilbertSpace

Hilbert space representing N identical copies of another space, with
correlations tracked up to a specified order.

Fields:
- `original_space::H` — the space being replicated
- `N::T` — total number of copies (Int or symbolic, used by QC for scaling)
- `order::Int` — number of distinct copies to track in the algebra
"""
struct ClusterSpace{H<:HilbertSpace,T} <: HilbertSpace
    original_space::H
    N::T
    order::Int
    function ClusterSpace(original_space::H, N::T, order::Int) where {H<:HilbertSpace,T}
        order >= 1 || throw(ArgumentError("Order must be >= 1, got $order"))
        return new{H,T}(original_space, N, order)
    end
end

Base.:(==)(a::ClusterSpace, b::ClusterSpace) = a.original_space == b.original_space && a.N == b.N && a.order == b.order
Base.hash(a::ClusterSpace, h::UInt) = hash(:ClusterSpace, hash(a.original_space, hash(a.N, hash(a.order, h))))

# ClusterSpace specialization (base method in hilbertspace.jl)
_unwrap_space(h::ClusterSpace) = h.original_space

"""
    has_cluster(h::HilbertSpace) -> Bool

Check if a Hilbert space contains any `ClusterSpace`.
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

Create `order` copies of `op` with `copy_index` set to `1, 2, ..., order`
and name suffixed `_1`, `_2`, etc.
"""
function cluster_expand(op::QSym, order::Int)
    return [_with_copy(op, Symbol(op.name, :_, i), i) for i in 1:order]
end

"""
    cluster_expand(op::QSym, h::ProductSpace) -> Vector{QSym}

Create copies of `op` using the order from the `ClusterSpace` at `op.space_index`.
"""
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
