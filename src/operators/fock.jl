"""
    Destroy(h::FockSpace, name::Symbol) -> Op

Bosonic annihilation operator ``a`` on a [`FockSpace`](@ref). Satisfies the
canonical commutation relation ``[a, a^\\dagger] = 1``. The adjoint `a'`
returns the corresponding [`Create`](@ref) operator. Returns an [`Op`](@ref)
tagged `OP_DESTROY`.

# Examples

```jldoctest
julia> h = FockSpace(:cavity);

julia> @qnumbers a::Destroy(h);

julia> a * a'
1 + a' * a
```

See also [`Create`](@ref), [`FockSpace`](@ref), [`@qnumbers`](@ref).
"""
Destroy(name::Symbol, si::Int, idx::Index) = Op(OP_DESTROY, _name_id(name), si, idx, 0, 0, 0, 0)
Destroy(name::Symbol, si::Int) = Destroy(name, si, NO_INDEX)

"""
    Create(h::FockSpace, name::Symbol) -> Op

Bosonic creation operator ``a^\\dagger`` on a [`FockSpace`](@ref). The
adjoint of [`Destroy`](@ref); typically obtained via `a'` rather than
constructed directly. Returns an [`Op`](@ref) tagged `OP_CREATE`.

See also [`Destroy`](@ref), [`FockSpace`](@ref).
"""
Create(name::Symbol, si::Int, idx::Index) = Op(OP_CREATE, _name_id(name), si, idx, 0, 0, 0, 0)
Create(name::Symbol, si::Int) = Create(name, si, NO_INDEX)

Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Destroy(name, idx)
end
function Create(h::ProductSpace, name::Symbol, idx::Int)
    1 <= idx <= length(h.spaces) || throw(ArgumentError("Index $idx out of range for ProductSpace with $(length(h.spaces)) spaces"))
    h.spaces[idx] isa FockSpace || throw(ArgumentError("Space at index $idx is not a FockSpace"))
    return Create(name, idx)
end

# Auto-detect subspace when the ProductSpace contains exactly one FockSpace.
Destroy(h::ProductSpace, name::Symbol) = Destroy(h, name, _unique_subspace_index(h, FockSpace))
Create(h::ProductSpace, name::Symbol) = Create(h, name, _unique_subspace_index(h, FockSpace))

Destroy(::HilbertSpace, name::AbstractString, args...) = _name_must_be_symbol(name)
Create(::HilbertSpace, name::AbstractString, args...) = _name_must_be_symbol(name)

"""
    IndexedOperator(op::Op, i::Index) -> Op

Return the indexed version of an operator.

`IndexedOperator` keeps the operator kind and quantum numbers, and only changes
the symbolic index label. Use it to build objects such as `a_i`, `σ_j₁₂`, or
`x_k` before forming sums with [`Σ`](@ref).

Use `IndexedOperator(op, NO_INDEX)` to remove an index.

# Examples
```jldoctest
julia> h = FockSpace(:cavity);

julia> @qnumbers a::Destroy(h);

julia> i = Index(h, :i, 5, h);

julia> IndexedOperator(a, i)
a_i
```

See also [`Index`](@ref), [`Σ`](@ref), [`change_index`](@ref).
"""
function IndexedOperator end
function IndexedOperator(op::Op, i::Index)
    is_collective_transition(op) && has_index(i) &&
        throw(ArgumentError("CollectiveTransition is already collective and cannot be indexed"))
    return Op(op.kind, op.name_id, op.space_index, i, op.l1, op.l2, op.g, op.nlev)
end
