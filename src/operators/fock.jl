"""
    Destroy(h::FockSpace, name::Symbol)
    Destroy(h::ProductSpace, name::Symbol, [subsystem::Int])

Annihilation operator of the bosonic mode `name`.

The operator ``a`` obeys ``[a, a^\\dagger] = 1``, and products come out normal
ordered. The adjoint `a'` is the [`Create`](@ref) operator. Modes with different
names are independent. On a [`ProductSpace`](@ref), `subsystem` selects the
[`FockSpace`](@ref) factor and may be omitted when there is only one.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> a = Destroy(h, :a);

julia> a * a'
1 + a' * a
```

See also [`Create`](@ref), [`FockSpace`](@ref), [`@qnumbers`](@ref).
"""
Destroy(name::Symbol, si::Int, idx::Index) = Op(OP_DESTROY, name_id(name), si, idx, 0, 0, 0, 0)
Destroy(name::Symbol, si::Int) = Destroy(name, si, NO_INDEX)

"""
    Create(h::FockSpace, name::Symbol)
    Create(h::ProductSpace, name::Symbol, [subsystem::Int])

Creation operator of the bosonic mode `name`.

The operator ``a^\\dagger`` is the adjoint of [`Destroy`](@ref), and usually written
as `a'` instead. Arguments are as for [`Destroy`](@ref).

# Examples

```jldoctest
julia> h = FockSpace(:cavity);

julia> a = Destroy(h, :a);

julia> Create(h, :a) == a'
true
```

See also [`Destroy`](@ref), [`FockSpace`](@ref).
"""
Create(name::Symbol, si::Int, idx::Index) = Op(OP_CREATE, name_id(name), si, idx, 0, 0, 0, 0)
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

Destroy(h::ProductSpace, name::Symbol) = Destroy(h, name, unique_subspace_index(h, FockSpace))
Create(h::ProductSpace, name::Symbol) = Create(h, name, unique_subspace_index(h, FockSpace))

Destroy(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
Create(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)

"""
    IndexedOperator(op, i::Index)

Copy of the operator `op` placed on site `i`.

Operators on different sites act independently; operators on the same site follow the
algebra of `op`. Any operator except a [`CollectiveTransition`](@ref) can be indexed.

# Examples
```jldoctest
julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, 10, h);

julia> σ = Transition(h, :σ, 2, 1);

julia> IndexedOperator(σ, i) * IndexedOperator(σ', i)
σ_i₂₂
```

See also [`Index`](@ref), [`Σ`](@ref), [`operator_index`](@ref).
"""
function IndexedOperator end
function IndexedOperator(op::Op, i::Index)
    is_collective_transition(op) && has_index(i) &&
        throw(ArgumentError("CollectiveTransition is already collective and cannot be indexed"))
    return Op(op.kind, op.name_id, op.space_index, i, op.l1, op.l2, op.g, op.nlev)
end
