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
Destroy(name::Symbol, si::Int) = Op(OP_DESTROY, name_id(name), si, NO_INDEX, 0, 0, 0, 0)

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
Create(name::Symbol, si::Int) = Op(OP_CREATE, name_id(name), si, NO_INDEX, 0, 0, 0, 0)

Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

Destroy(h::ProductSpace, name::Symbol, idx::Int = unique_subspace_index(h, FockSpace)) =
    on_factor(Destroy, h, FockSpace, idx, name)
Create(h::ProductSpace, name::Symbol, idx::Int = unique_subspace_index(h, FockSpace)) =
    on_factor(Create, h, FockSpace, idx, name)

Destroy(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
Create(::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)
