"""
    fundamental_operators(h::HilbertSpace; names = Symbol[])

Smallest set of operators from which every operator on `h` is built by products and adjoints.

A bosonic mode contributes its annihilation operator; a finite-level system its transitions
``\\sigma^{ij}`` with ``i < j`` and its projectors except the ground state; a
[`CollectiveNLevelSpace`](@ref) all of its collective transitions; Pauli and spin systems
their three components; and a pair of quadratures position and momentum. On a product space the sets of the factors are joined.
`names` sets the operator names, one per factor; otherwise they are generated.

# Examples
```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 3);

julia> fundamental_operators(h; names = [:b, :s])
6-element Vector{Op}:
 b
 s₁₂
 s₁₃
 s₂₂
 s₂₃
 s₃₃
```

See also [`find_operators`](@ref), [`unique_up_to_adjoint`](@ref).
"""
fundamental_operators(h::FockSpace, si::Int = 1; names::AbstractVector = Symbol[]) =
    Op[Destroy(factor_name(names, si, :a), si)]

function fundamental_operators(h::NLevelSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = factor_name(names, si, :σ)
    ops = Op[]
    for i in 1:(h.n)
        for j in i:(h.n)
            (i == j) && i == h.ground_state && continue
            push!(ops, Transition(name, i, j, si, NO_INDEX, h.ground_state, h.n))
        end
    end
    return ops
end

function fundamental_operators(h::CollectiveNLevelSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = factor_name(names, si, :S)
    return Op[CollectiveTransition(name, i, j, si) for i in 1:h.n for j in 1:h.n]
end

fundamental_operators(h::PauliSpace, si::Int = 1; names::AbstractVector = Symbol[]) =
    Op[Pauli(factor_name(names, si, :σ), axis, si) for axis in 1:3]

fundamental_operators(h::SpinSpace, si::Int = 1; names::AbstractVector = Symbol[]) =
    Op[Spin(factor_name(names, si, :S), axis, si) for axis in 1:3]

function fundamental_operators(h::PhaseSpace, si::Int = 1; names::AbstractVector = Symbol[])
    x, p = isempty(names) ? (:x, :p) : phase_space_names(names[si])
    return Op[Position(x, si), Momentum(p, si)]
end

phase_space_names(n::Union{Tuple, AbstractVector}) = (n[1], n[2])
phase_space_names(n) = throw(ArgumentError("a PhaseSpace names entry must be a (position, momentum) pair, got $n"))

function fundamental_operators(h::ProductSpace; names::AbstractVector = Symbol[])
    ops = Op[]
    for (i, space) in enumerate(h.spaces)
        append!(ops, fundamental_operators(space, i; names))
    end
    return ops
end

"""
    unique_up_to_adjoint(ops)

Operators of `ops` with duplicates removed, counting an operator and its adjoint as the same.

The first of each pair is kept. Equations for averages need only one of ``\\langle X \\rangle``
and ``\\langle X^\\dagger \\rangle``, since the other is its complex conjugate.

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> unique_up_to_adjoint([a, a', a' * a])
2-element Vector{SecondQuantizedAlgebra.QField}:
 a
 a' * a
```

See also [`find_operators`](@ref), [`fundamental_operators`](@ref).
"""
unique_up_to_adjoint(ops::AbstractVector{<:QField}) = unique_up_to_adjoint!(copy(ops))

function unique_up_to_adjoint!(ops::AbstractVector{<:QField})
    seen = Set{UInt}()
    j = 0
    for i in eachindex(ops)
        h = hash(ops[i])
        h_adj = hash(adjoint(ops[i]))
        if h ∉ seen && h_adj ∉ seen
            push!(seen, h)
            h == h_adj || push!(seen, h_adj)
            j += 1
            if j != i
                ops[j] = ops[i]
            end
        end
    end
    resize!(ops, j)
    return ops
end

"""
    find_operators(h::HilbertSpace, order::Int; names = Symbol[])

All distinct products of up to `order` operators on `h`, counting adjoints only once.

The products are built from [`fundamental_operators`](@ref) and their adjoints, brought into
canonical form, and filtered with [`unique_up_to_adjoint`](@ref). They are the averages a
set of moment or cumulant equations of that order is written in. `names` is passed on to
[`fundamental_operators`](@ref).

# Examples
```jldoctest
julia> find_operators(FockSpace(:cavity), 2)
3-element Vector{Union{Op, SecondQuantizedAlgebra.QAdd}}:
 a
 a * a
 a' * a
```

See also [`fundamental_operators`](@ref), [`unique_up_to_adjoint`](@ref).
"""
function find_operators(h::HilbertSpace, order::Int; names::AbstractVector = Symbol[])
    fund_ops = fundamental_operators(h; names = isempty(names) ? default_names(h) : names)
    fund_ops = unique(vcat(fund_ops, adjoint.(fund_ops)))

    all_ops = Union{Op, QAdd}[]
    for i in 1:order
        for c in with_replacement_combinations(fund_ops, i)
            product = prod(reverse(c))
            iszero(product) && continue
            if product isa QAdd && length(product.arguments) == 1
                (t, v), = product.arguments
                push!(all_ops, length(t.ops) == 1 && isone(v) ? only(t.ops) : product)
            elseif product isa QSym
                push!(all_ops, product)
            end
        end
    end

    unique_up_to_adjoint!(all_ops)
    return all_ops
end

default_names(::HilbertSpace) = Symbol[]
default_names(h::ProductSpace) =
    allunique(map(typeof, h.spaces)) ? Symbol[] : Symbol.(('a':'z')[1:length(h.spaces)])

factor_name(names::AbstractVector, si::Int, default::Symbol) = isempty(names) ? default : names[si]
