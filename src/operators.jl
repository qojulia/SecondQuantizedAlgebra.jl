"""
    fundamental_operators(h::HilbertSpace; names=nothing) -> Vector{QSym}

Return the fundamental (basis) operators for a Hilbert space `h`.

The operator set depends on the space type:
- [`FockSpace`](@ref): one [`Destroy`](@ref) operator
- [`NLevelSpace`](@ref) with `n` levels: `n(n+1)/2 - 1` [`Transition`](@ref) operators
  (diagonal + upper-triangular, excluding the ground-state projector)
- [`PauliSpace`](@ref): three [`Pauli`](@ref) operators (σx, σy, σz)
- [`SpinSpace`](@ref): three [`Spin`](@ref) operators (Sx, Sy, Sz)
- [`PhaseSpace`](@ref): [`Position`](@ref) and [`Momentum`](@ref)
- [`ProductSpace`](@ref): concatenation of fundamental operators from each subspace

# Keyword arguments
- `names=nothing` — custom operator names per subspace (default: auto-generated).

# Examples
```julia
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
ops = fundamental_operators(h)   # [Destroy(:a, 1), Transition(:σ, 1, 2, 2, NO_INDEX, 1, 2)]
```

See also [`find_operators`](@ref), [`unique_ops`](@ref).
"""
function fundamental_operators(h::FockSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :a : names[si]
    return QSym[Destroy(name, si)]
end

function fundamental_operators(h::NLevelSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :σ : names[si]
    ops = Transition[]
    for i in 1:(h.n)
        for j in i:(h.n)
            (i == j) && i == h.ground_state && continue
            push!(ops, Transition(name, i, j, si, NO_INDEX, h.ground_state, h.n))
        end
    end
    return ops
end

function fundamental_operators(h::PauliSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :σ : names[si]
    return QSym[Pauli(name, 1, si), Pauli(name, 2, si), Pauli(name, 3, si)]
end

function fundamental_operators(h::SpinSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :S : names[si]
    return QSym[Spin(name, 1, si), Spin(name, 2, si), Spin(name, 3, si)]
end

function fundamental_operators(h::PhaseSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name_pair = names === nothing ? (:x, :p) : names[si]
    return QSym[Position(name_pair[1], si), Momentum(name_pair[2], si)]
end

function fundamental_operators(h::ProductSpace; names::Union{Nothing, AbstractVector} = nothing)
    ops = QSym[]
    for (i, space) in enumerate(h.spaces)
        space_ops = fundamental_operators(space, i; names = names)
        append!(ops, space_ops)
    end
    return ops
end

"""
    unique_ops(ops) -> Vector

Return unique operators from `ops`, treating `op` and `op'` (adjoint) as the
same degree of freedom. Only the first occurrence of each operator/adjoint pair
is kept.

See also [`unique_ops!`](@ref), [`fundamental_operators`](@ref).
"""
function unique_ops(ops::AbstractVector)
    ops_ = copy(ops)
    unique_ops!(ops_)
    return ops_
end

"""
    unique_ops!(ops) -> Vector

In-place version of [`unique_ops`](@ref).
"""
function unique_ops!(ops::AbstractVector)
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
    find_operators(h::HilbertSpace, order::Int; names=nothing) -> Vector

Generate all unique operator products up to the given `order` for Hilbert space `h`.

Starting from [`fundamental_operators(h)`](@ref fundamental_operators) and their adjoints,
computes all products with up to `order` factors. Filters out zero terms (from algebraic
identities) and adjoint-duplicates (via [`unique_ops!`](@ref)).

Useful for constructing the operator basis needed for cumulant expansions or
moment equations.

# Keyword arguments
- `names=nothing` — passed to [`fundamental_operators`](@ref).

# Examples
```julia
h = NLevelSpace(:atom, 2)
ops = find_operators(h, 2)   # σ₁₂, σ₂₁, σ₂₂, σ₂₂σ₁₂, ...
```

See also [`fundamental_operators`](@ref), [`unique_ops`](@ref).
"""
function find_operators(h::HilbertSpace, order::Int; names::Union{Nothing, AbstractVector} = nothing)
    # Auto-generate names when ProductSpace has duplicate space types
    if names === nothing && h isa ProductSpace
        space_types = typeof.(h.spaces)
        if length(unique(space_types)) != length(space_types)
            alph = 'a':'z'
            names = Symbol.(alph[1:length(h.spaces)])
        end
    end
    fund_ops = fundamental_operators(h; names = names)
    fund_ops = unique(vcat(fund_ops, adjoint.(fund_ops)))

    all_ops = QField[]
    for i in 1:order
        for c in with_replacement_combinations(fund_ops, i)
            c_ = prod(reverse(c))
            iszero(c_) && continue
            # With eager ordering, c_ is QAdd. Keep only single-term products.
            if c_ isa QAdd && length(c_.arguments) == 1
                push!(all_ops, c_)
            elseif c_ isa QSym
                push!(all_ops, c_)
            end
        end
    end

    unique_ops!(all_ops)
    return all_ops
end


# --- Conjugation helpers ---

"""
    _conj(v)

Recursively conjugate a symbolic expression. For `BasicSymbolic` trees,
applies `conj` to leaves and reconstructs via `maketerm`.
"""
function _conj(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(_conj, SymbolicUtils.arguments(v))
        return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
_conj(x::Number) = conj(x)
_conj(x::QField) = adjoint(x)

"""
    _inconj(v)

Conjugate an expression containing averages, preserving operator ordering inside
averages (uses `adjoint` on the inner operator rather than naive `conj`).
"""
function _inconj(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(v)[1]
            inner = SymbolicUtils.isconst(arg) ? arg.val : arg
            return _average(adjoint(inner))
        elseif f === conj
            # conj(avg(...)) → recurse into the argument
            return _inconj(SymbolicUtils.arguments(v)[1])
        else
            args = map(_inconj, SymbolicUtils.arguments(v))
            return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
        end
    else
        return conj(v)
    end
end
_inconj(x::Number) = conj(x)
_inconj(x::Num) = _inconj(SymbolicUtils.unwrap(x))
_inconj(x::QField) = adjoint(x)

"""
    _adjoint(x)

Adjoint that works on both operators (uses `adjoint`) and symbolic numbers
(uses `_conj`).
"""
_adjoint(op::QField) = adjoint(op)
_adjoint(s::SymbolicUtils.BasicSymbolic) = _conj(s)
_adjoint(x::Number) = adjoint(x)
