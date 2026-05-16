"""
    fundamental_operators(h::HilbertSpace; names=nothing) -> Vector{QSym}

Return the minimal generating set of operators for Hilbert space `h`.

Returns one [`Destroy`](@ref) per [`FockSpace`](@ref); `n(n+1)/2 - 1`
[`Transition`](@ref) operators per [`NLevelSpace`](@ref) (upper-triangular plus
diagonals, excluding the ground-state projector); three [`Pauli`](@ref) or
[`Spin`](@ref) operators per [`PauliSpace`](@ref)/[`SpinSpace`](@ref); one
[`Position`](@ref) and one [`Momentum`](@ref) per [`PhaseSpace`](@ref); and the
concatenation of the above for a [`ProductSpace`](@ref). Pass `names` to override
the auto-generated operator names.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> length(fundamental_operators(h))
3
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

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> length(unique_ops([a, a', a]))
1
```

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

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> @qnumbers a::Destroy(h);

julia> v = [a, a'];

julia> SecondQuantizedAlgebra.unique_ops!(v);

julia> length(v)
1
```
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

Generate all unique operator products up to `order` factors for Hilbert space `h`.

Starts from [`fundamental_operators(h)`](@ref fundamental_operators) and their
adjoints, forms all products with up to `order` factors, then filters out zero
terms and adjoint-duplicates. Useful for constructing the operator basis needed
for cumulant expansions or moment equations. Pass `names` to override
auto-generated operator names.

# Examples

```jldoctest
julia> h = FockSpace(:f);

julia> length(find_operators(h, 1))
1
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


# --- Hermitian conjugation ---

"""
    qadjoint(x)

Hermitian conjugate that distributes through mixed operator/symbolic expressions.

On a [`QField`](@ref) returns `adjoint(x)`. On a `Number` returns `adjoint(x)`.
On a `SymbolicUtils.BasicSymbolic` tree, recurses into arguments so coefficients
distribute (e.g. `qadjoint(2im * a)` becomes `-2im * a'`) rather than producing
an opaque `conj(...)` wrapper.

Aliased as `qconj` and `dagger`.
"""
function qadjoint(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = map(qadjoint, SymbolicUtils.arguments(v))
        return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
qadjoint(x::Number) = adjoint(x)
qadjoint(x::QField) = adjoint(x)

const qconj = qadjoint
const dagger = qadjoint

"""
    inner_adjoint(x)

Push the adjoint inside `⟨...⟩` averages: rewrites `conj(⟨X⟩)` as `⟨X†⟩` so the
result stays expressed as an average of an operator rather than a `conj` wrapper
around one. Used when building equations of motion where both sides must share
the canonical "average-of-operator" form for substitution and hashing.

On non-average sub-expressions, behaves like [`qadjoint`](@ref).
"""
function inner_adjoint(v::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        if f isa AvgFunc
            arg = SymbolicUtils.arguments(v)[1]
            inner = SymbolicUtils.isconst(arg) ? arg.val : arg
            return _average(adjoint(inner))
        elseif f === conj
            return inner_adjoint(SymbolicUtils.arguments(v)[1])
        else
            args = map(inner_adjoint, SymbolicUtils.arguments(v))
            return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
        end
    else
        return conj(v)
    end
end
inner_adjoint(x::Number) = conj(x)
inner_adjoint(x::Num) = inner_adjoint(SymbolicUtils.unwrap(x))
inner_adjoint(x::QField) = adjoint(x)

# Total ordering across QSym concrete types — used by _site_compare cross-type fallback.
_type_order(::Type{Destroy}) = 0
_type_order(::Type{Create}) = 1
_type_order(::Type{Transition}) = 2
_type_order(::Type{Pauli}) = 3
_type_order(::Type{Spin}) = 4
_type_order(::Type{Position}) = 5
_type_order(::Type{Momentum}) = 6

# Generic fallbacks for cross-type operator pairs (always distinct sites).
_can_commute(::QSym, ::QSym)::Bool = true
_commute_pair(::QSym, ::QSym)::Tuple{QSym, QSym, CNum, Vector{QSym}} = error("unreachable: cross-type _commute_pair")
_reduce_pair(a::QSym, ::QSym)::Tuple{ReduceKind, QSym, CNum} = (NoReduction, a, _CNUM_ZERO)
_ground_state_expand(::QSym)::Tuple{Bool, Int, Int, Int} = (false, 0, 0, 0)

function _site_compare(a::QSym, b::QSym, ne::Vector{NonEqualPair})::SiteCmp
    ta = typeof(a); tb = typeof(b)
    ta === tb && error("unreachable: same-type _site_compare must be overridden by $ta")
    return _type_order(ta) < _type_order(tb) ? Less : Greater
end
