"""
    fundamental_operators(h::HilbertSpace; names = Symbol[])

Smallest set of operators from which every operator on `h` is built by products and adjoints.

A bosonic mode contributes its annihilation operator; a finite-level system its transitions
``\\sigma^{ij}`` with ``i < j`` and its projectors except the ground state; a
[`CollectiveNLevelSpace`](@ref) all of its collective transitions; Pauli and spin systems their three components; and a pair of
quadratures position and momentum. On a product space the sets of the factors are joined.
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
function fundamental_operators(h::FockSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = isempty(names) ? :a : names[si]
    return Op[Destroy(name, si)]
end

function fundamental_operators(h::NLevelSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = isempty(names) ? :σ : names[si]
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
    name = isempty(names) ? :S : names[si]
    return Op[CollectiveTransition(name, i, j, si) for i in 1:h.n for j in 1:h.n]
end

function fundamental_operators(h::PauliSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = isempty(names) ? :σ : names[si]
    return Op[Pauli(name, 1, si), Pauli(name, 2, si), Pauli(name, 3, si)]
end

function fundamental_operators(h::SpinSpace, si::Int = 1; names::AbstractVector = Symbol[])
    name = isempty(names) ? :S : names[si]
    return Op[Spin(name, 1, si), Spin(name, 2, si), Spin(name, 3, si)]
end

function fundamental_operators(h::PhaseSpace, si::Int = 1; names::AbstractVector = Symbol[])
    x, p = isempty(names) ? (:x, :p) : phase_space_names(names[si])
    return Op[Position(x, si), Momentum(p, si)]
end

phase_space_names(n::Union{Tuple, AbstractVector}) = (n[1], n[2])
phase_space_names(n) = throw(ArgumentError("a PhaseSpace names entry must be a (position, momentum) pair, got $n"))

function fundamental_operators(h::ProductSpace; names::AbstractVector = Symbol[])
    ops = Op[]
    for (i, space) in enumerate(h.spaces)
        space_ops = fundamental_operators(space, i; names = names)
        append!(ops, space_ops)
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
function unique_up_to_adjoint(ops::AbstractVector{<:QField})
    ops_ = copy(ops)
    unique_up_to_adjoint!(ops_)
    return ops_
end


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
    if isempty(names) && h isa ProductSpace
        space_types = typeof.(h.spaces)
        if length(unique(space_types)) != length(space_types)
            alph = 'a':'z'
            names = Symbol.(alph[1:length(h.spaces)])
        end
    end
    fund_ops = fundamental_operators(h; names = names)
    fund_ops = unique(vcat(fund_ops, adjoint.(fund_ops)))

    all_ops = Union{Op, QAdd}[]
    for i in 1:order
        for c in with_replacement_combinations(fund_ops, i)
            c_ = prod(reverse(c))
            iszero(c_) && continue
            if c_ isa QAdd && length(c_.arguments) == 1
                (t, v), = c_.arguments
                push!(all_ops, length(t.ops) == 1 && isone(v) ? only(t.ops) : c_)
            elseif c_ isa QSym
                push!(all_ops, c_)
            end
        end
    end

    unique_up_to_adjoint!(all_ops)
    return all_ops
end

"""
    qadjoint(x)
    qconj(x)

Hermitian conjugate of an operator expression or a coefficient.

For an operator expression this is `x'`. For a coefficient the conjugate is taken through
the expression tree, so that it reaches the individual parameters instead of wrapping the
whole expression in `conj`. An average becomes the average of the adjoint, as with
[`inner_adjoint`](@ref).

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g::Number;

julia> qadjoint(2im * a), qadjoint(g)
(-2im * a', conj(g))
```

See also [`dagger`](@ref), [`inner_adjoint`](@ref).
"""
function qadjoint(v::SymbolicUtils.BasicSymbolic)
    v === Symbolics.IM && return -Symbolics.IM
    (SymbolicUtils.hasmetadata(v, AverageOperator) || is_indexed_sum(v)) && return inner_adjoint(v)
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = SymbolicUtils.arguments(v)
        f isa AvgFunc && return inner_adjoint(v)
        f === conj && return args[1]
        f === expim && return expim_symbolic(-only(args))
        return TermInterface.maketerm(
            typeof(v), f, map(qadjoint, args), TermInterface.metadata(v),
        )
    else
        return conj(v)
    end
end
qadjoint(x::Number) = adjoint(x)
qadjoint(x::Num) = qadjoint(SymbolicUtils.unwrap(x))
qadjoint(x::QField) = adjoint(x)

const qconj = qadjoint

"""
    dagger(x::QField)

Hermitian conjugate of an operator expression, the same as `x'`.

Coefficients are conjugated along with the operators. The name matches `dagger`
from QuantumOptics.jl, so code can use it for symbolic and numeric operators
alike.

# Examples

```jldoctest
julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> @qnumbers a::Destroy(h, 1) σ::Transition(h, 1, 2, 2);

julia> dagger(2im * a' * σ)
-2im * a * σ₂₁
```

See also [`qadjoint`](@ref).
"""
QuantumInterface.dagger(x::QField) = qadjoint(x)

"""
    inner_adjoint(x)

Complex conjugate of an expression of averages, written again as averages.

The conjugate of ``\\langle X \\rangle`` is written as ``\\langle X^\\dagger \\rangle``
rather than as `conj(⟨X⟩)`, so that equations for averages stay expressed in averages of
operators. Coefficients are conjugated as by [`qadjoint`](@ref).

# Examples
```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> @variables g::Number;

julia> inner_adjoint(g * average(a))
⟨a'⟩*conj(g)
```

See also [`qadjoint`](@ref), [`average`](@ref).
"""
function inner_adjoint(v::SymbolicUtils.BasicSymbolic)
    v === Symbolics.IM && return -Symbolics.IM
    if SymbolicUtils.hasmetadata(v, AverageOperator)
        return adjoint_lifted(SymbolicUtils.getmetadata(v, AverageOperator), SymbolicUtils.arguments(v)[1])
    end
    if is_indexed_sum(v)
        return indexed_sum(inner_adjoint(sum_body(v)), sum_indices(v), sum_ne(v))
    end
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        args = SymbolicUtils.arguments(v)
        if f isa AvgFunc
            arg = args[1]
            return adjoint_average_leaf(SymbolicUtils.isconst(arg) ? arg.val : arg)
        elseif f === conj
            return inner_adjoint(args[1])
        elseif f === expim
            return expim_symbolic(-only(args))
        else
            return TermInterface.maketerm(
                typeof(v), f, map(inner_adjoint, args), TermInterface.metadata(v),
            )
        end
    else
        return conj(v)
    end
end
inner_adjoint(x::Number) = conj(x)

function adjoint_lifted(@nospecialize(op), @nospecialize(iv))
    op isa QAdd && iv isa RawExpression && return lift_average(average(adjoint(op)), iv)
    throw(ArgumentError("malformed time-dependent average with operator `$op` and argument `$iv`"))
end

function adjoint_average_leaf(@nospecialize(inner))
    inner isa Op && return average_leaf(adjoint(inner))
    inner isa QAdd && return average_leaf(adjoint(inner))
    throw(ArgumentError("cannot take the adjoint of the averaged argument `$inner`"))
end
inner_adjoint(x::Num) = inner_adjoint(SymbolicUtils.unwrap(x))
inner_adjoint(x::QField) = adjoint(x)

function Base.adjoint(o::Op)
    k = o.kind
    if k === OP_DESTROY
        return Op(OP_CREATE, o.name_id, o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
    elseif k === OP_CREATE
        return Op(OP_DESTROY, o.name_id, o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        return Op(k, o.name_id, o.space_index, o.index, o.l2, o.l1, o.g, o.nlev)
    else
        return o
    end
end

type_order(o::Op) = Int(o.kind)

"""
    order_key(op::Op) -> Tuple

Total ordering key of an operator that ties two operators exactly when they are `isequal`.

The tuple is `(space_index, kind, name rank, index, l1, l2, g, nlev)`: operators group by
subspace, then by kind, then alphabetically by name, and the packed level and axis fields
keep the remaining operators apart. Unlike `hash`, the key is reproducible across sessions
and Julia versions.
"""
order_key(o::Op) = (o.space_index, Int(o.kind), name_rank(o.name_id), index_key(o.index), Int(o.l1), Int(o.l2), Int(o.g), Int(o.nlev))

Base.isless(a::Op, b::Op) = isless(order_key(a), order_key(b))

"""
    ladder(op::Op)

Returns 1 for the lowering members of a canonical pair (annihilation `OP_DESTROY`,
momentum `OP_MOMENTUM`) and 0 otherwise. Used for canonical ordering within
operator product sequences.
"""
ladder(o::Op) = (o.kind === OP_DESTROY || o.kind === OP_MOMENTUM) ? 1 : 0

site_key(o::Op) =
    (o.space_index, o.index, (is_position(o) || is_momentum(o)) ? Int32(0) : o.name_id)

@inline site_family(k::OpKind) =
    (k === OP_DESTROY || k === OP_CREATE) ? 0x00 :
    (k === OP_POSITION || k === OP_MOMENTUM) ? 0x05 :
    UInt8(k)

function site_compare(a::Op, b::Op, ne::Vector{NonEqualPair})
    ka = a.kind; kb = b.kind
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    fa = site_family(ka)
    fa === site_family(kb) || return UInt8(ka) < UInt8(kb) ? Less : Greater
    if !(fa === 0x05 && ka !== kb)
        a.name_id == b.name_id || return name_rank(a.name_id) < name_rank(b.name_id) ? Less : Greater
    end
    a.index == b.index && return Equal
    fa === 0x05 && return Undetermined
    ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

function can_commute(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return false
    elseif ka === OP_TRANSITION && kb === OP_TRANSITION
        return false
    elseif ka === OP_PAULI && kb === OP_PAULI
        return false
    elseif ka === OP_SPIN && kb === OP_SPIN
        return a.l1 <= b.l1
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        a.name_id == b.name_id || return true
        return a.l1 > b.l1 || (a.l1 == b.l1 && a.l2 >= b.l2)
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return false
    else
        return true
    end
end

function commute_pair(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return (b, a, CNUM_ONE, EMPTY_OPS, CNUM_ZERO, EMPTY_OPS)
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return (b, a, CNUM_NEG_IM, EMPTY_OPS, CNUM_ZERO, EMPTY_OPS)
    elseif ka === OP_SPIN && kb === OP_SPIN
        eps = levi_civita[a.l1][b.l1]
        contracted = Op(OP_SPIN, a.name_id, a.space_index, a.index, 6 - a.l1 - b.l1, 0, 0, 0)
        return (b, a, levi_civita_coeff(eps), Op[contracted], CNUM_ZERO, EMPTY_OPS)
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        c1 = a.l2 == b.l1 ? CNUM_ONE : CNUM_ZERO
        c2 = b.l2 == a.l1 ? CNUM_NEG1 : CNUM_ZERO
        r1 = iszero(c1) ? EMPTY_OPS : Op[Op(ka, a.name_id, a.space_index, NO_INDEX, a.l1, b.l2, 0, 0)]
        r2 = iszero(c2) ? EMPTY_OPS : Op[Op(ka, a.name_id, a.space_index, NO_INDEX, b.l1, a.l2, 0, 0)]
        return (b, a, c1, r1, c2, r2)
    else
        return (b, a, CNUM_ZERO, EMPTY_OPS, CNUM_ZERO, EMPTY_OPS)
    end
end

@inline levi_civita_coeff(eps::Int) = eps == 1 ? CNUM_IM : (eps == -1 ? CNUM_NEG_IM : CNUM_ZERO)

may_reduce(a::Op, b::Op) =
    (a.kind === OP_TRANSITION && b.kind === OP_TRANSITION) ||
    (a.kind === OP_PAULI && b.kind === OP_PAULI)

function reduce_pair(a::Op, b::Op)
    ka = a.kind
    if ka === OP_TRANSITION && b.kind === OP_TRANSITION
        (a.name_id == b.name_id && a.space_index == b.space_index && a.index == b.index) ||
            return (NoReduction, a, CNUM_ZERO)
        if a.l2 == b.l1
            new = Op(OP_TRANSITION, a.name_id, a.space_index, a.index, a.l1, b.l2, a.g, a.nlev)
            return (OpReduction, new, CNUM_ONE)
        else
            return (ScalarReduction, a, CNUM_ZERO)
        end
    elseif ka === OP_PAULI && b.kind === OP_PAULI
        (a.name_id == b.name_id && a.space_index == b.space_index && a.index == b.index) ||
            return (NoReduction, a, CNUM_ZERO)
        if a.l1 == b.l1
            return (ScalarReduction, a, CNUM_ONE)
        else
            eps = levi_civita[a.l1][b.l1]
            new = Op(OP_PAULI, a.name_id, a.space_index, a.index, 6 - a.l1 - b.l1, 0, 0, 0)
            return (OpReduction, new, levi_civita_coeff(eps))
        end
    else
        return (NoReduction, a, CNUM_ZERO)
    end
end

function ground_state_expand(op::Op)
    op.kind === OP_TRANSITION || return (false, 0, 0, 0)
    (op.l1 == op.g && op.l2 == op.g) || return (false, 0, 0, 0)
    return (true, Int(op.g), Int(op.nlev), Int(op.space_index))
end
