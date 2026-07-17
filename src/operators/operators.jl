"""
    fundamental_operators(h::HilbertSpace; names=nothing) -> Vector{Op}

Return the minimal generating set of operators for Hilbert space `h`.

Returns one [`Destroy`](@ref) per [`FockSpace`](@ref); `n(n+1)/2 - 1`
[`Transition`](@ref) operators per [`NLevelSpace`](@ref) (upper-triangular plus
diagonals, excluding the ground-state projector); all `n²`
[`CollectiveTransition`](@ref) operators per [`CollectiveNLevelSpace`](@ref);
three [`Pauli`](@ref) or
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
    return Op[Destroy(name, si)]
end

function fundamental_operators(h::NLevelSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :σ : names[si]
    ops = Op[]
    for i in 1:(h.n)
        for j in i:(h.n)
            (i == j) && i == h.ground_state && continue
            push!(ops, Transition(name, i, j, si, NO_INDEX, h.ground_state, h.n))
        end
    end
    return ops
end

function fundamental_operators(h::CollectiveNLevelSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :S : names[si]
    return Op[CollectiveTransition(name, i, j, si) for i in 1:h.n for j in 1:h.n]
end

function fundamental_operators(h::PauliSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :σ : names[si]
    return Op[Pauli(name, 1, si), Pauli(name, 2, si), Pauli(name, 3, si)]
end

function fundamental_operators(h::SpinSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name = names === nothing ? :S : names[si]
    return Op[Spin(name, 1, si), Spin(name, 2, si), Spin(name, 3, si)]
end

function fundamental_operators(h::PhaseSpace, si::Int = 1; names::Union{Nothing, AbstractVector} = nothing)
    name_pair = names === nothing ? (:x, :p) : names[si]
    return Op[Position(name_pair[1], si), Momentum(name_pair[2], si)]
end

function fundamental_operators(h::ProductSpace; names::Union{Nothing, AbstractVector} = nothing)
    ops = Op[]
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


# Hermitian conjugation

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
    # TODO `Symbolics.IM` (used by `average(::QAdd)` to avoid `Complex{Num}`)
    # has no `conj` fold yet upstream; do it here.
    # https://github.com/JuliaSymbolics/SymbolicUtils.jl/pull/924
    v === Symbolics.IM && return -Symbolics.IM
    if SymbolicUtils.iscall(v)
        f = SymbolicUtils.operation(v)
        # adjoint of a scalar `conj(x)` is `x`; fold rather than nesting into
        # `conj(conj(x))`, which does not simplify and survives downstream (e.g.
        # leaving dead time-dependent terms after a coupling is substituted to 0).
        # Mirrors the `f === conj` case already handled in `inner_adjoint`.
        f === conj && return SymbolicUtils.arguments(v)[1]
        args = map(qadjoint, SymbolicUtils.arguments(v))
        return TermInterface.maketerm(typeof(v), f, args, TermInterface.metadata(v))
    else
        return conj(v)
    end
end
qadjoint(x::Number) = adjoint(x)
qadjoint(x::Num) = qadjoint(SymbolicUtils.unwrap(x))
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
    # TODO See `qadjoint`: upstream `conj(IM)` doesn't fold yet. https://github.com/JuliaSymbolics/SymbolicUtils.jl/pull/924
    v === Symbolics.IM && return -Symbolics.IM
    if SymbolicUtils.hasmetadata(v, AverageOperator)
        # Lifted average ⟨X⟩(iv): the inside-adjoint is ⟨X†⟩(iv), still a function of the
        # same iv. Re-lift the adjoint so time-dependence is preserved.
        op = SymbolicUtils.getmetadata(v, AverageOperator)
        iv = SymbolicUtils.arguments(v)[1]
        return _lift_average(_average(adjoint(op)), iv)
    end
    if is_indexed_sum(v)
        return _indexed_sum(inner_adjoint(_sum_body(v)), _sum_indices(v), _sum_ne(v))
    end
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


# Adjoint, ordering, ladder

function Base.adjoint(o::Op)
    k = o.kind
    if k === OP_DESTROY
        return Op(OP_CREATE, o.name_id, o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
    elseif k === OP_CREATE
        return Op(OP_DESTROY, o.name_id, o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        # |i⟩⟨j|† = |j⟩⟨i|: swap the packed levels.
        return Op(k, o.name_id, o.space_index, o.index, o.l2, o.l1, o.g, o.nlev)
    else
        return o   # Pauli, Spin, Position, Momentum are Hermitian
    end
end

# Instance form used by `_full_op_key`; the enum value is the cross-type order.
_type_order(o::Op) = Int(o.kind)

"""
    order_key(op::Op) -> Tuple

Total, identity-faithful ordering key: a comparable tuple that orders operators
reproducibly and ties two exactly when they are `isequal`. The leading
`space_index` groups operators by subspace; `kind` orders across families within
a subspace; the packed `l1..nlev` fields keep distinct levels/axes distinct.
"""
order_key(o::Op) = (o.space_index, Int(o.kind), _name_rank(o.name_id), _index_key(o.index), Int(o.l1), Int(o.l2), Int(o.g), Int(o.nlev))

Base.isless(a::Op, b::Op) = isless(order_key(a), order_key(b))

"""
    ladder(op::Op)

Returns 1 for the lowering members of a canonical pair (annihilation `OP_DESTROY`,
momentum `OP_MOMENTUM`) and 0 otherwise. Used for canonical ordering within
operator product sequences.
"""
ladder(o::Op) = (o.kind === OP_DESTROY || o.kind === OP_MOMENTUM) ? 1 : 0


# Operator hooks (single concrete methods branching on `kind`)

# Site family: Fock {Destroy,Create} and PhaseSpace {Position,Momentum} compare
# cross-role within the family; the others are singleton families. Distinct
# families fall back to the `kind` integer order (= the old `_type_order`).
@inline _site_family(k::OpKind) =
    (k === OP_DESTROY || k === OP_CREATE) ? 0x00 :
    (k === OP_POSITION || k === OP_MOMENTUM) ? 0x05 :
    UInt8(k)

function _site_compare(a::Op, b::Op, ne::Vector{NonEqualPair})
    ka = a.kind; kb = b.kind
    a.space_index == b.space_index || return a.space_index < b.space_index ? Less : Greater
    fa = _site_family(ka)
    fa === _site_family(kb) || return UInt8(ka) < UInt8(kb) ? Less : Greater
    # PhaseSpace x-vs-p ignores name (conjugate variables on one site); every
    # other same-family pair distinguishes by name.
    if !(fa === 0x05 && ka !== kb)
        a.name_id == b.name_id || return _name_rank(a.name_id) < _name_rank(b.name_id) ? Less : Greater
    end
    a.index == b.index && return Equal
    # PhaseSpace never consults `ne`; Fock/Transition/Pauli/Spin do.
    fa === 0x05 && return Undetermined
    _ne_contains(ne, a.index, b.index) && return a.index < b.index ? Less : Greater
    return Undetermined
end

function _can_commute(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return false                 # a·a† carries the identity residual
    elseif ka === OP_TRANSITION && kb === OP_TRANSITION
        return false                 # always compose
    elseif ka === OP_PAULI && kb === OP_PAULI
        return false                 # always compose
    elseif ka === OP_SPIN && kb === OP_SPIN
        return a.l1 <= b.l1          # commute only in ascending axis order
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        a.name_id == b.name_id || return true
        return a.l1 > b.l1 || (a.l1 == b.l1 && a.l2 >= b.l2)
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return false                 # P·X = X·P - i
    else
        return true
    end
end

# `_commute_pair` returns the swapped pair followed by two residual
# `(coefficient, operators)` slots. Existing roles use only the first slot.
function _commute_pair(a::Op, b::Op)
    ka = a.kind; kb = b.kind
    if ka === OP_DESTROY && kb === OP_CREATE
        return (b, a, _CNUM_ONE, _EMPTY_OPS, _CNUM_ZERO, _EMPTY_OPS) # aa† = a†a + 1
    elseif ka === OP_MOMENTUM && kb === OP_POSITION
        return (b, a, _to_cnum(-im), _EMPTY_OPS, _CNUM_ZERO, _EMPTY_OPS) # P·X = X·P - i·I
    elseif ka === OP_SPIN && kb === OP_SPIN
        # [Sj, Sk] = iϵⱼₖₗSl; the residual is the contracted spin on the third axis.
        eps = _levi_civita[a.l1][b.l1]
        contracted = Op(OP_SPIN, a.name_id, a.space_index, a.index, 6 - a.l1 - b.l1, 0, 0, 0)
        return (b, a, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE), Op[contracted], _CNUM_ZERO, _EMPTY_OPS)
    elseif ka === OP_COLLECTIVE_TRANSITION && kb === OP_COLLECTIVE_TRANSITION
        # [Sⁱʲ,Sᵏˡ] = δⱼₖSⁱˡ - δₗᵢSᵏʲ.
        c1 = a.l2 == b.l1 ? _CNUM_ONE : _CNUM_ZERO
        c2 = b.l2 == a.l1 ? _CNUM_NEG1 : _CNUM_ZERO
        r1 = _iszero_cnum(c1) ? _EMPTY_OPS : Op[Op(ka, a.name_id, a.space_index, NO_INDEX, a.l1, b.l2, 0, 0)]
        r2 = _iszero_cnum(c2) ? _EMPTY_OPS : Op[Op(ka, a.name_id, a.space_index, NO_INDEX, b.l1, a.l2, 0, 0)]
        return (b, a, c1, r1, c2, r2)
    else
        return (b, a, _CNUM_ZERO, _EMPTY_OPS, _CNUM_ZERO, _EMPTY_OPS)
    end
end

# Reduce-pass gate: only Transition·Transition and Pauli·Pauli compose locally,
# so non-reducing pairs skip the field checks below.
_may_reduce(a::Op, b::Op) =
    (a.kind === OP_TRANSITION && b.kind === OP_TRANSITION) ||
    (a.kind === OP_PAULI && b.kind === OP_PAULI)

function _reduce_pair(a::Op, b::Op)
    ka = a.kind
    if ka === OP_TRANSITION && b.kind === OP_TRANSITION
        # σⁱʲ · σᵏˡ = δⱼₖ σⁱˡ.
        (a.name_id == b.name_id && a.space_index == b.space_index && a.index == b.index) ||
            return (NoReduction, a, _CNUM_ZERO)
        if a.l2 == b.l1
            new = Op(OP_TRANSITION, a.name_id, a.space_index, a.index, a.l1, b.l2, a.g, a.nlev)
            return (OpReduction, new, _CNUM_ONE)
        else
            return (ScalarReduction, a, _CNUM_ZERO)    # δ_{j,k} = 0
        end
    elseif ka === OP_PAULI && b.kind === OP_PAULI
        # σⱼ·σₖ = δⱼₖI + iϵⱼₖₗσₗ.
        (a.name_id == b.name_id && a.space_index == b.space_index && a.index == b.index) ||
            return (NoReduction, a, _CNUM_ZERO)
        if a.l1 == b.l1
            return (ScalarReduction, a, _CNUM_ONE)     # σⱼ² = 1
        else
            eps = _levi_civita[a.l1][b.l1]
            new = Op(OP_PAULI, a.name_id, a.space_index, a.index, 6 - a.l1 - b.l1, 0, 0, 0)
            return (OpReduction, new, _mul_cnum(_to_cnum(im * eps), _CNUM_ONE))
        end
    else
        return (NoReduction, a, _CNUM_ZERO)
    end
end

# Ground-state expansion: σᵍᵍ -> 1 - Σ_{k≠g} σᵏᵏ (Transition only).
function _ground_state_expand(op::Op)
    op.kind === OP_TRANSITION || return (false, 0, 0, 0)
    (op.l1 == op.g && op.l2 == op.g) || return (false, 0, 0, 0)
    return (true, Int(op.g), Int(op.nlev), Int(op.space_index))
end
