"""
    OpKind

Kind of an operator, as returned by [`optype`](@ref).

The values are `OP_DESTROY`, `OP_CREATE`, `OP_TRANSITION`, `OP_COLLECTIVE_TRANSITION`,
`OP_PAULI`, `OP_SPIN`, `OP_POSITION` and `OP_MOMENTUM`. They are public but not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: optype

julia> optype(Destroy(FockSpace(:cavity), :a)') === SecondQuantizedAlgebra.OP_CREATE
true
```

See also [`optype`](@ref), [`is_destroy`](@ref).
"""
@enum OpKind::UInt8 OP_DESTROY OP_CREATE OP_TRANSITION OP_PAULI OP_SPIN OP_POSITION OP_MOMENTUM OP_COLLECTIVE_TRANSITION

"""
    Op <: QSym

Type of every single operator, whatever its kind.

Operators are made with [`Destroy`](@ref), [`Create`](@ref), [`Transition`](@ref),
[`CollectiveTransition`](@ref), [`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref)
and [`Momentum`](@ref), never with `Op` directly. Since all of them share this
type, use [`optype`](@ref) or the predicates such as [`is_destroy`](@ref) to tell
them apart instead of `isa`.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_destroy, is_create

julia> a = Destroy(FockSpace(:cavity), :a);

julia> a isa Op, is_destroy(a), is_create(a')
(true, true, true)
```

See also [`optype`](@ref), [`operator_name`](@ref), [`acts_on`](@ref).
"""
struct Op <: QSym
    kind::OpKind
    name_id::Int32
    space_index::Int32
    index::Index
    l1::Int32
    l2::Int32
    g::Int32
    nlev::Int32
end

const NOIDX_H = hash(NO_INDEX, zero(UInt))
function Base.hash(o::Op, h::UInt)
    h = hash(o.kind, h)
    h = hash(o.name_id, h)
    h = hash(o.space_index, h)
    h = hash(o.l1, h)
    h = hash(o.l2, h)
    h = hash(o.g, h)
    h = hash(o.nlev, h)
    return o.index === NO_INDEX ? hash(NOIDX_H, h) : hash(o.index, h)
end
Base.isequal(a::Op, b::Op) =
    a.kind === b.kind && a.name_id == b.name_id && a.space_index == b.space_index &&
    a.l1 == b.l1 && a.l2 == b.l2 && a.g == b.g && a.nlev == b.nlev && a.index == b.index

"""
    operator_name(op::Op) -> Symbol

Name of the operator `op`, as given at construction.

# Examples

```jldoctest
julia> a = Destroy(FockSpace(:cavity), :a);

julia> operator_name(a), operator_name(a')
(:a, :a)
```

See also [`operator_index`](@ref), [`acts_on`](@ref).
"""
operator_name(o::Op) = name_from_id(o.name_id)
Base.:(==)(a::Op, b::Op) = isequal(a, b)

"""
    is_destroy(o::Op) -> Bool

Whether `o` is an annihilation operator, see [`Destroy`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_destroy

julia> a = Destroy(FockSpace(:cavity), :a);

julia> is_destroy(a), is_destroy(a')
(true, false)
```

See also [`is_create`](@ref), [`is_transition`](@ref),
[`is_collective_transition`](@ref), [`is_pauli`](@ref), [`is_spin`](@ref),
[`is_position`](@ref), [`is_momentum`](@ref), [`optype`](@ref).
"""
is_destroy(o::Op) = o.kind === OP_DESTROY

"""
    is_create(o::Op) -> Bool

Whether `o` is a creation operator, see [`Create`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_create

julia> a = Destroy(FockSpace(:cavity), :a);

julia> is_create(a'), is_create(a)
(true, false)
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_create(o::Op) = o.kind === OP_CREATE

"""
    is_transition(o::Op) -> Bool

Whether `o` is a transition operator, see [`Transition`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_transition

julia> σ = Transition(NLevelSpace(:atom, 2), :σ, 1, 2);

julia> is_transition(σ)
true
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_transition(o::Op) = o.kind === OP_TRANSITION

"""
    is_collective_transition(o::Op) -> Bool

Whether `o` is a collective transition operator, see [`CollectiveTransition`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_collective_transition, is_transition

julia> S = CollectiveTransition(CollectiveNLevelSpace(:atoms, 2), :S, 1, 2);

julia> is_collective_transition(S), is_transition(S)
(true, false)
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_collective_transition(o::Op) = o.kind === OP_COLLECTIVE_TRANSITION

"""
    is_pauli(o::Op) -> Bool

Whether `o` is a Pauli operator, see [`Pauli`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_pauli

julia> σx = Pauli(PauliSpace(:qubit), :σ, 1);

julia> is_pauli(σx)
true
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_pauli(o::Op) = o.kind === OP_PAULI

"""
    is_spin(o::Op) -> Bool

Whether `o` is a spin operator, see [`Spin`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_spin

julia> Sx = Spin(SpinSpace(:S), :S, 1);

julia> is_spin(Sx)
true
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_spin(o::Op) = o.kind === OP_SPIN

"""
    is_position(o::Op) -> Bool

Whether `o` is a position quadrature, see [`Position`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_position, is_momentum

julia> x = Position(PhaseSpace(:osc), :x);

julia> is_position(x), is_momentum(x)
(true, false)
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_position(o::Op) = o.kind === OP_POSITION

"""
    is_momentum(o::Op) -> Bool

Whether `o` is a momentum quadrature, see [`Momentum`](@ref).

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: is_momentum

julia> p = Momentum(PhaseSpace(:osc), :p);

julia> is_momentum(p)
true
```

See also [`is_destroy`](@ref), [`optype`](@ref).
"""
is_momentum(o::Op) = o.kind === OP_MOMENTUM

"""
    optype(o::Op) -> OpKind

Kind of the operator `o`.

The result is an [`OpKind`](@ref) constant, public but not exported, so write it
qualified, e.g. `SecondQuantizedAlgebra.OP_DESTROY`. For a single check the predicates
such as [`is_destroy`](@ref) read better.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: optype

julia> σ = Transition(NLevelSpace(:atom, 2), :σ, 2, 1);

julia> optype(σ) == SecondQuantizedAlgebra.OP_TRANSITION
true
```

See also [`Op`](@ref), [`is_destroy`](@ref).
"""
optype(o::Op) = o.kind

"""
    operator_index(o::Op) -> Index

Site [`Index`](@ref) of the operator `o`. An operator without an index returns a
placeholder for which [`has_index`](@ref) is `false`.

# Examples

```jldoctest
julia> @variables N;

julia> h = NLevelSpace(:atom, 2); i = Index(h, :i, N, h);

julia> σ = Transition(h, :σ, 2, 1);

julia> operator_index(IndexedOperator(σ, i)) == i, has_index(operator_index(σ))
(true, false)
```

See also [`IndexedOperator`](@ref), [`operator_name`](@ref), [`acts_on`](@ref).
"""
operator_index(o::Op) = o.index

@noinline wrong_kind(f::Symbol, expected::String, o::Op) =
    throw(ArgumentError("`$f` is defined for $expected, got $(o.kind)"))

"""
    n_levels(o::Op) -> Int
    ground_state(o::Op) -> Int

The number of levels and the ground level of the [`NLevelSpace`](@ref) a [`Transition`](@ref)
acts on. Throws an `ArgumentError` for any other operator.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: n_levels, ground_state

julia> σ = Transition(NLevelSpace(:atom, 3), :σ, 3, 1);

julia> n_levels(σ), ground_state(σ)
(3, 1)
```

See also [`transition_levels`](@ref).
"""
function n_levels(o::Op)
    is_transition(o) || wrong_kind(:n_levels, "a Transition", o)
    return Int(o.nlev)
end
function ground_state(o::Op)
    is_transition(o) || wrong_kind(:ground_state, "a Transition", o)
    return Int(o.g)
end
@doc (@doc n_levels) ground_state

"""
    transition_levels(o::Op) -> Tuple{Int, Int}

The levels `(i, j)` of a [`Transition`](@ref) ``|i\\rangle\\langle j|``. Throws an
`ArgumentError` for any other operator.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: transition_levels

julia> transition_levels(Transition(NLevelSpace(:atom, 3), :σ, 3, 1))
(3, 1)
```

See also [`n_levels`](@ref).
"""
function transition_levels(o::Op)
    is_transition(o) || wrong_kind(:transition_levels, "a Transition", o)
    return (Int(o.l1), Int(o.l2))
end

"""
    operator_axis(o::Op) -> Int

The axis `1`, `2` or `3` of a [`Pauli`](@ref) or [`Spin`](@ref) operator. Throws an
`ArgumentError` for any other operator.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: operator_axis

julia> operator_axis(Pauli(PauliSpace(:q), :σ, 2))
2
```
"""
function operator_axis(o::Op)
    (is_pauli(o) || is_spin(o)) || wrong_kind(:operator_axis, "a Pauli or Spin operator", o)
    return Int(o.l1)
end

"""
    set_acts_on(op, subsystem::Int)

Copy of the operator `op` acting on another factor of a product space.

Name, kind, levels and index are kept. This reuses an operator built on a single space
inside a product space, or moves it to another factor. The function is not exported.

# Examples
```jldoctest
julia> using SecondQuantizedAlgebra: set_acts_on

julia> σ = Transition(NLevelSpace(:atom, 2), :σ, 2, 1);

julia> h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2);

julia> Destroy(h, :a, 1) * set_acts_on(σ, 2)
a * σ₂₁
```

See also [`acts_on`](@ref), [`rename`](@ref).
"""
set_acts_on(o::Op, aon::Integer) =
    Op(o.kind, o.name_id, Int32(aon), o.index, o.l1, o.l2, o.g, o.nlev)

rename(o::Op, name::Symbol) =
    Op(o.kind, name_id(name), o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
rename(o::Op, name::AbstractString) = name_must_be_symbol(name)

with_kind(o::Op, kind::OpKind) = Op(kind, o.name_id, o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)

"""
    with_levels(o::Op, l1, l2) -> Op

Copy of `o` with the level fields replaced: the levels of a transition, or the axis `l1`
(with `l2 = 0`) of a Pauli or spin operator.
"""
with_levels(o::Op, l1::Integer, l2::Integer) =
    Op(o.kind, o.name_id, o.space_index, o.index, Int32(l1), Int32(l2), o.g, o.nlev)

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
function IndexedOperator(op::Op, i::Index)
    is_collective_transition(op) && has_index(i) &&
        throw(ArgumentError("CollectiveTransition is already collective and cannot be indexed"))
    return Op(op.kind, op.name_id, op.space_index, i, op.l1, op.l2, op.g, op.nlev)
end

const EMPTY_OPS = Op[]
