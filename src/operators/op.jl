"""
    OpKind

Runtime tag distinguishing the eight operator roles carried by the single
concrete leaf type [`Op`](@ref). The integer values double as the cross-family
sort order (see `type_order`/`order_key`).
"""
@enum OpKind::UInt8 OP_DESTROY OP_CREATE OP_TRANSITION OP_PAULI OP_SPIN OP_POSITION OP_MOMENTUM OP_COLLECTIVE_TRANSITION

"""
    Op <: QSym

The single concrete leaf operator. A `kind::OpKind` tag selects the physical
role (annihilation, creation, transition, collective transition, Pauli, spin,
position, momentum); the
remaining fields are shared storage interpreted per `kind`:

- `name::Symbol`: display name
- `space_index::Int`: position in a [`ProductSpace`](@ref)
- `index::Index`: symbolic site index, or `NO_INDEX`
- `l1, l2, g, nlev::Int32`: packed level/axis data. Transition uses all four
  (`i`, `j`, ground state, number of levels); Pauli/Spin store the axis in `l1`;
  Fock/PhaseSpace leave them zero.

Construct via the role-named functions [`Destroy`](@ref), [`Create`](@ref),
[`Transition`](@ref), [`CollectiveTransition`](@ref), [`Pauli`](@ref),
[`Spin`](@ref), [`Position`](@ref), [`Momentum`](@ref); test the role via
[`is_destroy`](@ref) and siblings, or read
it with [`optype`](@ref). Collapsing the former per-type hierarchy into one
concrete struct makes the operator vector concrete-eltype, so the per-operator
hooks dispatch statically. See `docs/src/devdocs.md`.
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

The display name of `op` (recovered from the intern table).
"""
operator_name(o::Op) = name_from_id(o.name_id)
Base.:(==)(a::Op, b::Op) = isequal(a, b)

"""
    is_destroy(o::Op) -> Bool

True iff `o` is a [`Destroy`](@ref) (annihilation) operator. One of the role
predicates that replace `isa` on the collapsed [`Op`](@ref) type; see also
[`is_create`](@ref), [`is_transition`](@ref), [`is_pauli`](@ref),
[`is_spin`](@ref), [`is_position`](@ref), [`is_momentum`](@ref), [`optype`](@ref).
"""
is_destroy(o::Op) = o.kind === OP_DESTROY

"""
    is_create(o::Op) -> Bool

True iff `o` is a [`Create`](@ref) (creation) operator. See [`is_destroy`](@ref).
"""
is_create(o::Op) = o.kind === OP_CREATE

"""
    is_transition(o::Op) -> Bool

True iff `o` is a [`Transition`](@ref) operator. See [`is_destroy`](@ref).
"""
is_transition(o::Op) = o.kind === OP_TRANSITION

"""
    is_collective_transition(o::Op) -> Bool

True iff `o` is a [`CollectiveTransition`](@ref) operator. See [`is_destroy`](@ref).
"""
is_collective_transition(o::Op) = o.kind === OP_COLLECTIVE_TRANSITION

"""
    is_pauli(o::Op) -> Bool

True iff `o` is a [`Pauli`](@ref) operator. See [`is_destroy`](@ref).
"""
is_pauli(o::Op) = o.kind === OP_PAULI

"""
    is_spin(o::Op) -> Bool

True iff `o` is a [`Spin`](@ref) operator. See [`is_destroy`](@ref).
"""
is_spin(o::Op) = o.kind === OP_SPIN

"""
    is_position(o::Op) -> Bool

True iff `o` is a [`Position`](@ref) operator. See [`is_destroy`](@ref).
"""
is_position(o::Op) = o.kind === OP_POSITION

"""
    is_momentum(o::Op) -> Bool

True iff `o` is a [`Momentum`](@ref) operator. See [`is_destroy`](@ref).
"""
is_momentum(o::Op) = o.kind === OP_MOMENTUM

"""
    optype(o::Op) -> OpKind

Return the [`OpKind`](@ref) tag of `o`.
"""
optype(o::Op) = o.kind

"""
    operator_index(o::Op) -> Index

The site [`Index`](@ref) attached to `o`, or the sentinel `NO_INDEX` when `o`
carries no index (test with [`has_index`](@ref)). Public accessor for the
`index` field; pairs with [`operator_name`](@ref) and [`acts_on`](@ref).
"""
operator_index(o::Op) = o.index

"""
    n_levels(o::Op) -> Int
    ground_state(o::Op) -> Int

The number of levels and the ground level of the [`NLevelSpace`](@ref) a [`Transition`](@ref)
acts on. Throw an `ArgumentError` for any other operator.
"""
function n_levels(o::Op)
    is_transition(o) || throw(ArgumentError("`n_levels` is defined for a Transition, got $(o.kind)"))
    return Int(o.nlev)
end
function ground_state(o::Op)
    is_transition(o) ||
        throw(ArgumentError("`ground_state` is defined for a Transition, got $(o.kind)"))
    return Int(o.g)
end
@doc (@doc n_levels) ground_state

"""
    transition_levels(o::Op) -> Tuple{Int, Int}

The levels `(i, j)` of a [`Transition`](@ref) ``|i\\rangle\\langle j|``. Throws an
`ArgumentError` for any other operator.
"""
function transition_levels(o::Op)
    is_transition(o) ||
        throw(ArgumentError("`transition_levels` is defined for a Transition, got $(o.kind)"))
    return (Int(o.l1), Int(o.l2))
end

"""
    operator_axis(o::Op) -> Int

The axis `1`, `2` or `3` of a [`Pauli`](@ref) or [`Spin`](@ref) operator. Throws an
`ArgumentError` for any other operator.
"""
function operator_axis(o::Op)
    (is_pauli(o) || is_spin(o)) ||
        throw(ArgumentError("`operator_axis` is defined for a Pauli or Spin operator, got $(o.kind)"))
    return Int(o.l1)
end

"""
    set_acts_on(o::Op, aon::Integer) -> Op

Rebind `o` to subspace `aon`, returning a copy with its `space_index` set to
`aon` and every other field (role, name, index, levels) preserved. Use it to
embed a single-site operator into a different factor of a [`ProductSpace`](@ref).
The site index, if any, is left untouched; rename it separately with
[`rename`](@ref) when the new subspace needs a matching index.

See also [`acts_on`](@ref), [`operator_index`](@ref).
"""
set_acts_on(o::Op, aon::Integer) =
    Op(o.kind, o.name_id, Int32(aon), o.index, o.l1, o.l2, o.g, o.nlev)

"""
    rename(o::Op, name::Symbol) -> Op

Return a copy of `o` with its display name replaced by `name`, preserving role,
subspace, index, and levels. The `Op` method of [`rename`](@ref rename(::Index, ::Symbol)).
"""
rename(o::Op, name::Symbol) =
    Op(o.kind, intern_name(name), o.space_index, o.index, o.l1, o.l2, o.g, o.nlev)
rename(o::Op, name::AbstractString) = name_must_be_symbol(name)

const EMPTY_OPS = Op[]
