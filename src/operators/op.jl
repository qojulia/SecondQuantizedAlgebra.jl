"""
    OpKind

Runtime tag distinguishing the seven operator roles carried by the single
concrete leaf type [`Op`](@ref). The integer values double as the cross-family
sort order (see `_type_order`/`order_key`).
"""
@enum OpKind::UInt8 OP_DESTROY OP_CREATE OP_TRANSITION OP_PAULI OP_SPIN OP_POSITION OP_MOMENTUM

"""
    Op <: QSym

The single concrete leaf operator. A `kind::OpKind` tag selects the physical
role (annihilation, creation, transition, Pauli, spin, position, momentum); the
remaining fields are shared storage interpreted per `kind`:

- `name::Symbol`: display name
- `space_index::Int`: position in a [`ProductSpace`](@ref)
- `index::Index`: symbolic site index, or `NO_INDEX`
- `l1, l2, g, nlev::Int32`: packed level/axis data. Transition uses all four
  (`i`, `j`, ground state, number of levels); Pauli/Spin store the axis in `l1`;
  Fock/PhaseSpace leave them zero.

Construct via the role-named functions [`Destroy`](@ref), [`Create`](@ref),
[`Transition`](@ref), [`Pauli`](@ref), [`Spin`](@ref), [`Position`](@ref),
[`Momentum`](@ref); test the role via [`is_destroy`](@ref) and siblings, or read
it with [`optype`](@ref). Collapsing the former per-type hierarchy into one
concrete struct makes the operator vector concrete-eltype, so the per-operator
hooks dispatch statically. See `docs/src/devdocs.md`.
"""
struct Op <: QSym
    kind::OpKind
    name_id::Int32        # interned name (see intern.jl)
    space_index::Int32
    index::Index
    l1::Int32
    l2::Int32
    g::Int32
    nlev::Int32
end

# All fields are isbits, so `Op` is isbits and `Vector{Op}` stores inline with no
# GC scanning. hash/isequal compare the interned ids as pure integers. The custom
# methods are kept for explicit field ordering and the `Index` `===` short-circuit
# (NO_INDEX is a shared const).
Base.hash(o::Op, h::UInt) =
    hash(o.index, hash(o.nlev, hash(o.g, hash(o.l2, hash(o.l1, hash(o.space_index, hash(o.name_id, hash(o.kind, h))))))))
Base.isequal(a::Op, b::Op) =
    a.kind === b.kind && a.name_id == b.name_id && a.space_index == b.space_index &&
    a.l1 == b.l1 && a.l2 == b.l2 && a.g == b.g && a.nlev == b.nlev && a.index == b.index

"""
    operator_name(op::Op) -> Symbol

The display name of `op` (recovered from the intern table).
"""
operator_name(o::Op)::Symbol = _name_from_id(o.name_id)
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

# Shared sentinel for empty operator vectors on hot paths. Never mutated.
const _EMPTY_OPS = Op[]
