# ClusterSpace Design

## Summary

Add `copy_index::Int` to all `QSym` subtypes and introduce `ClusterSpace{T}` for representing N identical quantum subsystems with correlations tracked up to order M. Operators are expanded into M copies via an explicit `cluster_expand` function. This provides the algebra-level foundation for QuantumCumulants.jl's scaling and cumulant expansion.

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Copy index representation | `copy_index::Int` field on every `QSym` | Uniform, concrete, no new index types or Union dispatch |
| ClusterSpace parameterization | `ClusterSpace{T}` with `N::T` | N can be `Int` or symbolic (`Num`) for QC equation generation |
| Operator expansion | Explicit `cluster_expand(op, ...)` function | Constructor always returns one operator — predictable, type-stable |
| ClusterSpace location | In SecondQuantizedAlgebra | Operator construction, ordering, and commutation are algebra concerns |

## 1. `copy_index` on QSym subtypes

Every `QSym` subtype gains a `copy_index::Int` field as the last field, defaulting to `1`.

### Affected structs

```julia
struct Destroy <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end

struct Create <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end

struct Transition <: QSym
    name::Symbol
    i::Int
    j::Int
    space_index::Int
    copy_index::Int
end

struct Pauli <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
end

struct Spin <: QSym
    name::Symbol
    axis::Int
    space_index::Int
    copy_index::Int
end

struct Position <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end

struct Momentum <: QSym
    name::Symbol
    space_index::Int
    copy_index::Int
end
```

### Constructor changes

All existing constructors from HilbertSpaces pass `copy_index=1`. The inner constructors accept all fields; outer convenience constructors append `1`:

```julia
Destroy(name::Symbol, space_index::Int) = Destroy(name, space_index, 1)
```

Equality, hashing, and adjoint all include `copy_index`.

## 2. Canonical ordering

Sort key in `qmul.jl` changes from `op -> op.space_index` to:

```julia
sort!(args_nc; by = op -> (op.space_index, op.copy_index))
```

Tuple comparison is lexicographic — operators sort first by space, then by copy within the same space.

## 3. Commutation rules

**Key principle:** Operators with same `space_index` but different `copy_index` commute (they act on independent copies).

### `simplify.jl`

Adjacent-pair checks for Fock, Spin, PhaseSpace, Transition, and Pauli already compare `space_index` and `name`. Add `copy_index` to these checks:

```julia
# Example: Fock swap rule
if a isa Destroy && b isa Create &&
    a.space_index == b.space_index &&
    a.copy_index == b.copy_index &&
    a.name == b.name
```

### `commutator.jl`

`_has_space` and `_has_shared_space` compare `(space_index, copy_index)` pairs:

```julia
function _has_space(ops::Vector{QSym}, si::Int, ci::Int)
    for op in ops
        op.space_index == si && op.copy_index == ci && return true
    end
    return false
end

function _has_shared_space(a::Vector{QSym}, b::Vector{QSym})
    for x in a, y in b
        x.space_index == y.space_index &&
            x.copy_index == y.copy_index && return true
    end
    return false
end
```

`commutator(a::QSym, b::QSym)` short-circuit also checks both:

```julia
function commutator(a::QSym, b::QSym)
    (a.space_index == b.space_index && a.copy_index == b.copy_index) || return _ZERO_QADD
    isequal(a, b) && return _ZERO_QADD
    return simplify(a * b - b * a)
end
```

## 4. `ClusterSpace{T}`

```julia
struct ClusterSpace{T} <: HilbertSpace
    original_space::HilbertSpace
    N::T
    order::Int
end
```

- `original_space` — the Hilbert space being replicated (FockSpace, NLevelSpace, etc.)
- `N` — total number of identical copies (Int or symbolic). Used by QC for scaling factors.
- `order` — number of distinct copies to track in the algebra (determines `cluster_expand` output size)

Note: `original_space::HilbertSpace` is abstract-typed. `ClusterSpace` is a construction-time wrapper — never stored on operators or used in hot paths. Skip in `CheckConcreteStructs` test.

### Equality and hashing

```julia
Base.:(==)(a::ClusterSpace, b::ClusterSpace) =
    a.original_space == b.original_space && a.N == b.N && a.order == b.order
Base.hash(a::ClusterSpace, h::UInt) =
    hash(:ClusterSpace, hash(a.original_space, hash(a.N, hash(a.order, h))))
```

## 5. ProductSpace integration

`ClusterSpace` participates in `⊗`:

```julia
h = FockSpace(:c) ⊗ ClusterSpace(NLevelSpace(:atom, 2, 1), N, 2)
```

This creates a `ProductSpace` where `h.spaces[2]` is a `ClusterSpace`.

When constructing operators on a ProductSpace, if `h.spaces[idx]` is a `ClusterSpace`, validate against `h.spaces[idx].original_space`:

```julia
function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    space = h.spaces[idx]
    actual = space isa ClusterSpace ? space.original_space : space
    actual isa FockSpace || throw(ArgumentError(...))
    return Destroy(name, idx, 1)
end
```

Same pattern for `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`.

## 6. `cluster_expand`

```julia
cluster_expand(op::QSym, order::Int) -> Vector{QSym}
```

Returns `order` copies of `op` with `copy_index` set to `1, 2, ..., order` and name suffixed `_1`, `_2`, etc.

```julia
function cluster_expand(op::Destroy, order::Int)
    return [Destroy(Symbol(op.name, :_, i), op.space_index, i) for i in 1:order]
end
```

Same pattern for each QSym subtype. Uses a `_with_copy` internal function to avoid repetition:

```julia
_with_copy(op::Destroy, name::Symbol, ci::Int) = Destroy(name, op.space_index, ci)
_with_copy(op::Create, name::Symbol, ci::Int) = Create(name, op.space_index, ci)
_with_copy(op::Transition, name::Symbol, ci::Int) = Transition(name, op.i, op.j, op.space_index, ci)
_with_copy(op::Pauli, name::Symbol, ci::Int) = Pauli(name, op.axis, op.space_index, ci)
_with_copy(op::Spin, name::Symbol, ci::Int) = Spin(name, op.axis, op.space_index, ci)
_with_copy(op::Position, name::Symbol, ci::Int) = Position(name, op.space_index, ci)
_with_copy(op::Momentum, name::Symbol, ci::Int) = Momentum(name, op.space_index, ci)

function cluster_expand(op::QSym, order::Int)
    return [_with_copy(op, Symbol(op.name, :_, i), i) for i in 1:order]
end
```

ProductSpace convenience:

```julia
function cluster_expand(op::QSym, h::ProductSpace)
    space = h.spaces[op.space_index]
    space isa ClusterSpace || throw(ArgumentError("Space at index $(op.space_index) is not a ClusterSpace"))
    return cluster_expand(op, space.order)
end
```

## 7. `has_cluster`

```julia
has_cluster(::HilbertSpace) = false
has_cluster(::ClusterSpace) = true
function has_cluster(h::ProductSpace)
    for space in h.spaces
        space isa ClusterSpace && return true
    end
    return false
end
```

Exported for QC to detect cluster systems and trigger scaling.

## 8. Printing

When `copy_index > 1`, display it as a subscript suffix. E.g., `a_2` for `Destroy(:a, 1, 2)`. When `copy_index == 1`, display is unchanged (no suffix).

## 9. QC integration surface

What QuantumCumulants.jl can build on:

| SQA provides | QC uses it for |
|---|---|
| `op.copy_index` | Identifying which cluster copy an operator belongs to |
| `ClusterSpace.N` | Scaling factors in equations (`N * g * ⟨...⟩`) |
| `ClusterSpace.order` | Number of independent copies to track |
| `cluster_expand(op, h)` | Generating the M operator copies |
| `has_cluster(h)` | Detecting when to apply scaling |
| `ClusterSpace.original_space` | Validating operator construction |

## 10. Exports

Add to export list: `ClusterSpace`, `cluster_expand`, `has_cluster`

## 11. File layout

| File | Change |
|---|---|
| `src/fock.jl` | Add `copy_index` field, update constructors/equality/hash/adjoint |
| `src/nlevel.jl` | Same |
| `src/pauli.jl` | Same |
| `src/spin.jl` | Same |
| `src/phase_space.jl` | Same |
| `src/qmul.jl` | Update sort key to `(space_index, copy_index)` |
| `src/simplify.jl` | Add `copy_index` to adjacent-pair checks |
| `src/commutator.jl` | Update `_has_space`, `_has_shared_space`, QSym short-circuit |
| `src/printing.jl` | Show copy_index suffix when > 1 |
| `src/latexify_recipes.jl` | Same |
| `src/hilbertspace.jl` | `⊗` support for ClusterSpace |
| New: `src/cluster.jl` | `ClusterSpace`, `cluster_expand`, `has_cluster`, `_with_copy` |
| `src/SecondQuantizedAlgebra.jl` | Include `cluster.jl`, update exports |
| New: `test/cluster_test.jl` | Tests for all of the above |
| `test/concrete_test.jl` | Skip `ClusterSpace` in concrete check |
