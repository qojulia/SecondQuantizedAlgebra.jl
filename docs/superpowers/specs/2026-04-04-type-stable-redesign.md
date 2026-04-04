# SecondQuantizedAlgebra.jl v2 — Type-Stable Redesign

**Date:** 2026-04-04
**Scope:** Fock space only (Destroy/Create) — minimal core with full type stability
**Inspiration:** KeldyshContraction.jl's algebra layer

## Motivation

The current package is pervasively type-unstable:
- `QMul.arg_c` is untyped (`Any`)
- `QMul.args_nc` is `Vector{Any}`
- `QAdd.arguments` is `Vector{Any}`
- `*` eagerly applies commutation relations, returning `Union{Int, QSym, QMul, QAdd}` — the core source of type instability
- Operators carry parametric types for Hilbert space, name, acts-on, and metadata (`Destroy{H,S,A,M}`)

This redesign replaces all of that with a concrete, static, type-stable architecture.

## Type Hierarchy

```
QField (abstract)
├── QSym (abstract) — fundamental operators (leaves)
│   ├── Destroy (concrete struct)
│   └── Create (concrete struct)
└── QTerm (abstract) — compound expressions
    ├── QMul{T<:Number} — product with typed prefactor
    └── QAdd{T<:Number} — sum of QMul{T}
```

## Core Structs

### Operators

```julia
struct Destroy <: QSym
    name::Symbol
    space_index::Int
end

struct Create <: QSym
    name::Symbol
    space_index::Int
end
```

No type parameters. No metadata. No Hilbert space reference.
`space_index` defaults to `1` for single-space systems.
`adjoint(::Destroy)` returns `Create` (same name, same index) and vice versa.

### Compound Expressions

```julia
struct QMul{T<:Number} <: QTerm
    arg_c::T
    args_nc::Vector{QSym}
end

struct QAdd{T<:Number} <: QTerm
    arguments::Vector{QMul{T}}
end
```

- A bare scalar is `QMul{Int}(5, QSym[])` (empty `args_nc`).
- A bare operator in a sum is wrapped as `QMul{Int}(1, [op])`.
- `QMul(0, [...])` stays as-is — zero-pruning happens in `simplify()`, not at construction.
- `QMul(1, [op])` stays as-is — no unwrapping to bare `QSym`.

### Canonical Ordering

Operators in `args_nc` are sorted by `(ladder, space_index, name)` where:
- `ladder(::Create) = 0`
- `ladder(::Destroy) = 1`

Creation operators come first, then sorted by space index, then by name.

## Hilbert Spaces (Construction-Time Only)

```julia
abstract type HilbertSpace end

struct FockSpace <: HilbertSpace
    name::Symbol
end

struct ProductSpace{T<:Tuple{Vararg{HilbertSpace}}} <: HilbertSpace
    spaces::T
end
```

`ProductSpace` uses a `Tuple` for fully concrete storage. Hilbert spaces are only used at operator construction time for validation — they are not stored on operators.

### Tensor Product

```julia
⊗(a::FockSpace, b::FockSpace) = ProductSpace((a, b))
⊗(a::ProductSpace, b::FockSpace) = ProductSpace((a.spaces..., b))
⊗(a::FockSpace, b::ProductSpace) = ProductSpace((a, b.spaces...))
⊗(a::ProductSpace, b::ProductSpace) = ProductSpace((a.spaces..., b.spaces...))
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a ⊗ b, c...)
tensor(args::Vararg{HilbertSpace}) = ⊗(args...)
```

### Construction API

```julia
# Single space
Destroy(h::FockSpace, name::Symbol) = Destroy(name, 1)
Create(h::FockSpace, name::Symbol) = Create(name, 1)

# Product space — user specifies subspace index
function Destroy(h::ProductSpace, name::Symbol, idx::Int)
    @assert 1 <= idx <= length(h.spaces)
    @assert h.spaces[idx] isa FockSpace
    Destroy(name, idx)
end
```

### @qnumbers Macro

```julia
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
# expands to: a = Destroy(h, :a)
```

Same user-facing syntax as current package.

## Arithmetic Operations

All operations have predictable, concrete return types.

### Multiplication (`*` → always `QMul`)

| Operation | Return Type |
|---|---|
| `QSym * QSym` | `QMul{Int}` |
| `QSym * Number` | `QMul` |
| `Number * QSym` | `QMul` |
| `QMul * QSym` | `QMul` |
| `QMul * QMul` | `QMul` |
| `QMul * Number` | `QMul` |
| `QAdd * Number` | `QAdd` |
| `QAdd * QSym` | `QAdd` |
| `QAdd * QMul` | `QAdd` |
| `QAdd * QAdd` | `QAdd` |

**Key rule:** `*` never applies commutation relations. It collects operators, sorts by canonical order, and multiplies prefactors.

### Addition (`+` → always `QAdd`)

| Operation | Return Type |
|---|---|
| `QSym + QSym` | `QAdd{Int}` |
| `QMul + QMul` | `QAdd` |
| `QMul + QSym` | `QAdd` |
| `QAdd + QMul` | `QAdd` |
| `QAdd + QAdd` | `QAdd` |
| `QField + Number` | `QAdd` |

Bare `QSym` operands are wrapped as `QMul{Int}(1, [op])`. Scalar operands become `QMul{T}(val, QSym[])`.

### Other

- `a / b::Number` → `(1/b) * a`
- `a ^ n::Integer` → repeated `*` (returns `QMul`)
- `-a` → `QMul(-1, ...)` or `QMul(-a.arg_c, a.args_nc)`
- `a - b` → `a + (-b)`

## Normal Ordering & Simplification

### `normal_order()`

Applies the bosonic commutation relation `[a, a†] = 1`. Always returns `QAdd`.

```julia
normal_order(::QSym) → QAdd
normal_order(::QMul) → QAdd
normal_order(::QAdd) → QAdd
```

Algorithm: walk through `args_nc`, whenever `Destroy` appears left of `Create` on the same space, swap them and add the identity term. Repeat until fully ordered. Does **not** collect like terms.

### `simplify()`

Groups `QMul` terms with identical `args_nc`, sums their `arg_c` prefactors (works with both plain `Number` and `BasicSymbolic`). Removes zero-prefactor terms. Returns `QAdd`.

```julia
simplify(QAdd([QMul(g, [a†,a]), QMul(h, [a†,a])]))
→ QAdd([QMul(g + h, [a†,a])])
```

### SymbolicUtils / Symbolics Integration

- `SymbolicUtils.simplify(::QField)` — calls `simplify()`, then applies `SymbolicUtils.simplify` to each `arg_c` prefactor.
- `Symbolics.expand(::QField)` — distributes products, then applies `Symbolics.expand` to each `arg_c` prefactor.

## TermInterface Integration

### QSym (leaves)

```julia
SymbolicUtils.iscall(::QSym) = false
TermInterface.head(::QField) = :call
TermInterface.metadata(::QSym) = nothing
```

### QMul (products)

```julia
SymbolicUtils.iscall(::QMul) = true
SymbolicUtils.operation(::QMul) = (*)
SymbolicUtils.arguments(a::QMul) = Any[a.arg_c, a.args_nc...]  # heterogeneous by TermInterface design
TermInterface.metadata(::QMul) = nothing

function TermInterface.maketerm(::Type{<:QMul}, ::typeof(*), args, metadata)
    args_c = filter(x -> !(x isa QField), args)
    args_nc = filter(x -> x isa QField, args)
    arg_c = isempty(args_c) ? 1 : *(args_c...)
    return QMul(arg_c, QSym[args_nc...])
end
```

Note: `SymbolicUtils.arguments` returns `Vector{Any}` because TermInterface inherently mixes c-numbers and operators. This is confined to the TermInterface boundary — it does not leak into the internal representation.

### QAdd (sums)

```julia
SymbolicUtils.iscall(::QAdd) = true
SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
TermInterface.metadata(::QAdd) = nothing
TermInterface.maketerm(::Type{<:QAdd}, ::typeof(+), args, metadata) = QAdd(args)
```

### Type Promotion

```julia
SymbolicUtils.symtype(x::T) where {T<:QField} = T
SymbolicUtils.promote_symtype(f, T::Type{<:QField}, S::Type{<:Number}) = T
```

No metadata stored anywhere — always returns `nothing`.

## C-Numbers (Symbolic Parameters)

**Drop all custom c-number types.** Use Symbolics.jl directly:

```julia
using Symbolics
@variables g h       # complex symbolic parameters
@variables ω::Real   # real symbolic parameter
```

These work as `arg_c` in `QMul{T}` where `T` is `Num` or `BasicSymbolic`.

**Removed:** `Parameter`, `RealParameter`, `CNumber`, `RNumber`, `@cnumbers`, `@rnumbers`, `cnumber()`, `rnumber()`.

## Numerical Conversion (QuantumOpticsBase)

### `to_numeric()`

```julia
to_numeric(op::Destroy, basis::FockBasis) → destroy(basis)
to_numeric(op::Create, basis::FockBasis)  → create(basis)
to_numeric(op::QMul, bases...)  → arg_c * prod(to_numeric.(args_nc, ...))
to_numeric(op::QAdd, bases...)  → sum(to_numeric.(arguments, ...))
```

For product spaces, `space_index` determines which basis to use. Operators are embedded into the composite space via `LazyTensor`.

### `numeric_average()`

Computes `tr(ρ * to_numeric(op))` for a given density matrix.

### Basis Mapping

| Symbolic | QuantumOpticsBase |
|---|---|
| `FockSpace` | `FockBasis(N)` |
| `ProductSpace` | `CompositeBasis(bases...)` |
| `Destroy` | `destroy(basis)` |
| `Create` | `create(basis)` |

User provides the basis with cutoff dimension.

## Printing & LaTeX

### REPL Display

```
Destroy(:a, 1)  → "a"
Create(:a, 1)   → "a†"
QMul(3, [Create(:a,1), Destroy(:a,1)])  → "3 * a† * a"
QAdd([QMul(1, [Create(:a,1), Destroy(:a,1)]), QMul(1, QSym[])])  → "a† * a + 1"
```

- Scalar `QMul` (empty `args_nc`) prints just the prefactor.
- Unit prefactor (`arg_c == 1`) is omitted.
- Single-space systems omit the space index.
- Smart sign handling: `QMul(-1, [a])` prints as `-a`.

### LaTeX via Latexify.jl

`@latexrecipe` for each type:
- `Destroy` → `a`
- `Create` → `a^{\dagger}`
- `QMul` → `\text{prefactor} \cdot op_1 op_2 \ldots`
- `QAdd` → `term_1 + term_2 + \ldots`

### Hilbert Space Display

```
FockSpace(:cavity)  → "ℋ(cavity)"
ProductSpace((FockSpace(:a), FockSpace(:b)))  → "ℋ(a) ⊗ ℋ(b)"
```

## Dependencies

| Dependency | Purpose |
|---|---|
| SymbolicUtils | TermInterface integration, symbolic type promotion |
| Symbolics | `@variables` for c-numbers, `expand` |
| TermInterface | Expression tree interface |
| QuantumOpticsBase | Numerical matrix conversion |
| SciMLBase | Required by QuantumOpticsBase |
| Latexify | LaTeX rendering |
| LaTeXStrings | LaTeX string type |
| MacroTools | `@qnumbers` macro parsing |
| Combinatorics | Utility combinatorics (kept for future use) |

### Removed Dependencies

None — all current deps are retained. Custom c-number types are removed but their deps remain for other uses.

## What Is Deferred (Future Scope)

The following features from the current package are **not** in this initial redesign. They will be added incrementally on top of the stable core:

- NLevelSpace / Transition operators
- PauliSpace / Pauli operators
- SpinSpace / Spin operators
- PhaseSpace / Position / Momentum operators
- ClusterSpace / cluster expansion
- Averages (`average()`, `undo_average()`)
- Symbolic indexing (Index, SingleSum, DoubleSum)
- NumberedOperator / IndexedOperator
- IndexedVariable / DoubleIndexedVariable
- IndexedAverageSum / IndexedAverageDoubleSum
- `commutator()` function
- `find_operators()` / `fundamental_operators()`
- ModelingToolkit weak dependency extension

## Design Invariants

These must hold at all times and should be enforced by tests:

1. **Type stability:** every arithmetic operation returns a concrete type — `*` returns `QMul{T}`, `+` returns `QAdd{T}`.
2. **No `Vector{Any}`:** all containers use concrete element types.
3. **No abstract field types:** all struct fields are concrete.
4. **Canonical ordering:** `args_nc` is always sorted by `(ladder, space_index, name)`.
5. **Lazy multiplication:** `*` never applies commutation relations.
6. **Closed algebra:** `normal_order()` and `simplify()` return `QAdd`.
