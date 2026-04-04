# NLevel, Pauli, and Spin Operators — Design Spec

**Date:** 2026-04-04
**Scope:** Add NLevelSpace/Transition, PauliSpace/Pauli, SpinSpace/Spin to the type-stable core
**Depends on:** v2 redesign (Fock-only core on `redesign-v2` branch)

## Motivation

The Fock-only core provides the type-stable foundation. This extension adds the
three remaining operator families needed for quantum optics: N-level atoms,
Pauli two-level systems, and collective spin operators. All follow the same
architectural pattern: new `HilbertSpace` subtype, new `QSym` subtypes, algebra
rules applied lazily in `normal_order()`.

## New Hilbert Spaces

```julia
struct NLevelSpace <: HilbertSpace
    name::Symbol
    n::Int              # number of levels
    ground_state::Int   # ground state index for population conservation
end

struct PauliSpace <: HilbertSpace
    name::Symbol
end

struct SpinSpace <: HilbertSpace
    name::Symbol
    spin::Rational{Int}  # spin quantum number: 1//2, 1, 3//2, ...
end
```

Equality and hashing follow the same pattern as `FockSpace`. All fields are
concrete. `ProductSpace{Tuple}` already accepts any `HilbertSpace` subtype.

## New Operators

```julia
struct Transition <: QSym
    name::Symbol
    i::Int              # bra level  (represents |i⟩⟨j|)
    j::Int              # ket level
    space_index::Int
end

struct Pauli <: QSym
    name::Symbol
    axis::Int           # 1=x, 2=y, 3=z
    space_index::Int
end

struct Spin <: QSym
    name::Symbol
    axis::Int           # 1=x, 2=y, 3=z
    space_index::Int
end
```

All structs are fully concrete with no type parameters.

### Construction from Hilbert Spaces

```julia
# NLevel
Transition(h::NLevelSpace, name::Symbol, i::Int, j::Int) = Transition(name, i, j, 1)
Transition(h::ProductSpace, name::Symbol, i::Int, j::Int, idx::Int)  # validates spaces[idx] isa NLevelSpace

# Pauli
Pauli(h::PauliSpace, name::Symbol, axis::Int) = Pauli(name, axis, 1)
Pauli(h::ProductSpace, name::Symbol, axis::Int, idx::Int)  # validates spaces[idx] isa PauliSpace

# Spin
Spin(h::SpinSpace, name::Symbol, axis::Int) = Spin(name, axis, 1)
Spin(h::ProductSpace, name::Symbol, axis::Int, idx::Int)  # validates spaces[idx] isa SpinSpace
```

The `@qnumbers` macro works unchanged:
```julia
h = NLevelSpace(:atom, 3, 1)
@qnumbers σ::Transition(h, 1, 2)
```

### Adjoints

- `Transition(:σ, i, j)' → Transition(:σ, j, i)` — swap bra/ket
- `Pauli(:σ, ax)' → Pauli(:σ, ax)` — Hermitian
- `Spin(:S, ax)' → Spin(:S, ax)` — Hermitian

### Canonical Ordering

Same as Fock: `canonical_lt` sorts by `space_index` only. Operators on the
same space preserve insertion order. The `ladder` function returns a default
for non-Fock operators (they don't participate in Fock-style normal ordering).

```julia
ladder(::Transition) = 0
ladder(::Pauli) = 0
ladder(::Spin) = 0
```

## Algebra Rules — All Lazy

`*` never applies any algebra rules for any operator type. It collects
operators into `QMul`, sorted by `space_index`.

`normal_order()` applies all simplification rules for adjacent operators on
the same space (same `space_index` and same `name`):

### Fock (existing)

`Destroy(a) · Create(a)` → swap + identity: `Create(a) · Destroy(a) + 1`

### Transition

`Transition(i,j) · Transition(k,l)` on same space:
- `j == k` → `Transition(i, l)` (composition: `|i⟩⟨j| · |j⟩⟨k| = |i⟩⟨k|`)
- `j ≠ k` → `0` (orthogonality: `|i⟩⟨j| · |k⟩⟨l| = 0` when `j ≠ k`)

Ground state rewriting: when `normal_order()` encounters `Transition(g, g)`
where `g` is the ground state of the associated `NLevelSpace`, it rewrites:
`|g⟩⟨g| = 1 - Σ_{k≠g} |k⟩⟨k|`

This requires `normal_order` to know the `NLevelSpace` associated with a
`Transition`. Since operators don't store the Hilbert space, the user must
pass it: `normal_order(expr, spaces)` where `spaces` maps `space_index` to
the `HilbertSpace`. For single-space expressions, `normal_order(expr, h)`.
The existing no-argument `normal_order(expr)` skips ground state rewriting
(Fock-only behavior is unchanged).

### Pauli

`Pauli(j) · Pauli(k)` on same space:
- `j == k` → `QMul(1, QSym[])` (identity: `σⱼ² = I`)
- `j ≠ k` → `QMul(im * ϵⱼₖₗ, QSym[Pauli(l)])` (product rule: `σⱼσₖ = iϵⱼₖₗσₗ`)

Where `ϵⱼₖₗ` is the Levi-Civita symbol computed via `Combinatorics.levicivita`
or a simple inline implementation.

### Spin

`Spin(j) · Spin(k)` on same space (with `j > k`, i.e. out of canonical axis order):
Apply commutation relation `[Sⱼ, Sₖ] = iϵⱼₖₗSₗ`:
- Swap to `Spin(k) · Spin(j)` and add commutator term `iϵⱼₖₗ · Spin(l)`

This is analogous to Fock normal ordering but for axis indices instead of
ladder indices.

## Numerical Conversion

Extend `to_numeric`:

| Symbolic | QuantumOpticsBase |
|---|---|
| `Transition(:σ, i, j)` | `transition(basis, i, j)` |
| `Pauli(:σ, 1)` | `sigmax(basis)` |
| `Pauli(:σ, 2)` | `sigmay(basis)` |
| `Pauli(:σ, 3)` | `sigmaz(basis)` |
| `Spin(:S, 1)` | `0.5 * sigmax(basis)` |
| `Spin(:S, 2)` | `0.5 * sigmay(basis)` |
| `Spin(:S, 3)` | `0.5 * sigmaz(basis)` |

Basis validation:

| HilbertSpace | QuantumOpticsBase Basis |
|---|---|
| `NLevelSpace(:a, n, g)` | `NLevelBasis(n)` |
| `PauliSpace(:p)` | `SpinBasis(1//2)` |
| `SpinSpace(:s, S)` | `SpinBasis(S)` |

Composite basis embedding via `LazyTensor` uses `space_index` to select
the sub-basis, same as Fock.

## Printing

### REPL

```
Transition(:σ, 1, 2, 1)  → "σ₁₂"
Pauli(:σ, 1, 1)          → "σx"
Pauli(:σ, 2, 1)          → "σy"
Pauli(:σ, 3, 1)          → "σz"
Spin(:S, 1, 1)           → "Sx"
Spin(:S, 2, 1)           → "Sy"
Spin(:S, 3, 1)           → "Sz"
NLevelSpace(:a, 3, 1)    → "ℋ(a)"
PauliSpace(:p)            → "ℋ(p)"
SpinSpace(:s, 1//2)       → "ℋ(s)"
```

### LaTeX

```
Transition(:σ, 1, 2)  → σ^{12}    (or σ_{12} with toggle)
Pauli(:σ, 1)           → σ_x
Spin(:S, 3)            → S_z
```

The `transition_superscript(::Bool)` function toggles between superscript
and subscript for Transition level indices. Default is superscript (`^`).
Stored as a module-level `Ref{Symbol}`.

## Changes to Existing Code

### `canonical_lt` (src/fock.jl)

No change needed — already sorts by `space_index` only. New operator types
automatically work.

### `normal_order` (src/normal_order.jl)

Extended to handle new operator pairs. The scanning loop in
`_normal_order_qmul` checks for:
1. Destroy/Create pairs (existing)
2. Transition/Transition pairs → compose or zero
3. Pauli/Pauli pairs → identity or `iϵ·σ`
4. Spin/Spin pairs out of axis order → swap + commutator

A new `normal_order(expr, spaces)` overload enables ground state rewriting
for Transitions.

### `@qnumbers` (src/macros.jl)

No change — the macro already forwards all constructor arguments. `@qnumbers σ::Transition(h, 1, 2)` expands to `Transition(h, :σ, 1, 2)`.

### `to_numeric` (src/numeric.jl)

New methods added for each operator type and basis combination.

### Printing/LaTeX (src/printing.jl, src/latexify_recipes.jl)

New `show` methods and `@latexrecipe` definitions for each type.

## Design Invariants (unchanged from v2)

1. **Type stability:** `*` returns `QMul{T}`, `+` returns `QAdd{T}`.
2. **No `Vector{Any}`:** all containers use concrete element types.
3. **No abstract field types:** all struct fields are concrete.
4. **Canonical ordering:** `args_nc` sorted by `space_index`, order within space preserved.
5. **Lazy multiplication:** `*` never applies algebra rules.
6. **Closed algebra:** `normal_order()` and `simplify()` return `QAdd`.

## Dependencies

`Combinatorics` needs to be re-added to `[deps]` for `levicivita` (Levi-Civita
symbol used in Pauli and Spin algebra). Alternatively, implement a simple
3-element Levi-Civita inline to avoid the dependency.
