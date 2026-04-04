# Commutator Feature Design

## Summary

Implement a type-stable `commutator(a, b)` function that computes `[a, b] = a*b - b*a`,
applies `simplify` to resolve commutation relations, and always returns `QAdd{CT}`.

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Return type | Always `QAdd{CT}` | Type stability; matches `simplify` contract |
| Simplification | Always simplify internally | Physics expectations (`[a, a†] = 1`); keeps nested commutators compact |
| Zero representation | `QAdd{Complex{Int}}` | `Complex{Int}` is the natural bottom type; promotes cleanly with everything; matches `simplify(::QSym)` |
| Short-circuits | Keep space-index checks | Cheap `Int` comparisons avoid allocations; compounds in nested use |
| Bilinearity promotion | Promote across sub-results via `reduce(+, ...)` | Uses existing `QAdd + QAdd` promotion for free |

## Public API

```julia
commutator(a, b) -> QAdd
```

Single function. Exported from `SecondQuantizedAlgebra`.

## Canonical Zero

```julia
const _ZERO_QADD = QAdd(QMul{Complex{Int}}[QMul(zero(Complex{Int}), QSym[])])
```

Pre-allocated constant returned by all short-circuit paths. No per-call allocation for trivially commuting operands.

## Dispatch Table

### Scalar methods (return `_ZERO_QADD`)

```julia
commutator(::Number, ::Number)
commutator(::Number, ::QField)
commutator(::QField, ::Number)
```

Scalars commute with everything. Three methods to avoid `::Any` ambiguity.

### QSym, QSym

```julia
function commutator(a::QSym, b::QSym) -> QAdd
```

- `a.space_index != b.space_index` → `_ZERO_QADD`
- `isequal(a, b)` → `_ZERO_QADD`
- else → `simplify(a * b - b * a)`

### QMul, QSym (and symmetric)

```julia
function commutator(a::QMul, b::QSym) -> QAdd
function commutator(a::QSym, b::QMul) -> QAdd
```

- Check if any operator in `args_nc` shares `b.space_index` (resp. `a.space_index`)
- No match → `_ZERO_QADD`
- else → `simplify(a * b - b * a)`

### QMul, QMul

```julia
function commutator(a::QMul, b::QMul) -> QAdd
```

- Compute space indices of both sides, check intersection
- Empty intersection → `_ZERO_QADD`
- else → `simplify(a * b - b * a)`

### QAdd bilinearity (3 methods)

```julia
function commutator(a::QAdd, b::QField) -> QAdd
function commutator(a::QField, b::QAdd) -> QAdd
function commutator(a::QAdd, b::QAdd) -> QAdd
```

Distribute over `QAdd` arguments:
- Call `commutator` on each (term, operand) pair
- Combine results with `reduce(+, sub_results)`
- Uses existing `QAdd + QAdd` which handles type promotion automatically
- If all sub-commutators are zero, the sum of `_ZERO_QADD`s is a valid `QAdd{Complex{Int}}`

## File Layout

- **New file:** `src/commutator.jl`
- **Included in** `SecondQuantizedAlgebra.jl` after `simplify.jl` (depends on `simplify`)
- **Exported:** `commutator` added to the `export` list

## Type Flow Example

```
commutator(a::Destroy, a'::Create)
  → a.space_index == a'.space_index, !isequal(a, a')
  → simplify(a * a' - a' * a)
  → simplify returns QAdd{Complex{Int}}
  → result: QAdd{Complex{Int}} with term QMul(1+0im, QSym[])  (i.e., the identity = 1)

commutator(a::Destroy, b::Destroy)  where a.space_index != b.space_index
  → short-circuit: _ZERO_QADD  (QAdd{Complex{Int}})

commutator(a + a', b::Destroy)
  → commutator(a_term, b) + commutator(a'_term, b)
  → QAdd{Complex{Int}} + QAdd{Complex{Int}}
  → QAdd{Complex{Int}}
```

## Dependencies

- `simplify` from `src/simplify.jl` — called for non-trivial commutators
- `QAdd`, `QMul`, `QSym` from `src/qadd.jl`, `src/qmul.jl`, `src/types.jl`
- Existing `QAdd + QAdd` promotion from `src/qadd.jl`
