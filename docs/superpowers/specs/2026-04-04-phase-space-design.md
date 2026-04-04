# PhaseSpace / Position / Momentum Operators Design

## Goal

Add position and momentum (quadrature) operators to SecondQuantizedAlgebra.jl, following the established type-stable, concrete operator pattern.

## Architecture

Position and Momentum are Hermitian operators on a `PhaseSpace` Hilbert space. They follow the Destroy/Create pattern (two separate types, no axis field). The canonical commutation relation `[X, P] = i` is applied during `simplify` via the existing ordering dispatch. Numeric conversion maps to QuantumOpticsBase ladder operators on a `FockBasis`.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Hilbert space | Separate `PhaseSpace <: HilbertSpace` | Consistent with package pattern; prevents mixing with Fock |
| Operator types | Two structs: `Position`, `Momentum` | Follows Destroy/Create pattern; they are not axis rotations |
| Commutation | `[X, P] = i` in simplify | Standard convention (ℏ=1); consistent with Fock/Spin rules |
| Canonical order | Position before Momentum | "Coordinate before conjugate" convention |
| Numeric basis | `FockBasis` from QuantumOpticsBase | X = (a+a†)/√2, P = i(a†-a)/√2; no new basis type needed |

## Types

### PhaseSpace

```julia
struct PhaseSpace <: HilbertSpace
    name::Symbol
end
```

Minimal, like `PauliSpace`. No additional fields needed.

### Position

```julia
struct Position <: QSym
    name::Symbol
    space_index::Int
end
```

- Hermitian: `adjoint(::Position) = self`
- Construction from `PhaseSpace` (single space, index=1) or `ProductSpace` (explicit index with validation)

### Momentum

```julia
struct Momentum <: QSym
    name::Symbol
    space_index::Int
end
```

- Hermitian: `adjoint(::Momentum) = self`
- Same construction pattern as Position

## Construction

```julia
# Single space
h = PhaseSpace(:q)
x = Position(h, :x)    # space_index = 1
p = Momentum(h, :p)    # space_index = 1

# Product space
hp = FockSpace(:c) ⊗ PhaseSpace(:q)
x = Position(hp, :x)     # auto-finds PhaseSpace at index 2
x = Position(hp, :x, 2)  # explicit index
```

Validation:
- Single `PhaseSpace`: always valid, `space_index = 1`
- `ProductSpace` without index: find unique `PhaseSpace`, error if zero or multiple
- `ProductSpace` with index: validate `spaces[i] isa PhaseSpace`, else `ArgumentError`

## Algebra Rules

### Ordering-dependent (in `_apply_ordering_rule` for `NormalOrder`)

When `Momentum` before `Position` on the same space (same `space_index` and `name`):

```
P·X → X·P - i
```

Implementation: add a branch in `_apply_ordering_rule(a, b, same_space, i, arg_c, ops, ::NormalOrder)`:

```julia
if a isa Momentum && b isa Position && same_space
    # Swap: P·X → X·P - i
    swapped = QSym[ops[1:(i-1)]..., b, a, ops[(i+2):end]...]
    sort!(swapped; lt=canonical_lt)
    term1 = _simplify_qmul(arg_c, swapped, NormalOrder())
    contracted = QSym[ops[1:(i-1)]..., ops[(i+2):end]...]
    sort!(contracted; lt=canonical_lt)
    term2 = _simplify_qmul(-im * arg_c, contracted, NormalOrder())
    all_terms = QMul[term1.arguments..., term2.arguments...]
    return _collect_qadd(all_terms)
end
```

This follows the exact same pattern as the Fock `[a, a†] = 1` rule.

### No ordering-independent reductions

Unlike Pauli (σ²=I), Position and Momentum have no product simplification rules.

## Numeric Conversion

```julia
function to_numeric(op::Position, b::FockBasis; kwargs...)
    return (QuantumOpticsBase.destroy(b) + QuantumOpticsBase.create(b)) / sqrt(2)
end

function to_numeric(op::Momentum, b::FockBasis; kwargs...)
    return im * (QuantumOpticsBase.create(b) - QuantumOpticsBase.destroy(b)) / sqrt(2)
end
```

CompositeBasis embedding uses the existing `space_index` → `LazyTensor` mechanism (no changes needed).

## Ground State Rewriting

No-op for PhaseSpace (continuous spectrum, no ground state projection):

```julia
_apply_ground_state(expr::QAdd, ::PhaseSpace) = expr
```

## Printing

### REPL

Position and Momentum print as their name, like Destroy/Create:

```julia
Base.show(io::IO, x::Position) = write(io, string(x.name))
Base.show(io::IO, x::Momentum) = write(io, string(x.name))
```

### LaTeX

Hat notation for quadrature operators:

```julia
# Position: \hat{x}
@latexrecipe function f(x::Position)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end

# Momentum: \hat{p}
@latexrecipe function f(x::Momentum)
    return Expr(:latexifymerge, "\\hat{$(x.name)}")
end
```

## Files

| File | Action | What |
|---|---|---|
| `src/phase_space.jl` | Create | PhaseSpace, Position, Momentum structs + constructors |
| `src/simplify.jl` | Modify | Add `_apply_ordering_rule` branch for P·X → X·P - i |
| `src/normal_order.jl` | Modify | Add `_apply_ground_state` no-op for PhaseSpace |
| `src/numeric.jl` | Modify | Add `to_numeric` for Position/Momentum on FockBasis |
| `src/printing.jl` | Modify | Add `show` methods |
| `src/latexify_recipes.jl` | Modify | Add LaTeX recipes |
| `src/SecondQuantizedAlgebra.jl` | Modify | Add include + exports |
| `test/phase_space_test.jl` | Create | Tests for all of the above |

## Testing

- Construction: single space, product space, auto-find, validation errors
- Adjoint: `x' == x`, `p' == p`
- Equality/hashing: `Position(:x, 1) == Position(:x, 1)`, hash consistency
- Lazy multiplication: `x * p` returns `QMul`, no commutation applied
- Simplify: `simplify(p * x)` applies `[X, P] = i` → `x·p - i`
- Simplify: `simplify(x * p)` already ordered → `x·p` (no swap)
- Numeric: `to_numeric` on FockBasis matches `(destroy+create)/√2` and `im*(create-destroy)/√2`
- Printing: REPL shows name
- LaTeX: renders with hat notation
- Mixed spaces: Position on space 1, Destroy on space 2 — no interaction
