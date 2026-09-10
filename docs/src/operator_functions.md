```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Formal Operator Functions

SecondQuantizedAlgebra represents polynomial operator expressions canonically as [`QAdd`](@ref). Functions such as `sin(A)`, `cos(A)`, and ``e^{iA}`` are generally non-polynomial, so they use a separate cold-path [`QExpr`](@ref) representation instead of changing the polynomial term storage.

## Exact formal expressions

```@example operator-functions
using SecondQuantizedAlgebra

h = FockSpace(:f)
@qnumbers a::Destroy(h)
A = a + a'

C = cos(A)
S = sin(A)
U = expim(A)
nothing # hide
```

These are formal operator functions. Constructing them does **not** expand a power series. Ordinary polynomial arithmetic remains on the existing `QAdd` path, while expressions that contain a formal function retain the cold `QExpr` outer type.

Products remain ordered:

```@example operator-functions
left = a * cos(A)
right = cos(A) * a
@assert left != right
nothing # hide
```

`expim(A)` denotes the formal operator ``e^{iA}``. It deliberately requires `A` to be provably Hermitian, matching the package's existing unit-phase semantics. It should not be confused with [`UnitaryTransform`](@ref):

- `expim(A)` is the operator itself;
- `UnitaryTransform` stores an exact compiled adjoint/change-of-frame action.

## Explicit Taylor lowering

When a finite polynomial is required, use the `taylor` generic re-exported from Symbolics:

```@example operator-functions
taylor(cos(A), 0:4)
```

For operator functions, the second argument must currently be a prefix range `0:n`. The result is a canonical `QAdd` with exact coefficients. For example,

```@example operator-functions
@assert taylor(cos(A), 0:4) == 1 - (1 // 2) * A^2 + (1 // 24) * A^4
nothing # hide
```

The two-argument operator method is distinct from the scalar Symbolics API `taylor(f, x, ns)`: an operator expression is not treated as a scalar differentiation variable.

The first implementation intentionally rejects nested formal functions such as `cos(sin(A))` during Taylor lowering. It also rejects a formal function whose argument carries a bound [`Σ`](@ref) scope, because powers such as ``(\sum_i A_i)^2`` require fresh dummy indices rather than silently reusing `i`.

## Structural transformations

Exact algebraic transformations recurse through formal functions. This includes `adjoint`, [`substitute`](@ref), [`change_index`](@ref), [`normal_order`](@ref), [`simplify`](@ref), [`expand`](@ref), [`expand_completeness`](@ref), [`assume_distinct_index`](@ref), and exact [`conjugate`](@ref) operations.

For an exact unitary transformation `T`, for example,

```@example operator-functions
using Symbolics: @variables
@variables θ::Real
T = Rotation(a, θ)
@assert conjugate(cos(A), T) == cos(conjugate(A, T))
nothing # hide
```

No BCH or perturbative truncation is introduced by this operation.

## Numeric conversion

Direct matrix functional calculus for `QExpr` is intentionally not part of this first API. Convert explicitly through a chosen Taylor order:

```@example operator-functions
using QuantumOpticsBase
b = FockBasis(8)
C4 = taylor(cos(A), 0:4)
C4_num = to_numeric(C4, b)
nothing # hide
```

Calling `to_numeric`, `numeric_average`, or `expect` directly on an unlowered formal operator function raises an `ArgumentError` rather than choosing a hidden series order. Direct finite-dimensional matrix `sin`, `cos`, and `exp(im*A)` are tracked separately in issue #268.

## Current boundaries

Formal functions inside [`average`](@ref), symbolic [`Σ`](@ref) over a formal body, and normal/symmetric (Weyl) ordering conversion of an unlowered `QExpr` are outside the first implementation. Lower to a polynomial first when those operations are required. Perturbative BCH, Schrieffer–Wolff, and van Vleck transformations are a separate effective-Hamiltonian/Lie-transform problem rather than part of the `QExpr` representation.
