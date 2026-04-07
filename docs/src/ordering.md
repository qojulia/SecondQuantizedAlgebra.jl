```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Ordering Conventions

SecondQuantizedAlgebra.jl supports two operator ordering conventions that control when commutation relations are applied: **eager** (at multiplication time) and **lazy** (deferred to an explicit call). This page explains the trade-offs and shows why eager ordering is the default.

## Eager ordering (`NormalOrder`)

With [`NormalOrder`](@ref) (the default), every `*` call immediately applies commutation relations and returns a fully normal-ordered [`QAdd`](@ref):

```@example ordering
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
a * a'   # [a, a†] = 1 applied immediately → a†a + 1
```

This means expressions are always in canonical form — no intermediate un-ordered products exist. Subsequent operations (addition, further multiplication) work on already-simplified terms, so like terms are collected early and the expression stays compact.

## Lazy ordering (`LazyOrder`)

With [`LazyOrder`](@ref), multiplication concatenates operators without applying any commutation rules:

```@example ordering
set_ordering!(LazyOrder())
expr = a * a'   # stored as a·a† — no reordering
```

The expression must be explicitly normal-ordered later:

```@example ordering
set_ordering!(NormalOrder())  # hide
normal_order(expr)   # now applies [a, a†] = 1 → a†a + 1
```

```@example ordering
set_ordering!(NormalOrder())  # restore default
nothing # hide
```

## Why eager ordering is the default

In a lazy scheme, intermediate expressions can grow large because like terms are not collected until normal ordering is applied at the end. Consider computing ``H^3`` for a Jaynes–Cummings Hamiltonian: each multiplication doubles or triples the number of terms, and without eager simplification the intermediate ``H^2`` carries many redundant terms that cancel only after final ordering.

With eager ordering, commutation rules fire at each multiplication step. This collects like terms early, keeping intermediate results compact. The result: **eager ordering is consistently faster** for typical quantum optics workflows.

### Benchmark results

The plot below compares eager vs lazy ordering across a range of systems, from simple Fock powers ``(a a^\dagger)^n`` to multi-mode hopping Hamiltonians and lambda systems. Each benchmark measures the total time to construct and order the expression.

![Eager vs Lazy benchmark](assets/eager_vs_lazy.svg)

Key observations:
- **Eager is faster in all cases**, typically by 1.5–5×.
- The advantage grows with expression complexity: more terms means more opportunity for early cancellation.
- For simple cases (e.g. single Fock mode), the difference is small since there are few terms to collect.

The benchmark can be reproduced with:
```
julia --project=benchmark benchmark/eager_vs_lazy.jl
```

### When to use `LazyOrder`

Despite the performance advantage of eager ordering, [`LazyOrder`](@ref) can be useful when:
- You want to inspect **un-ordered operator products** before applying commutation rules.
- You are building expressions where the ordering convention should be chosen later.
- You need to apply [`simplify`](@ref) (algebraic identities only) without commutation swaps.

Switch between conventions at any time:
```julia
set_ordering!(LazyOrder())      # disable eager ordering
set_ordering!(NormalOrder())    # restore eager ordering (default)
```
