```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Symbolic Sums and Indices

Many physical systems contain multiple identical elements with different parameters. Symbolic summation lets you write compact Hamiltonians with indexed operators and derive equations once, then instantiate for specific values later.

A canonical example is the Tavis-Cummings Hamiltonian describing ``N`` two-level atoms coupled to a cavity mode:

```math
H_\mathrm{TC} = \omega_c a^\dagger a + \sum_i^N \omega_i \sigma_i^{22} + \sum_i^N g_i \left(a^\dagger \sigma_i^{12} + a\, \sigma_i^{21}\right).
```

!!! note "Replaces ClusterSpace"
    Earlier versions of SecondQuantizedAlgebra.jl shipped a `ClusterSpace` type — inherited from QuantumCumulants.jl, from which this package was refactored — that wrapped a single-site Hilbert space with a fixed number of identical copies and a correlation-truncation order. QuantumCumulants.jl has since deprecated that workflow in favor of indexed operators (compare its old [cluster-based superradiant laser tutorial](https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant-laser/) with the [indexed superradiant laser tutorial](https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/)), and we have followed suit: `ClusterSpace`, `cluster_expand`, the integer-instantiation helper `insert_index`, and the cluster-aware `to_numeric(op, b, ranges)` overload have all been removed. The recommended way to handle many identical subsystems is to keep a plain Hilbert space and use a symbolic [`Index`](@ref) of range ``N`` together with [`Σ`](@ref). The number of copies stays symbolic through equation derivation, and the cumulant truncation order is supplied to the downstream solver (`meanfield(...; order=k)`) rather than being baked into the Hilbert space.


## Index

An [`Index`](@ref) represents a symbolic summation variable. It is constructed from a Hilbert space, a name, a range (symbolic or numeric), and the subspace it acts on:

```@example sums
using SecondQuantizedAlgebra

@variables N

ha = NLevelSpace(:atoms, 2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

i = Index(h, :i, N, ha)       # index over the NLevelSpace
j = Index(h, :j, N, 2)        # equivalent: specify subspace by position
nothing # hide
```


## Indexed operators

Associate an operator with an [`Index`](@ref) using [`IndexedOperator`](@ref):

```@example sums
σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y, 2), z)
σ(2, 1, i)
```

Similarly, [`IndexedVariable`](@ref) creates symbolic parameters with an index:

```@example sums
gi = IndexedVariable(:g, i)
```


## Concrete slots

After an index has been unrolled (e.g. by `QuantumCumulants.evaluate`), you
often want to refer to a specific position. Calling an `Index` value with an
`Integer` returns a fresh per-slot `Index` named `Symbol(i.name, "_", k)` with
the same range and subspace. Pass it through [`IndexedOperator`](@ref) to build
the concrete-site operator that matches the post-unroll naming convention:

```@example sums
σ(2, 1, i(3))
```


## Summations

Use [`Σ`](@ref) (or equivalently [`∑`](@ref)) to create a symbolic sum over an index:

```@example sums
Σ(σ(2, 2, i), i)
```

A third argument specifies indices that are not equal to the summation index:

```@example sums
Σ(σ(2, 2, i), i, [j])
```


### Diagonal splitting

When two indices acting on the same Hilbert space meet inside a product, the diagonal term (where both indices are equal) is automatically separated:

```@example sums
k = Index(h, :k, N, ha)
l = Index(h, :l, N, ha)

Σ(σ(2, 1, k) * σ(1, 2, l), k, l)
```

This also happens when a sum is multiplied by an indexed operator on the same space:

```@example sums
Σ(σ(2, 2, k), k) * σ(2, 1, l)
```

The diagonal term may produce a ground-state projector — for instance, the diagonal of ``\sigma^{12}_k \cdot \sigma^{21}_l`` at ``k=l`` is ``\sigma^{11}_l``. Same-site composition that produces a ground-state projector triggers eager completeness expansion via ``\sigma^{gg} = 1 - \sum_{k \neq g}\sigma^{kk}``, so the canonical basis stays ``\{\sigma^{ij} : (i,j) \neq (g,g)\} \cup \{1\}``. Standalone ``\sigma^{gg}`` built directly with `Transition(h, :σ, g, g)` stays atomic until it enters a `*`; you can also trigger the rewrite by hand via [`expand_completeness`](@ref).


## Free indices and `assume_distinct_index`

Two indexed operators on the same Hilbert subspace whose indices are *both free* (neither bound by a [`Σ`](@ref)) have an undetermined site relationship — the algebra cannot tell whether the user means "distinct atomic sites" or "two index variables that may coincide". The conservative default is the second: the pair stays in its physical order and no same-site collapse fires.

When you do want to assert that two free indices denote distinct sites, use [`assume_distinct_index`](@ref). It augments every term's per-term inequality constraints with the supplied pairs, re-canonicalizes so the resolved pairs sort deterministically, and runs [`expand_completeness`](@ref) on any ground-state projectors that emerge from same-site composition under the new constraint:

```@example sums
ki = Index(h, :ki, N, ha)
li = Index(h, :li, N, ha)

assume_distinct_index(σ(1, 2, ki) * σ(2, 1, li), [(ki, li)])
```

When at least one index is bound by a `Σ`, the diagonal-splitting machinery shown above handles the boundary case automatically and no explicit assumption is needed.


## Example: Tavis-Cummings model

Putting it all together for ``N`` two-level atoms in a cavity:

```@example sums
@variables Δ κ

@qnumbers a::Destroy(h, 1)

H = Δ * a' * a + Σ(gi * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
```
