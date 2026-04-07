# Symbolic Sums and Indices

Many physical systems contain multiple identical elements with different parameters. Symbolic summation lets you write compact Hamiltonians with indexed operators and derive equations once, then instantiate for specific values later.

A canonical example is the Tavis-Cummings Hamiltonian describing ``N`` two-level atoms coupled to a cavity mode:

```math
H_\mathrm{TC} = \omega_c a^\dagger a + \sum_i^N \omega_i \sigma_i^{22} + \sum_i^N g_i \left(a^\dagger \sigma_i^{12} + a\, \sigma_i^{21}\right).
```


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


## Example: Tavis-Cummings model

Putting it all together for ``N`` two-level atoms in a cavity:

```@example sums
@variables Δ κ

@qnumbers a::Destroy(h, 1)

H = Δ * a' * a + Σ(gi * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
```
