```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Collective Atoms — Two Computational Approaches

A system of `N` identical atoms in collective dynamics — Dicke models, superradiance, Tavis–Cummings — can be tackled in two complementary ways. Both describe the same physics; they differ in the *computational strategy* and the kind of result you get out at the end.

| | Indexed `Σ` (cumulant approach) | [`CollectiveTransition`](@ref) (exact symmetric subspace) |
|---|---|---|
| Numerical strategy | Cumulant / mean-field truncation | Exact diagonalization on `ManyBodyBasis` |
| Hilbert space at numeric layer | Single representative atom + symbolic site index | Full symmetric subspace, dim ``\binom{N+n-1}{n-1}`` |
| Site-dependent parameters (e.g. ``g_i``) | Yes, via [`IndexedVariable`](@ref) | No (atoms identical by construction) |
| `N` scaling | Large `N` (truncation order chosen at solve time) | Modest `N` (``\le`` thousands for two-level atoms) |
| Output | Mean-field / cumulant equations of motion | State vector or density matrix on the symmetric subspace |
| Algebra at the symbolic level | Per-atom transition rules + diagonal split + completeness | Closed ``\mathfrak{su}(N)`` Lie algebra |

Pick the one that matches your computational strategy, and don't mix them on the same `NLevelSpace`.

## Cumulant approach — indexed `Σ`

Use this when you want **mean-field-style equations of motion** and `N` is large enough that exact simulation is infeasible. Atoms can carry site-dependent symbolic parameters via [`IndexedVariable`](@ref). The cumulant truncation order is supplied to the downstream solver — the algebra layer keeps `N` symbolic.

```julia
using SecondQuantizedAlgebra
@variables N

ha = NLevelSpace(:atom, 2)          # single representative atom
hc = FockSpace(:cavity)
h  = hc ⊗ ha

i  = Index(h, :i, N, ha)
σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
gi = IndexedVariable(:g, i)

@qnumbers a::Destroy(h, 1)
H_TC = Σ(gi * (a' * σ(1, 2) + a * σ(2, 1)), i)
```

Heisenberg equations follow from `commutator(H, op)`; the diagonal-split and completeness machinery handles the per-atom algebra correctly. See [Symbolic Sums and Indices](@ref) for the full set of rules.

## Exact symmetric subspace — `CollectiveTransition`

Use this when you want **exact dynamics in the symmetric subspace** and `N` is small enough that the polynomial-size basis is tractable. Atoms are identical by construction (no site labels); the operator algebra is the closed ``\mathfrak{su}(N)`` Lie algebra; numeric simulation goes through `QuantumOpticsBase.ManyBodyBasis`.

```math
\dim \mathcal{H}_\text{sym} = \binom{N + n - 1}{n - 1}
```

For `N = 100` atoms with `n = 3` levels: full product space ``\sim 10^{47}`` (intractable), symmetric subspace ``5151`` (very tractable).

```julia
h = NLevelSpace(:atom, 3)
S(i, j) = CollectiveTransition(h, :S, i, j)

S(1, 2) * S(2, 3)         # eagerly normal-ordered: S(1, 3) + S(2, 3) * S(1, 2)
```

The collective operators satisfy ``[S^{ij}, S^{kl}] = \delta_{jk} S^{il} - \delta_{li} S^{kj}``, and the swap rule fires eagerly under [`NormalOrder`](@ref) so derivations stay in the closed collective form rather than expanding into per-atom sums.

### Numeric conversion

Build a `ManyBodyBasis` over the single-atom basis with `bosonstates`:

```julia
using QuantumOpticsBase
n_atoms = 100
b1 = NLevelBasis(3)
b  = ManyBodyBasis(b1, bosonstates(b1, n_atoms))    # symmetric subspace

to_numeric(S(1, 2), b)   # many-body operator on the symmetric subspace
```

The number of atoms `N` enters only here, not at the algebra layer.

### What the algebra does *not* do

- **Single-atom completeness is not applied** to collective operators. Per-atom completeness ``\sigma^{gg} = 1 - \sum_{k \neq g}\sigma^{kk}`` does not hold for ``S^{gg}``: collectively, ``\sum_j S^{jj} = N \cdot I``. Since `N` lives at the numeric layer, the rewrite is deferred to it. This is why `simplify(expr, h)` does *not* expand ``S^{gg}`` even though it expands ``\sigma^{gg}``.
- **No mixing with [`Transition`](@ref) on the same `NLevelSpace`.** `Transition` is a per-atom operator (the "single representative atom" of the cumulant approach); `CollectiveTransition` is a sum-over-atoms operator on the symmetric subspace. These live in different conceptual spaces, and the package will not throw if you write `Transition(h, :σ, 1, 2) * CollectiveTransition(h, :S, 2, 1)`, but the result is meaningless.
