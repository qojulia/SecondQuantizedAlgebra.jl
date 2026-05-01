```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Collective Atoms — Two Regimes

A system of `N` atoms can be modeled in two qualitatively different ways depending on whether the atoms are *individually addressable*. SecondQuantizedAlgebra.jl supports both — pick the one that matches your physics, and don't mix them within the same atomic Hilbert space.

## Distinguishable atoms — indexed `Σ`

Use this when atoms have site-resolved properties: position-dependent couplings ``g_i``, individual addressing, disorder, etc. Each atom lives in its own copy of an [`NLevelSpace`](@ref). Symbolic sums over a site index keep the bookkeeping compact at the algebra level.

```julia
using SecondQuantizedAlgebra
@variables N

ha = NLevelSpace(:atom, 2)
hc = FockSpace(:cavity)
h  = hc ⊗ ha

i  = Index(h, :i, N, ha)
σ(α, β) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
gi = IndexedVariable(:g, i)

@qnumbers a::Destroy(h, 1)
H_TC = Σ(gi * (a' * σ(1, 2) + a * σ(2, 1)), i)
```

The full Hilbert space is ``\dim(\mathcal{H}_{\text{atom}})^{N}`` — exponential in `N`. Numeric simulation is tractable only for small `N` (typically `N ≲ 10` for `n=2` levels). The diagonal-split machinery, completeness, and per-atom σᵍᵍ rewriting are documented in [Symbolic Sums and Indices](@ref) and [Ordering Conventions](@ref).

## Indistinguishable atoms — `CollectiveTransition`

Use this when the dynamics preserves permutation symmetry: Dicke physics, single-mode driving with no atom-resolved coupling, superradiance with identical atoms, etc. The state stays in the *symmetric (bosonic) subspace*, whose dimension is polynomial in `N`:

```math
\dim \mathcal{H}_\text{sym} = \binom{N + n - 1}{n - 1}
```

For `N = 100` atoms with `n = 3` levels: full Hilbert space ``\sim 10^{47}``, symmetric subspace ``5151``. The latter is the only thing you can actually simulate.

[`CollectiveTransition`](@ref) ``S^{ij} = \sum_k \sigma_k^{ij}`` is the operator on this subspace:

```julia
h = NLevelSpace(:atom, 3)
S(i, j) = CollectiveTransition(h, :S, i, j)

S(1, 2) * S(2, 3)         # eagerly normal-ordered: S(1, 3) + S(2, 3) * S(1, 2)
```

These satisfy the ``\mathfrak{su}(N)`` Lie algebra ``[S^{ij}, S^{kl}] = \delta_{jk} S^{il} - \delta_{li} S^{kj}``, and the swap rule fires eagerly under [`NormalOrder`](@ref) so derivations stay in the closed collective form.

### Numeric conversion

Numeric simulation requires a `QuantumOpticsBase.ManyBodyBasis` over the single-atom basis with `bosonstates`:

```julia
using QuantumOpticsBase
n_atoms = 100
b1 = NLevelBasis(3)
b  = ManyBodyBasis(b1, bosonstates(b1, n_atoms))    # symmetric subspace

to_numeric(S(1, 2), b)   # many-body operator on the symmetric subspace
```

The number of atoms `N` enters only here, not at the algebra layer.

## Comparison

| | Distinguishable atoms | Indistinguishable atoms |
|---|---|---|
| Representation | [`Σ`](@ref) of [`IndexedOperator`](@ref) on a [`ProductSpace`](@ref) | [`CollectiveTransition`](@ref) on a single [`NLevelSpace`](@ref) |
| Atoms | individually addressable | symmetric subspace, no atom labels |
| Hilbert-space dim | ``n^N`` | ``\binom{N+n-1}{n-1}`` |
| Numeric scaling | exponential in `N` | polynomial in `N` |
| Algebra | per-atom Transition rules + diagonal split | ``\mathfrak{su}(N)`` Lie algebra (closed) |
| Ground-state completeness | applies (``\sigma^{gg} = 1 - \sum_{k \neq g}\sigma^{kk}``) | does not apply (``\sum_j S^{jj} = N \cdot I``, deferred to numeric layer) |
| Typical use case | Tavis–Cummings with disorder, lattice models, individual addressing | Dicke model, superradiant laser with identical atoms |

## Don't mix

[`Transition`](@ref) and [`CollectiveTransition`](@ref) on the *same* [`NLevelSpace`](@ref) are physically incompatible representations — pick one regime per atomic subspace. The package will not throw if you write `Transition(h, :σ, 1, 2) * CollectiveTransition(h, :S, 2, 1)`, but the result is meaningless: `Transition` lives in the full Hilbert space and `CollectiveTransition` in the symmetric subspace, and there is no canonical embedding between them at the algebra level.
