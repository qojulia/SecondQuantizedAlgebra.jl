```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Many-body systems

An ensemble of identical atoms or an array of coupled modes is written once per site and summed. Consider ``N`` two-level atoms in a cavity, the Tavis–Cummings model

```math
H = \Delta a^\dagger a + \sum_{i=1}^{N} g_i \left( a^\dagger \sigma_i^{12} + a \sigma_i^{21} \right).
```

The site label ``i`` and the number of atoms ``N`` both stay symbolic, so equations are derived once for any ``N`` and any distribution of couplings. When the atoms are indistinguishable, collective operators describe the same ensemble without any site label.

```@contents
Pages = ["many_body.md"]
Depth = 2:3
```

## Indices

A site label is an [`Index`](@ref). It belongs to one subsystem of the product space and runs over a range that may be a number or a symbol. For the Tavis–Cummings model the index runs over the atoms, the second factor.

```@example tc
using SecondQuantizedAlgebra
@variables N Δ
h = FockSpace(:cavity) ⊗ NLevelSpace(:atoms, 2)
i = Index(h, :i, N, 2)
j = Index(h, :j, N, 2)
nothing # hide
```

```@docs
Index
```

!!! note "Replaces `ClusterSpace`"
    Earlier versions inherited `ClusterSpace` from QuantumCumulants.jl, which wrapped a fixed number of identical copies of a space. It has been removed, together with `cluster_expand` and `insert_index`, in favor of an [`Index`](@ref) with a symbolic range. The number of copies then stays symbolic, and the truncation order of a cumulant expansion is set in the solver rather than in the Hilbert space.

## Indexed operators and parameters

An operator is placed on a site with [`IndexedOperator`](@ref). Wrapping it in a function of the levels and the index keeps the notation close to ``\sigma_i^{kl}``.

```@example tc
@qnumbers a::Destroy(h, 1)
σ(k, l, m) = IndexedOperator(Transition(h, :σ, k, l), m)
σ(2, 1, i)
```

Operators on different sites commute, while operators on the same site obey the local operator algebra. To retain individual atom labels, use indexed [`Transition`](@ref) operators as above and [`Σ`](@ref) to form their collective sums.

A [`CollectiveTransition`](@ref) instead represents an entire ensemble; see [Collective systems](@ref). An index on such an operator could label separate ensembles, each with its own collective algebra, but `IndexedOperator` does not currently support this. For a fixed number of ensembles, give each its own [`CollectiveNLevelSpace`](@ref) in a product space.

Site-dependent parameters such as the couplings ``g_i`` are made with [`IndexedVariable`](@ref), and pair parameters with [`DoubleIndexedVariable`](@ref).

```@example tc
g = IndexedVariable(:g, i)
```

```@docs
IndexedOperator
IndexedVariable
DoubleIndexedVariable
```

## Sums

With these pieces the Hamiltonian is a single [`Σ`](@ref).

```@example tc
H = Δ * a' * a + Σ(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
```

A sum takes part in the algebra like any other expression. The Heisenberg equation of the cavity field collects the contribution of every atom.

```@example tc
1im * commutator(H, a)
```

```@docs
Σ
```

### Coinciding sites

An index under a sum also runs over the site of every other index on the same subsystem. On that site the two operators act on one atom and compose, so the term is split off and the remaining sum excludes it.

```@example tc
Σ(σ(2, 2, i), i) * σ(2, 1, j)
```

A ground-state projector produced on the shared site stays as it is, as everywhere else.

```@example tc
Σ(σ(1, 2, i), i) * σ(2, 1, j)
```

### Free indices

Two indices that are not summed over may or may not label the same site. Their operators are therefore left in the order they were written, and no rule is applied between them.

```@example tc
σ(2, 2, j) * σ(1, 2, i)
```

When the indices are known to label different sites, [`assume_distinct_index`](@ref) says so, and the product is brought into canonical form.

```@example tc
assume_distinct_index(σ(2, 2, j) * σ(1, 2, i), [(i, j)])
```

```@docs
assume_distinct_index
```

### Renaming indices

An equation derived for one site holds for every other site under a new index name, and a correlation between two sites is obtained from a single-site expression by renaming. [`change_index`](@ref) renames an index throughout and brings the result back into canonical form. The coherence of atom ``i`` becomes that of atom ``j``.

```@example tc
change_index(a' * σ(1, 2, i), i, j)
```

```@docs
change_index
```

### Inspecting indices

The indices of an expression and their properties are read off directly.

```@example tc
get_indices(a' * σ(1, 2, i) * σ(2, 1, j))
```

```@example tc
index_name(i), index_range(i)
```

```@docs
get_indices
has_index
index_name
index_range
index_slot
index_sym
rename
constraint_pairs
```

## Explicit sites

An index treats all of its sites alike. A lattice with nearest-neighbour couplings, open boundaries or a definite number of sites is written instead with one subsystem per site, as in the driven Bose–Hubbard chain of ``M`` coupled cavities

```math
H = \sum_{k=1}^{M} \left[ -\Delta\, a_k^\dagger a_k + \frac{U}{2}\, a_k^\dagger a_k^\dagger a_k a_k + \Omega\,(a_k + a_k^\dagger) \right] - J \sum_{k=1}^{M-1} \left( a_k^\dagger a_{k+1} + a_{k+1}^\dagger a_k \right).
```

Each site gets its own [`FockSpace`](@ref), and both sums are a `sum` over a vector of terms.

```@example chain
using SecondQuantizedAlgebra
@variables Δ U Ω J
M = 6
hc = reduce(⊗, [FockSpace(Symbol(:c, k)) for k in 1:M])
a = [Destroy(hc, Symbol(:a, k), k) for k in 1:M]

onsite = sum([-Δ * a[k]' * a[k] + (U / 2) * a[k]' * a[k]' * a[k] * a[k] + Ω * (a[k] + a[k]') for k in 1:M])
hopping = sum([-J * (a[k]' * a[k+1] + a[k+1]' * a[k]) for k in 1:(M-1)])
H = onsite + hopping
length(H)
```

Summing a vector builds the result in place. A chain of `+` gives the same expression, but each `+` copies the sum accumulated so far, so ``n`` terms cost ``O(n^2)`` instead of ``O(n)``.

```@example chain
terms = [Ω * (a[k] + a[k]') for k in 1:M]
sum(terms) == foldl(+, terms)
```

For a few sites the difference is negligible. For long chains the copies outweigh the work of forming the terms themselves:

| Sites ``M`` | Terms | Forming the terms | `foldl(+, terms)` | `sum(terms)` |
|:-:|:-:|:-:|:-:|:-:|
| 24  | 142  | 0.4 ms | 0.18 ms | 0.008 ms |
| 400 | 2398 | 6.3 ms | 63 ms   | 0.22 ms  |

!!! warning "Use brackets, not a bare generator or a `+=` loop"
    The fast path is taken by `sum`, `reduce(+, …)` and `MutableArithmetics.operate(sum, …)` over an array. A bare generator and an accumulation loop both fall back to repeated `+`:

    ```julia
    H = sum([Ω * (a[k] + a[k]') for k in 1:M])            # fast
    H = sum(Ω * (a[k] + a[k]') for k in 1:M)              # slow: generator
    H = 0; for k in 1:M; H += Ω * (a[k] + a[k]'); end     # slow: copies on every +=
    ```

    Wrapping the generator in brackets, or collecting the terms into a vector first, is enough. The extra vector of length ``M`` is negligible next to the quadratic saving.

`MutableArithmetics.@rewrite` gives the same result, but its macro expansion adds overhead and it does not take the fast path for generators. It suits code that already composes a larger MutableArithmetics expression.

```julia
import MutableArithmetics as MA
fast = sum([Ω * (a[k] + a[k]') for k in 1:M])
fast == MA.operate(sum, [Ω * (a[k] + a[k]') for k in 1:M]) == MA.@rewrite(sum(Ω * (a[k] + a[k]') for k in 1:M))
```

The binary `+`, `-` and scalar `*` copy because the caller still holds both operands, so neither may be modified in place. Only a reduction over a whole collection can build its result in one buffer. The same holds wherever many terms are added up:

- equations of motion and cumulant systems, where each derived ``\frac{d}{dt}\langle O \rangle`` is a sum of many contributions;
- ensembles written with explicit sites, for instance collective operators ``S_\pm = \sum_i \sigma_i^\pm``;
- Lindbladians with one dissipator per channel;
- truncated Magnus, Dyson or BCH expansions.

Products are not affected; their cost lies in the reordering itself. A symbolic [`Σ`](@ref) stays a single sum however many sites it covers, so the question does not arise there.

## Collective systems

Collective operators ``S^{ij} = \sum_k \sigma_k^{ij}`` act on the ensemble as a whole and form a closed algebra without individual site labels. This algebra holds on the full atomic Hilbert space as well as within its permutation-symmetric sector. The numerical example below chooses a basis for that sector; this restriction belongs to the numerical representation, not to the symbolic commutation relations. The two descriptions suit different problems:

| | Indexed sites | Collective operators |
| --- | --- | --- |
| Operators | ``\sum_k \sigma_k^{ij}`` with [`Index`](@ref) and [`Σ`](@ref) | [`CollectiveTransition`](@ref) ``S^{ij}`` |
| Number of atoms | symbolic | fixed by the numerical basis |
| Site-dependent parameters | yes | no |
| Site-resolved moments and correlations | yes | no |
| Numerical basis | full product space for a concrete ``N`` | permutation-symmetric sector in the example below |

They encode different physical roles and live on different Hilbert-space types, so they do not mix on one subsystem.

```@docs
CollectiveNLevelSpace
CollectiveTransition
```

The number of systems enters only through the numerical basis. A bosonic many-body basis over an `NLevelBasis` with ``N`` particles is the permutation-symmetric sector, of dimension ``\binom{N+n-1}{n-1}`` for ``n`` levels, so ``N + 1`` for two.

```@example collective
using SecondQuantizedAlgebra, QuantumOpticsBase
h = CollectiveNLevelSpace(:atoms, (:g, :e))
S(k, l) = CollectiveTransition(h, :S, k, l)

N = 10
b1 = NLevelBasis(2)
b = ManyBodyBasis(b1, bosonstates(b1, N))
length(b)
```

In this basis the collective projectors add up to the number of atoms.

```@example collective
to_numeric(S(:g, :g) + S(:e, :e), b) ≈ N * identityoperator(b)
```

For two levels this is the collective spin of length ``N/2``. With ``e`` as the upper level, ``S^{eg} = S_+``, ``S^{ge} = S_-`` and ``S_z = (S^{ee} - S^{gg})/2``. In the ordering of the `SpinBasis` of QuantumOpticsBase, where level 1 is spin up, the same identification reads ``S^{12} = S_+``.


`get_sum_indices(expr)` reads only the bound summation indices. It returns an empty vector for an unsummed expression.

```@docs
get_sum_indices(::SecondQuantizedAlgebra.QAdd)
```
