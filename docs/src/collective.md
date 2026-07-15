# Collective N-Level Systems

There are two complementary ways to describe many identical N-level systems.
Choose between them based on the approximation and Hilbert space you need.

| Approach | Symbolic operator | Numeric representation | Best suited to |
|---|---|---|---|
| Indexed sites | ``\sum_k \sigma_k^{ij}`` with [`Index`](@ref) and [`Σ`](@ref) | Product-space or cumulant workflow | site-dependent parameters, correlations, and large-``N`` mean-field/cumulant models |
| Collective symmetric sector | [`CollectiveTransition`](@ref) ``S^{ij}`` | `QuantumOpticsBase.ManyBodyBasis` | exact permutation-symmetric dynamics at finite ``N`` |

## Collective transition algebra

`CollectiveTransition` represents

```math
S^{ij}=\sum_{k=1}^{N}|i\rangle_k\langle j|,
```

without expanding the site sum. Its closed Lie algebra is

```math
[S^{ij},S^{kl}]=\delta_{jk}S^{il}-\delta_{li}S^{kj}.
```

```jldoctest collective
julia> using SecondQuantizedAlgebra

julia> h = CollectiveNLevelSpace(:atoms, (:g, :e));

julia> S(i, j) = CollectiveTransition(h, :S, i, j);

julia> commutator(S(:g, :e), S(:e, :g))
S₁₁ - S₂₂
```

Unlike a single-site [`Transition`](@ref), a product such as ``S^{ge}S^{eg}``
does not reduce to a projector: it contains contributions from distinct atoms.
There is likewise no algebra-layer completeness rule. For a fixed-``N``
symmetric numeric basis,

```math
\sum_i S^{ii}=N I,
```

not ``I``.

## Exact symmetric-subspace conversion

Use a bosonic `ManyBodyBasis` whose one-body basis is an `NLevelBasis`. Fixing
the total occupation to ``N`` selects the permutation-symmetric sector, whose
dimension is ``\binom{N+n-1}{n-1}`` for ``n`` internal levels.

```jldoctest collective_numeric
julia> using QuantumOpticsBase

julia> using SecondQuantizedAlgebra

julia> h = CollectiveNLevelSpace(:atoms, (:g, :e));

julia> S(i, j) = CollectiveTransition(h, :S, i, j);

julia> N = 10; b1 = NLevelBasis(2);

julia> b = ManyBodyBasis(b1, bosonstates(b1, N));

julia> size(to_numeric(S(:g, :e), b).data)
(11, 11)

julia> sum(to_numeric(S(x, x), b) for x in (:g, :e)) ≈ N * one(b)
true
```

For two levels this is the usual collective-spin representation with spin
``N/2``: ``S^{12}=S_+``, ``S^{21}=S_-``, and
``(S^{11}-S^{22})/2=S_z`` in the QuantumOpticsBase basis convention.

## Indexed alternative

Keep the site label when atoms are distinguishable or when the calculation
needs site-resolved moments:

```julia
N = 100
h = NLevelSpace(:atom, (:g, :e))
i = Index(h, :i, N, h)
σge = IndexedOperator(Transition(h, :σ, :g, :e), i)
collective_sum = Σ(σge, i)
```

This form and `CollectiveTransition` intentionally do not mix on one subspace:
they encode different physical roles and live on different Hilbert-space types.
