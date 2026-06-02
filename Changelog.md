# Changelog

All notable changes to [`SecondQuantizedAlgebra.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Documentation

- Document why alpha-equivalent sums (`Σᵢ σᵢ + Σⱼ σⱼ`) are not auto-collected and clarify the bound-index naming policy in the developer docs; add a "read the devdocs first" note for contributors ([#156]). Resolves [[#134](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/134)].

## [v0.5.0]

Breaking release: redesign of the algebraic core. Most public entry points keep their names and meaning. The substantive changes fall into three groups: direct renames, constructs replaced by more general machinery, and behavioural changes in shared API names. This entry is the migration reference for users with code written against **SecondQuantizedAlgebra.jl v0.4** or **QuantumCumulants.jl**, which shares the v0.4 API surface.

### Renamed

One-to-one renames. Replace the left column with the right column anywhere it appears in your code.

| v0.4 / QuantumCumulants.jl | v0.5 |
|---|---|
| `reorder(expr, [(α, β), …])` | `assume_distinct_index(expr, [(α, β), …])` |
| `simplify(expr, h)` | `expand_completeness(simplify(expr))` |
| `normal_order(expr, h)` | `expand_completeness(normal_order(expr))` |
| `@cnumbers a b c` | `@variables a b c` *(reexported from Symbolics)* |
| `@rnumbers a b c` | `@variables a::Real b::Real c::Real` |
| `Parameter(:x)` | `@variables x` and use `x` |

### Removed

A handful of constructs that lived at the type level in v0.4 have been replaced by mechanisms that scale better and avoid bespoke types.

**`ClusterSpace`.** Indexed families of identical subsystems are now expressed through the `Index` and `Σ` machinery. Where v0.4 wrote

```julia
hc = ClusterSpace(NLevelSpace(:atom, 2), N, k)
σ = Transition(hc, :σ, 1, 2)
H = Δ * sum(σ' * σ)
```

v0.5 writes

```julia
ha = NLevelSpace(:atom, 2)
i  = Index(ha, :i, N, ha)
σ(i) = IndexedOperator(Transition(ha, :σ, 1, 2), i)
H = Δ * Σ(σ(i)' * σ(i), i)
```

`N` stays symbolic through equation derivation. The cumulant truncation order, which used to be the `k` parameter on `ClusterSpace`, is now a property of the moment expansion and is passed to the downstream solver (e.g. `meanfield(…; order=k)` in QuantumCumulants.jl) rather than baked into the Hilbert space.

**`SingleSum` / `DoubleSum` / `SpecialIndexedTerm`.** These wrapper types are gone. `Σ(expr, i)` constructs a `QAdd` directly, with the summation index recorded on the `QAdd.indices` field and per-term inequality constraints recorded on each term's `ne` vector. The `SpecialIndexedTerm` role (carrying `(α, β) ∈ ne` constraints attached to a single product) is filled by the per-term `ne` field. Code that pattern-matched on these types should iterate the resulting `QAdd` and inspect `term.ne` directly.

**`NumberedOperator` / `insert_index`.** v0.5 keeps indices symbolic throughout. Numeric basis instantiation for an indexed family of operators is the responsibility of the downstream numeric layer, not the algebra.

**`IndexedAverageSum` / `IndexedAverageDoubleSum`.** Averaging an indexed `QAdd` now lifts the `indices` and per-term `ne` onto the resulting `BasicSymbolic` as metadata via `SymbolicUtils.setmetadata`. There is no dedicated indexed-average type; the same averaging machinery handles indexed and non-indexed cases uniformly.

**`QMul`.** All arithmetic returns `QAdd`. Code that branched on `QMul`-vs-`QAdd` can drop the branch: multiplication is uniform-return-type. Iteration over a `QAdd` yields `Pair{QTerm, CNum}`; reach into `term.ops` and `term.ne` directly for the operator string and its constraints.

**`order_by_index`.** Operators with symbolic indices on the same Hilbert subspace are sorted automatically when `_partial_sort!` can resolve their relationship (sum bound, declared inequality). When the relationship is undetermined, declare the pairs with `assume_distinct_index` to trigger the sort.

### Changed

A few API names survive but their behaviour has shifted. The differences are intentional and pay off in correctness or performance, so be aware of them when porting tests.

**`σᵍᵍ` no longer expands automatically.** In v0.4, every product passing through `*` rewrote ground-state projectors via `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ`. In v0.5, `σᵍᵍ` is a legitimate canonical-form atom and stays atomic through `*`, `normal_order`, and `simplify`. Use `expand_completeness(expr)` to materialise the identity explicitly. The change preserves like-term collection across operations and avoids an `n_levels^k` blow-up for products containing many ground-state factors.

**`simplify` now includes coefficient simplification.** In v0.4, `simplify` applied algebraic identities to operators only. In v0.5, `simplify` is `normal_order` followed by `Symbolics.simplify` on each surviving coefficient and a drop of summation indices that no surviving term depends on. The expensive symbolic simplification deliberately runs once per surviving term rather than on every internal dict insertion, so the cost shows up at the `simplify` call site rather than inside `*`.

**`Transition` carries `ground_state` and `n_levels`.** The constructor signature `Transition(h, :σ, i, j)` still works and infers these from the surrounding `NLevelSpace`, but the fields now live on the operator itself rather than being looked up on demand. Downstream code that previously had to consult `h` to expand completeness no longer needs the Hilbert space.

**Free indices stay in physical order under `*`.** Two operators carrying different symbolic indices on the same Hilbert subspace, with no `Σ` binding either and no declared inequality, are stored in their physical order rather than being sorted alphabetically by name. Declare the pairs with `assume_distinct_index(q, [(j, k)])` to trigger sorting (and any same-site collapse that follows).

### Migration

Side-by-side, an indexed Tavis-Cummings Hamiltonian.

```julia
# v0.4 / QuantumCumulants.jl
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h  = hc ⊗ ha
@qnumbers a::Destroy(h, 1)
σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
@cnumbers N Δ
g(i) = IndexedVariable(:g, i)
i = Index(h, :i, N, ha)
H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

# v0.5
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h  = hc ⊗ ha
@qnumbers a::Destroy(h, 1)
σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j, 2), k)
@variables N Δ
g(i) = IndexedVariable(:g, i)
i = Index(h, :i, N, ha)
H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
```

The differences are local: `@cnumbers` becomes `@variables`, and `Transition` takes an explicit `n_levels` argument (2 for a two-level atom). Everything downstream (commutators, averaging, equation derivation) uses the same names with the same meaning.

### Unchanged

These names keep their meaning across the migration. Code that only uses them should port without changes:

`FockSpace`, `NLevelSpace`, `PauliSpace`, `SpinSpace`, `PhaseSpace`, `ProductSpace`, `⊗`, `tensor`, `@qnumbers`, `Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`, `Index`, `IndexedOperator`, `IndexedVariable`, `DoubleIndexedVariable`, `Σ`, `∑`, `change_index`, `get_indices`, `commutator`, `anticommutator`, `acts_on`, `find_operators`, `fundamental_operators`, `unique_ops`, `prefactor`, `operators`, `substitute`, `expand`, `normal_order`, `simplify`, `normal_to_symmetric`, `symmetric_to_normal`, `average`, `undo_average`, `is_average`, `to_numeric`, `numeric_average`.


<!-- Links generated by Changelog.jl -->

[v0.5.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.0
[#156]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/156
