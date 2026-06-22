# Changelog

All notable changes to [`SecondQuantizedAlgebra.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.8.0]

### Changed

- **Breaking: the seven operator types collapsed into one concrete `Op`.** `Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, and `Momentum` were separate `QSym` subtypes; they are now a single concrete `struct Op <: QSym` carrying a `kind::OpKind` tag plus shared packed fields. The role names stay as constructor functions, so every existing `Destroy(h, :a)` or `Transition(h, :σ, 1, 2)` call still builds the right operator with no source change. What breaks is type-level dispatch: code that branched on `isa Destroy` (or used `::Transition` method signatures, or `typeof(op) == Pauli`) must switch to the exported predicates `is_destroy`, `is_create`, `is_transition`, `is_pauli`, `is_spin`, `is_position`, `is_momentum`, or read the role with `optype(op)::OpKind`. Direct field reads also moved: a `Transition`'s `i`/`j`/`ground_state`/`n_levels` are now `l1`/`l2`/`g`/`nlev`, and a `Pauli`/`Spin` `axis` is `l1`. `Op` and the predicates are exported. Because operator products are now stored as a concrete `Vector{Op}` (`QTerm.ops`), the per-operator hooks (`_site_compare`, `_can_commute`, `_commute_pair`, `_reduce_pair`, `hash`) dispatch and inline statically instead of through the former abstract-`QSym` dynamic dispatch, and `_commute_pair`/`_reduce_pair` return concrete `Tuple{Op, …}` with no tuple boxing; custom `hash`/`isequal` on `Op` exclude the index's `Num` fields. Measured speedups over the v0.7.1 baseline from the operator collapse alone: Fock `(a·a†)^10` about 2.7×, single-mode `H^4` about 2.1×, many-mode `H^2` (M=8) about 1.7×, with roughly 30% fewer allocations on the Fock-dominated cases ([#164](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/164)).
- Parameter-polynomial coefficient arithmetic (the `Poly` tail of `Coeff`) allocates less: single-term times single-term (the common `ω·J` case) skips the sort/compact pass, term canonicalisation compacts its owned buffer in place instead of allocating a second vector, and a scalar-times-term product shares the other operand's factor vectors instead of re-merging. Results are byte-identical (canonical form unchanged) and the arithmetic stays fully type-stable. This cuts about 18% of allocations on many-mode `H^2` (M=8) and 11% on single-mode `H^4`, bringing the combined speedup over v0.7.1 to about 1.9× and 2.2× on those cases ([#164](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/164)).

## [v0.7.2]

### Changed

- Complex parameters and irreducible scalar couplings now stay on the native `Poly` coefficient path instead of escalating to the `Complex{Num}` symbolic tail, where every later product/sum round-tripped through SymbolicUtils hashconsing. `_is_atom`/`_rec` now recognize any non-algebraic one-argument call on an atom (`real(g)`, `imag(g)`, `exp`, `sin`, ...) as a single opaque polynomial atom; a `@variables g::Complex` parameter (stored as `Complex(real(g), imag(g))`) therefore no longer poisons downstream coefficient arithmetic. Conjugation of a `Poly` coefficient is now native (`_conj_poly`: conjugate scalars and atoms, re-canonicalize) rather than materializing each term through `to_num`. On the downstream `QuantumInputOutput.jl` SLH cascade benchmark this lowers a 3-cavity `hamiltonian`/`lindblad` build from about 140 µs to about 7.6 µs (roughly 18×).
- Radicals of a single parameter now canonicalize on the native `Poly` path: `Monomial` exponents are `Rational{Int}`, so `sqrt(p)` is stored as `p^(1//2)` (likewise `cbrt(p)` as `p^(1//3)`, and `p^(m//n)`). Building `sqrt(p)` and squaring it folds to `p` natively, matching the `sqrt(p)^2 -> p` fold Symbolics performs at the `Num` level; previously `sqrt(p)` was an opaque integer-exponent atom, so the two construction paths produced distinct, non-`isequal` coefficients and `QAdd` like-term dedup could miss merges. A radical of a non-atom (`sqrt(g*κ)`, `(g+κ)^(1//2)`) is not distributed (unsound in general) and stays a `Complex{Num}` symbolic leaf, consistent with how Symbolics treats it.

## [v0.7.1]

### changed

- The per-term operator machinery is faster on product and power workloads, lowering the floor that remained after the `Coeff` rework. `QTerm` now caches its hash in a field computed once at construction, so the repeated dict probes, inserts, and growth-rehashes that key every canonicalization no longer re-walk the operator vector (whose abstract `QSym` eltype forces a dynamic dispatch per element); `isequal` also short-circuits on the cached hash before comparing operators. The reduce pass gained a type-domain `_may_reduce` gate consulted before `_reduce_pair`: only same-type `Transition` and `Pauli` pairs can compose, so every Fock, Spin, and `PhaseSpace` pair (the bulk of bosonic products) now skips the dynamically dispatched `_reduce_pair` call whose `(kind, op, factor)` tuple would otherwise box on each adjacent pair. Measured speedups over the post-`Coeff` baseline: Fock `(a·a†)^10` about 1.74×, many-mode `H^2` (M=8) about 1.81×, single-mode `H^4` about 1.54×, with roughly half the allocations ([#141](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/141), [#164](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/164)).

## [v0.7.0]

### Fixed

- Averaging an indexed sum whose coefficient depends on the summation index no longer drops the scope or lets the index dangle outside the `Σ`. `average` previously stamped `SumIndices`/`SumNonEqual` metadata on the average leaf and emitted `coeff * leaf`; SymbolicUtils discards metadata on composite/numerically-scaled nodes, so `average(Σ(u(kk,k), k))` collapsed to a bare `u(kk,k)` and an index-dependent coefficient was hoisted outside the sum. Index-dependent terms are now wrapped in a dedicated moment-layer node (`SumFunc`/`sym_sum`) carrying the summation scope as a `SumScope` argument, so the whole averaged body stays inside the sum and the representation survives `Add`/`Mul` canonicalization. Because the scope rides as an argument (not metadata, which `isequal`/`hash` ignore), differently-scoped sums over the same body no longer wrongly cancel in a subtraction. New `is_indexed_sum` predicate; `has_sum_metadata`/`get_sum_indices`/`get_sum_non_equal` retained and now read the node ([#175](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/175)).

### Changed

- `make_time_dependent` on an averaged indexed sum now yields `Σ(i) ⟨a_i⟩(t)` (per-site time-dependent moments under the sum) instead of `⟨Σ_i a_i⟩(t)` (one collective lumped variable), matching the non-lifted display `Σ(i) ⟨a_i⟩` and giving indexable per-site unknowns for indexed equations.
- Operator prefactors are stored as a concrete `Coeff` with three forms (a native `ComplexF64` fast path, a sparse parameter polynomial for products and sums of named parameters, and a `Complex{Num}` fallback) instead of always `Complex{Num}`. Numeric and parameter-polynomial coefficient arithmetic stays native and never routes through SymbolicUtils hashconsing; a coefficient lowers to `Complex{Num}` only at the symbolic boundaries (`substitute`/`average`/printing/`prefactor`). The polynomial arithmetic is fully type-stable, with factor identity via `objectid`/`===` (which assumes SymbolicUtils hashconsing is enabled, the default). Polynomial coefficients are kept in canonical expanded form, so `(g+h)^2` is stored as `g^2 + 2*g*h + h^2`. Measured speedups over `Complex{Num}`: numeric power expansion about 2.1×, single-mode `H^4` about 2.65×, many-mode `H^2` about 3.5×, nested commutator about 2.4× ([#164](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/164), [#183](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/183)).


## [v0.6.5]

### Added

- `@variables η::Number` declares a complex coefficient that stays a single atomic symbol, conjugating to the symbolic `conj(η)`. Unlike `::Complex` (which splits into independent real/imaginary unknowns and pays `(a+bi)(c+di)` expansion on every product), `::Number` keeps coefficient arithmetic on one symbol (`η * η → η²`). Adjoint routes coefficients through a symtype-aware `_conj_cnum`, so `(η a)†(η a)` correctly carries `|η|² = conj(η)·η`, and `to_numeric` reduces such coefficients (e.g. an unfolded `conj` of a complex literal) via `build_function`.

## [v0.6.4]

### Fixed

- `change_index` now correctly zeros `DoubleIndexedVariable` nodes tagged `identical=false` even when they appear nested inside a larger product or sum. Previously `_check_not_identical` only examined the root of the substituted expression, so a `J(i,j)` factor inside `J(i,j) * x` would survive a `j → i` substitution intact instead of collapsing to zero. The fix walks the full expression tree, collects every `NotIdentical`-tagged node whose two index arguments became equal, and replaces them all via a single `Symbolics.substitute` pass.

## [v0.6.3]

### Changed

- Lifted time-dependent averages now display as `⟨op⟩(t)` in both the REPL and `latexify` output, instead of exposing the internal `_avg_…` variable name. The `AverageOperator` metadata drives the rendering through new `show_metadata`/`_toexpr_metadata` hooks, so a moment's MTK unknown reads the same as its `average` ([#173](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/173)).

## [v0.6.2]

### Added

- `index_slot` (public): recover the concrete position `k` of a per-slot index symbol minted by `(i::Index)(k)`, read from metadata rather than parsed out of the symbol's name. `(i::Index)(k)` now stamps `k` as metadata; since symbol equality and hashing ignore metadata, a per-slot index still dedup-equals its name-only counterpart ([#170](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/170)).

## [v0.6.1]

### Added

- `order_key`, `term_order_key`, `qadd_order_key` (public): a total, identity-faithful structural ordering for operators, products, and sums, built from one per-type `order_key` method. The key ties two operators exactly when they are `isequal`, giving downstream packages a reproducible way to pick canonical representatives and compare expressions without round-tripping through `show` or `hash` ([#169](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/169)).

## [v0.6.0]

### Changed

- Average expressions now use `Number` as their Symbolics `symtype`, and the new `make_time_dependent` helper lifts them into ModelingToolkit-style time-dependent unknowns while preserving average metadata for round-tripping.

## [v0.5.2]

### Fixed

- Commuting commutators with symbolic rational coefficients no longer bloat derived expressions ([#162](https://github.com/qojulia/SecondQuantizedAlgebra.jl/pull/162)). Deriving `[Hₖ, op]` for an `Hₖ` that commutes with `op` should vanish, but with a rational coupling such as `γ/|xᵢ-xⱼ|³` the two halves land as `γ/D + (-γ)/D`, which Symbolics leaves un-combined and `_iszero_cnum` does not see as zero; the spurious term then survives and inflates every downstream RHS. Fixed at two levels:
  - `_addto_key!` cancels exact-negation prefactor pairs for symbolic coefficients (the numeric path stays allocation-free).
  - `commutator` distributes over terms and skips disjoint-subspace pairs (`[aₖ, bₗ] = 0`), so the cancelling products are never formed.

## [v0.5.1]


### Fixed

- `Σ` now **propagates index-inequality constraints onto the diagonal** when a summed index collapses onto an external one. When splitting the diagonal of `Σ_{i≠j} σᵢ²¹ σₖ¹²`, the `i = k` contribution now inherits `k ≠ j` instead of dropping it. Previously the constraint was discarded, letting a later sum re-admit the `k = j` point and double-count the diagonal of nested and ollective double sums ([#161](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/161)).

### Changed

- LaTeX rendering for averaged expressions now uses `\langle \cdots \rangle`, and averaged indexed sums preserve their `\sum` prefix in `latexify` output via the new Symbolics LaTeX hooks ([#153](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/153)).

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
[v0.5.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.1
[v0.5.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.5.2
[v0.6.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.0
[v0.6.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.1
[v0.6.2]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.2
[v0.6.3]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.3
[v0.6.4]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.4
[v0.6.5]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.6.5
[v0.7.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.7.0
[v0.7.1]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.7.1
[v0.8.0]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/releases/tag/v0.8.0
[#156]: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/156
