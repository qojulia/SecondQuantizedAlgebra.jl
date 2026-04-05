# TODOs

## Ready to implement
- [ ] **SymmetricOrder** — `struct SymmetricOrder <: OrderingConvention` for TWA (Truncated Wigner Approximation). Spec in `src/simplify.jl:195` and `src/types.jl:37`.
- [ ] **Dict-based QAdd** — Replace `Vector{QMul}` storage with `Dict{Vector{QSym}, CNum}` for auto like-term collection, order-independent equality, and elimination of `_collect_like_terms`. Resolves audit issues #2 (fragile Dict keys), #4 (order-dependent equality), #5 (`_ZERO_QADD` mutable singleton), #10 (redundant Dict allocation in simplify). Plan: `docs/superpowers/plans/2026-04-05-dict-based-qadd.md`.

## Legacy features to port
- [ ] **Average system** — `average()`, `undo_average()`, symbolic `<a'a>` expressions. Legacy: `src/legacy/average.jl`. Related: qojulia/SecondQuantizedAlgebra.jl#57 (double average bug), qojulia/SecondQuantizedAlgebra.jl#61 (change average representation to MTK-compatible time-dependent variables).
- [x] **Indexed operators** — `Index`, `IndexedOperator`, `IndexedVariable`, `DoubleIndexedVariable`, `Σ`/`∑`, `change_index`, `expand_sums`, `get_indices`, commutator diagonal collapse. All operator types carry `index::Index` field. Tested in `test/indexing_test.jl`.
- [x] **Cluster spaces** — `ClusterSpace` for mean-field cluster expansions. Legacy: `src/legacy/cluster.jl`.

## QC integration hooks (SQA-side complete, QC must build on top)
- [ ] **Averaging types** — QC needs `IndexedAverageSum`, `IndexedAverageDoubleSum` wrapping `average()` on `QAdd` with indices. SQA exports: `Index`, `has_index`, `QAdd.indices`, `QAdd.non_equal`.
- [ ] **NumberedOperator** — QC needs concrete integer-indexed operators for numeric evaluation (`σ_i` → `σ₃` → matrix). Lives in QC, not SQA.
- [ ] **`scale()` orchestration** — QC calls `expand_sums()` for diagonal splitting, then does QC-specific post-processing (N-dependent factors, averaging).
- [ ] **`value_map`** — QC substitutes concrete index values for numeric solving.

## Issues from GitHub
- [ ] **expand() and simplify() for QNumbers** — qojulia/SecondQuantizedAlgebra.jl#70. Current impl goes through `average`/`undo_average` round-trip, causes problems.
- [ ] **Operator functions** — qojulia/SecondQuantizedAlgebra.jl#64. Support `exp()`, trig functions of operators (needed for superconducting circuits).
- [ ] **Move to BasicSymbolic** — qojulia/SecondQuantizedAlgebra.jl#60. Use `BasicSymbolic` from SymbolicUtils directly for c-number representation.
- [ ] **Use MTK parameters macro** — qojulia/SecondQuantizedAlgebra.jl#62. Replace `@cnumbers` with `@parameters` from ModelingToolkit once complex number support lands.
- [ ] **Test downstream packages** — qojulia/SecondQuantizedAlgebra.jl#16. CI for QuantumCumulants.jl and other dependents.

## Done (in redesign-v2)
- [x] **Concretely typed structs** — qojulia/SecondQuantizedAlgebra.jl#25.
- [x] **Symbolics v7 / SymbolicUtils v4 migration** — qojulia/SecondQuantizedAlgebra.jl#82.
- [x] **`QMul`/`QAdd` use `Complex{Num}` prefactors** — dropped type parameters, always use `CNum = Complex{Num}`.
- [x] **`simplify` extends `SymbolicUtils.simplify`** — public API delegates to internal `_qsimplify`.
- [x] **Export pruning** — minimal public API, no internal types exported.
- [x] **Indexed operators** — ground-up redesign with `index::Index` field on all `QSym` subtypes.

## Minor
- [ ] **Symbol-level NLevelSpace API** — `NLevelSpace(:atom, (:g,:e))` and `Transition(h, :s, :g, :e)` convenience constructors.
- [x] **`order_by_index`** — canonical ordering of operators by index name within `QMul.args_nc` via `_sort_key(op) = (space_index, copy_index, index.name)`.
