# TODOs

## Ready to implement
- [ ] **SymmetricOrder** — `struct SymmetricOrder <: OrderingConvention` for TWA (Truncated Wigner Approximation). Spec in `src/simplify.jl:195` and `src/types.jl:37`.

## Legacy features to port
- [ ] **Average system** — `average()`, `undo_average()`, symbolic `<a'a>` expressions. Legacy: `src/legacy/average.jl`. Related: qojulia/SecondQuantizedAlgebra.jl#57 (double average bug), qojulia/SecondQuantizedAlgebra.jl#61 (change average representation to MTK-compatible time-dependent variables).
- [ ] **Indexed operators** — `Index`, `IndexedOperator`, `SingleSum`, `DoubleSum`, symbolic summations. Legacy: `src/legacy/indexing.jl`, `index_average.jl`, `index_double_sums.jl`.
- [ ] **Cluster spaces** — `ClusterSpace` for mean-field cluster expansions. Legacy: `src/legacy/cluster.jl`.

## Issues from GitHub
- [ ] **expand() and simplify() for QNumbers** — qojulia/SecondQuantizedAlgebra.jl#70. Current impl goes through `average`/`undo_average` round-trip, causes problems.
- [ ] **Operator functions** — qojulia/SecondQuantizedAlgebra.jl#64. Support `exp()`, trig functions of operators (needed for superconducting circuits).
- [ ] **Move to BasicSymbolic** — qojulia/SecondQuantizedAlgebra.jl#60. Use `BasicSymbolic` from SymbolicUtils directly for c-number representation.
- [ ] **Use MTK parameters macro** — qojulia/SecondQuantizedAlgebra.jl#62. Replace `@cnumbers` with `@parameters` from ModelingToolkit once complex number support lands.
- [ ] **Test downstream packages** — qojulia/SecondQuantizedAlgebra.jl#16. CI for QuantumCumulants.jl and other dependents.

## Done (in redesign-v2)
- [x] **Concretely typed structs** — qojulia/SecondQuantizedAlgebra.jl#25.
- [x] **Symbolics v7 / SymbolicUtils v4 migration** — qojulia/SecondQuantizedAlgebra.jl#82.

## Minor
- [ ] **Symbol-level NLevelSpace API** — `NLevelSpace(:atom, (:g,:e))` and `Transition(h, :s, :g, :e)` convenience constructors.
