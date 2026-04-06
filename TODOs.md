# TODOs

## Ready to implement
- [ ] **SymmetricOrder** — `struct SymmetricOrder <: OrderingConvention` for TWA (Truncated Wigner Approximation). Spec in `src/simplify.jl:195` and `src/types.jl:37`.
- [ ] **Dict-based QAdd** — Replace `Vector{QMul}` storage with `Dict{Vector{QSym}, CNum}` for auto like-term collection, order-independent equality, and elimination of `_collect_like_terms`. Resolves audit issues #2 (fragile Dict keys), #4 (order-dependent equality), #5 (`_ZERO_QADD` mutable singleton), #10 (redundant Dict allocation in simplify). Plan: `docs/superpowers/plans/2026-04-05-dict-based-qadd.md`.

## Legacy features to port
- [x] **Average system** — `average()`, `undo_average()`, symbolic `<a'a>` expressions. Implemented in `src/average.jl`. `IndexedAverageSum`/`IndexedAverageDoubleSum` eliminated by redesign (summation metadata on SymbolicUtils nodes). LaTeX rendering pending upstream hook (see Upstream issues).
- [x] **Indexed operators** — `Index`, `IndexedOperator`, `IndexedVariable`, `DoubleIndexedVariable`, `Σ`/`∑`, `change_index`, `expand_sums`, `get_indices`, commutator diagonal collapse. All operator types carry `index::Index` field. Tested in `test/indexing_test.jl`.
- [x] **Cluster spaces** — `ClusterSpace` for mean-field cluster expansions. Legacy: `src/legacy/cluster.jl`.
- [ ] **insert_index** — Substitute a concrete integer for an `Index` to create numbered operators (e.g. `insert_index(σ(2,2,ind(:j)), ind(:j), 1)` → `σ(2,2,1)`). Needed for expanding sums into explicit terms and for indexed numeric conversion.
- [ ] **to_numeric with `ranges` kwarg** — Convert indexed/numbered operators to numeric form on composite bases built from repeated single-site bases (e.g. `to_numeric(σ(1,2,1), basis; ranges=[N_atoms, N_modes])`). Essential for multi-atom simulations (superradiant pulse example).
- [ ] **substitute** — Symbolic parameter substitution in operator expressions (e.g. `substitute(x*a, Dict(x=>0))`, `substitute(H, Dict(N=>10))`). Needed for plugging in parameter values before numeric evaluation.
- [ ] **LazyKet support in numeric_average** — `numeric_average` should work with `QuantumOpticsBase.LazyKet` (lazy tensor product states), not just regular `Ket`.
- [ ] **expand_sums diagonal splitting** — Richer `expand_sums` that splits nested indexed sums into diagonal (`i==j`) and off-diagonal (`i≠j`) terms, sum multiplication with explicit index (`*(sum1, sum2; ind=j)`), and commutators with sums involving symbolic N. Legacy `DoubleSum`/`SingleSum` behavior.
- [ ] **Symbolic NLevel levels with `level_map`** — `NLevelSpace(:atom, (:g,:e,:a))` with `to_numeric(σ(:e,:g), basis; level_map=Dict(:g=>1, :e=>2, :a=>3))`.
- [x] **get_indices for average expressions** — Fixed: added `get_indices(::BasicSymbolic)` in `average.jl` that traverses symbolic trees, unwrapping `AvgFunc` nodes.
- [ ] **NLevel completeness relation in simplify** — `simplify(σ_gg)` should reduce ground-state projectors via `Σ_i |i⟩⟨i| = 1` (e.g. `σ₁₁ = 1 - σ₂₂` for 2-level). Blocked: `Transition` doesn't carry a reference to its `NLevelSpace` (no `n`/`ground_state` info). Needs either fields on `Transition` or a context-aware `simplify(expr, h::HilbertSpace)` API.
- [x] **undo_average with symbolic prefactors** — Fixed: only recurse into `+`/`*` nodes in `undo_average`, convert other `BasicSymbolic` c-number nodes (e.g. `complex(0, η)`) to `CNum` directly. Also added `_to_cnum(::BasicSymbolic)` to split `complex(a,b)` nodes.
- [x] **_to_number for BasicSymbolic** — Fixed: added `_to_number(::BasicSymbolic)` in `numeric.jl` that unwraps `Const` nodes and converts via `Num`.
- [x] **Σ(constant, i) simplification** — Fixed: `_depends_on_index(::QMul, ::Index)` checks both operators and prefactor for index dependence. `Σ` multiplies by `range` when summand is independent.
- [x] **numeric_average product-of-averages** — Fixed: `numeric_average` now checks `operation isa AvgFunc` instead of `is_average` (which matched products of averages via symtype).
- [x] **Σ(0, i) simplification** — `Σ(0, i)` now simplifies via `_depends_on_index` (scalar 0 is index-independent → `range * 0 = 0`). `Σ(σ₂₁*σ₂₁, i)` still needs normal_order first (lazy product, not zero until simplified).
- [ ] **change_index with identical=false propagation** — `change_index(Ω(i,j), i, j)` where `Ω` has `identical=false` should produce 0. Currently returns `Ω(j,j)` instead.
- [ ] **create_index_arrays** — Utility to generate Cartesian product of index ranges, needed for `expand_sums`. Legacy: `sqa.create_index_arrays([i, j], [1:10, 1:5])`.
- [x] **conj propagation into indexed variables** — `conj(g_i)` works (Real-typed `IndexedVariable` is identity under `conj`). `adjoint` of sums with indexed variable prefactors works. Tests added in `indexing_test.jl`.
- [x] **LaTeX rendering for indexed sums** — Fixed `\ne` → `\neq` in `latexify_recipes.jl`. Test added in `latexify_test.jl`.
- [x] **one/zero for symbolic variables** — `one(@variables ω::Real)` and `zero(@variables G::Complex)` already work via Symbolics. Tests added in `operators_test.jl`.
- [x] **_conj/_adjoint on symbolic expressions** — `_conj(3*G) == 3*conj(G)` works. Tests added in `operators_test.jl`. Note: `exp(im*ϕ)` case fails structurally (`cos(ϕ)-im*sin(ϕ)` vs `cos(-ϕ)+im*sin(-ϕ)`) — Symbolics limitation, not ours.
- [ ] **acts_on for indexed operators** — `acts_on(IndexedOperator(σ, i))` should return the space index consistently. Legacy returned `Int`, current returns `Vector{Int}`.

## Upstream issues
- [ ] **Printing/LaTeX hooks for custom `BasicSymbolic` nodes** — Two related display issues caused by SymbolicUtils v4 not providing dispatch points for custom `Term`/`AddMul` rendering:
  1. **LaTeX `\langle...\rangle`** — `latexify(average(a))` renders as `\mathrm{avg}(a)` instead of `\langle a \rangle`. Unicode printing works (`⟨a⟩`) via `SymbolicUtils.show_call` hook on our `AvgFunc` type, but the Symbolics LaTeX pipeline (`_toexpr` in `SymbolicsLatexifyExt`) has no equivalent hook.
  2. **Summation `Σ` prefix on averaged sums** — `average(Σ(b'σ_i, i))` prints as `⟨b† * σ_i₁₂⟩` instead of `Σ(i=1:N)⟨b† * σ_i₁₂⟩`. The summation metadata is preserved (recoverable via `get_sum_indices`) but invisible in both unicode and LaTeX display. The legacy code solved this with dedicated `IndexedAverageSum` types that owned their `show` method; we eliminated those types for cleaner algebra.
  
  **Root cause:** `BasicSymbolic{T <: SymVariant}` constrains `T` to SymbolicUtils' own variants, so downstream packages can't dispatch `show`/`latexify` on the type parameter. The only non-piracy hook is `SymbolicUtils.show_call(io, f, x)` for `Term` nodes, which works for unicode but has no LaTeX equivalent.
  
  **Fix:** Request hooks in SymbolicUtils/Symbolics — e.g. a `_toexpr(f::CustomOp, args)` dispatch point for LaTeX, and a metadata-aware `show` that renders `Σ` prefixes for nodes carrying `SumIndices`. Filed as: (TODO: open issue on JuliaSymbolics/Symbolics.jl or JuliaSymbolics/SymbolicUtils.jl).

## Issues from GitHub
- [x] **expand() and simplify() for QNumbers** — qojulia/SecondQuantizedAlgebra.jl#70. Current impl goes through `average`/`undo_average` round-trip, causes problems.
- [ ] **Operator functions** — qojulia/SecondQuantizedAlgebra.jl#64. Support `exp()`, trig functions of operators (needed for superconducting circuits).
- [x] **Move to BasicSymbolic** — qojulia/SecondQuantizedAlgebra.jl#60. Use `BasicSymbolic` from SymbolicUtils directly for c-number representation.
- [x] **Use MTK parameters macro** — qojulia/SecondQuantizedAlgebra.jl#62. Replace `@cnumbers` with `@parameters` from ModelingToolkit once complex number support lands.

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
