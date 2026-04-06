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
- [ ] **get_indices for average expressions** — `get_indices` should traverse `BasicSymbolic` trees to extract indices from inner operators (currently throws MethodError on averaged sums).
- [ ] **NLevel completeness relation in simplify** — `simplify(σ_gg)` should reduce ground-state projectors via `Σ_i |i⟩⟨i| = 1` (e.g. `σ₁₁ = 1 - σ₂₂` for 2-level).
- [ ] **undo_average with symbolic prefactors** — `undo_average` on expressions with `@variables` prefactors hits `complex(BasicSymbolic, BasicSymbolic)` MethodError. Fix CNum reconstruction in `average.jl`.
- [ ] **_to_number for BasicSymbolic** — `numeric_average(average(a)^2, ψ)` fails because `^` produces a `BasicSymbolic` exponent. Add `_to_number(::BasicSymbolic)` method in `numeric.jl`.
- [ ] **Σ(constant, i) simplification** — `Σ(g(j), i)` where summand is independent of index `i` should simplify to `N * g(j)`. Currently stays as lazy sum.
- [ ] **Σ(0, i) simplification** — `Σ(0, i)` and `Σ(σ₂₁*σ₂₁, i)` (zero product) should simplify to 0. Currently stays as lazy QAdd wrapping zero.
- [ ] **change_index with identical=false propagation** — `change_index(Ω(i,j), i, j)` where `Ω` has `identical=false` should produce 0. Currently returns `Ω(j,j)` instead.
- [ ] **create_index_arrays** — Utility to generate Cartesian product of index ranges, needed for `expand_sums`. Legacy: `sqa.create_index_arrays([i, j], [1:10, 1:5])`.
- [ ] **conj propagation into indexed variables** — `conj(g_ik)` and adjoint of sums containing indexed variables should correctly conjugate the variable prefactors (legacy: `Σ(conj(g_ik)*a'*σ, i) == Ssum1'`).
- [ ] **LaTeX rendering for indexed sums** — LaTeX output should use `\neq` not `≠` for non-equal index constraints. Legacy: `!contains(latex_str, "≠")`.
- [ ] **one/zero for symbolic variables** — `one(@variables ω::Real)` and `zero(@variables G::Complex)` should work. Legacy: `one(ω) == 1.0`, `zero(G) == 0.0 + 0.0im`.
- [ ] **_conj/_adjoint on symbolic expressions** — `_conj(3*G) == 3*conj(G)`, `_adjoint(exp(im*ϕ)*r) == exp(-im*ϕ)*r` for `@variables` parameters. Partially works but not tested.
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
