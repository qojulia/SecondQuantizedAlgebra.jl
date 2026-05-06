# TODOs

## Codebase improvements

Independently-shippable items from the post–PR-#99 codebase audit. Each is a discrete
PR with no ripple effect; suggested merge order is roughly the order listed.

### Top priority

- [ ] **Drop the type-piratical `expect` overloads on `BasicSymbolic`/`Num`**
      ([numeric.jl:271](src/numeric.jl#L271)). The file itself flags the piracy and
      `aqua_test.jl` suppresses it via `treat_as_own = [expect]`. Either delete the
      methods (callers use `numeric_average` directly — same length) or wrap averaged
      operands behind an SQA-owned type. Recommendation: **delete**. Removes four
      methods and the Aqua suppression.

### Worth doing if you're touching that area

- [x] **Unify the `numeric_average` signatures** ([numeric.jl](src/numeric.jl)).
      Eight methods cover `(QField | BasicSymbolic | Number | Num) × (state |
      states) × (with d | without d)` — many are thin pass-throughs. Collapse to
      two entry points (`numeric_average(op, state, d=Dict())` and the `states`
      vector overload) plus dispatch internals.

- [x] **`change_index` short-circuit on no-op substitutions**
      ([index.jl:82](src/index.jl#L82)). For prefactors that don't reference the
      `from` index — most of them — `Symbolics.substitute` still builds a dict
      and walks the tree. A cheap `Symbolics.get_variables(x)` membership check
      upfront skips the common case. Measurable on commutator-heavy workloads.

- [x] **Tail-recursive worklist in `_apply_ordering`**
      ([ordering.jl:144](src/ordering.jl#L144)). Today every rule firing pushes
      back via `copy(ops)` — fine for small products, but multi-level
      Schrieffer–Wolff–style products allocate hundreds of intermediate terms.
      Mutate in place and recurse on the leftover. Same logical complexity, fewer
      allocations.

- [ ] **Lazy `Index.sym` field** ([index_types.jl](src/index_types.jl)). Every
      `Index` constructor allocates a `Num` for `sym`; only `change_index`
      substitutions actually need it. Compute on demand via a small cache. Modest
      gain — only worth it if `Index` construction shows up in profiles.

### Architectural — discuss before doing

- [ ] **Reconsider eager ground-state expansion**
      ([ordering.jl:108](src/ordering.jl#L108)). `_push_ground_state_expansion!`
      pushes `n_levels` terms per σᵍᵍ reduction, multiplied across multiple σᵍᵍs
      in one product — exponential by construction. Two paths:
      1. **Defer GS expansion** to a single post-pass after all ordering finishes,
         so like-term collection trims duplicates as they form.
      2. **Keep `σ_gg` as a first-class operator**; add an explicit
         `apply_completeness` pass. The dev docs already describe expansion as a
         canonical-form *choice* — making it opt-in (the way `LazyOrder` already
         is for the user) would be a cleaner separation.

      Worth a design discussion before any code moves.

## Documentation

- [ ] **Update `devdocs.md` "Diagonal splitting" section**
      ([docs/src/devdocs.md:179](docs/src/devdocs.md#L179)). The text currently
      describes pre-sort substitution for `sum_idx` pairs only. The rule is more
      general: *any distinct-index pair on the same Hilbert subspace whose
      pre-sort vs post-sort substitution disagree gets a constraint and a
      corrected diagonal at multiplication time.* State as one rule.

- [ ] **Add a "what a dict entry means" paragraph to `devdocs.md`**. Once stated,
      the rest of the design follows from it: *a dict entry `[ops] => (c, ne)`
      represents the term `c·ops` valid for any index assignment satisfying the
      pairwise constraints in `ne`.*

- [ ] **Update the `QAdd` field listing in `devdocs.md`** to reflect that
      `non_equal` is no longer a separate field — constraints live per-term in
      the dict value.

- [ ] **Add a perf note on `LazyOrder`** for users building large symbolic
      expressions where they only need canonicalization at the end.

## Anti-patterns — document as "don't do"

These are tempting answers to specific bugs that would undo the rewrite's design.
Worth recording so the next reader doesn't re-litigate them.

- **Don't reintroduce per-term operator structure** (a `QMul`-like compound
  type). Even if it solves a class of bugs, it undoes the single-`QAdd`,
  type-stable-`*` simplification at the heart of the rewrite. The Patch-1
  (PR #99) approach handles the same bugs without splitting the type hierarchy.

- **Don't specialize `*` per `QSym` subtype** (`*(::Destroy, ::Create)` etc.).
  The dispatch table grows N², and the rules already live in `_apply_ordering`
  cleanly.

- **Don't add hash invariants on `Vector{QSym}` keys**. Dict already does the
  right thing; custom hash invariants are silent-corruption surface.

## Upstream issues
- [ ] **Printing/LaTeX hooks for custom `BasicSymbolic` nodes** — Two related display issues caused by SymbolicUtils v4 not providing dispatch points for custom `Term`/`AddMul` rendering:
  1. **LaTeX `\langle...\rangle`** — `latexify(average(a))` renders as `\mathrm{avg}(a)` instead of `\langle a \rangle`. Unicode printing works (`⟨a⟩`) via `SymbolicUtils.show_call` hook on our `AvgFunc` type, but the Symbolics LaTeX pipeline (`_toexpr` in `SymbolicsLatexifyExt`) has no equivalent hook.
  2. **Summation `Σ` prefix on averaged sums** — `average(Σ(b'σ_i, i))` prints as `⟨b† * σ_i₁₂⟩` instead of `Σ(i=1:N)⟨b† * σ_i₁₂⟩`. The summation metadata is preserved (recoverable via `get_sum_indices`) but invisible in both unicode and LaTeX display. The legacy code solved this with dedicated `IndexedAverageSum` types that owned their `show` method; we eliminated those types for cleaner algebra.
  
  **Root cause:** `BasicSymbolic{T <: SymVariant}` constrains `T` to SymbolicUtils' own variants, so downstream packages can't dispatch `show`/`latexify` on the type parameter. The only non-piracy hook is `SymbolicUtils.show_call(io, f, x)` for `Term` nodes, which works for unicode but has no LaTeX equivalent.
  
  **Fix:** Request hooks in SymbolicUtils/Symbolics — e.g. a `_toexpr(f::CustomOp, args)` dispatch point for LaTeX, and a metadata-aware `show` that renders `Σ` prefixes for nodes carrying `SumIndices`. Filed as: (TODO: open issue on JuliaSymbolics/Symbolics.jl or JuliaSymbolics/SymbolicUtils.jl).
