# TODOs

## Ready to implement
- [ ] **SymmetricOrder** — `struct SymmetricOrder <: OrderingConvention` for TWA (Truncated Wigner Approximation). Spec in `src/simplify.jl:195` and `src/types.jl:37`.

## Upstream issues
- [ ] **Printing/LaTeX hooks for custom `BasicSymbolic` nodes** — Two related display issues caused by SymbolicUtils v4 not providing dispatch points for custom `Term`/`AddMul` rendering:
  1. **LaTeX `\langle...\rangle`** — `latexify(average(a))` renders as `\mathrm{avg}(a)` instead of `\langle a \rangle`. Unicode printing works (`⟨a⟩`) via `SymbolicUtils.show_call` hook on our `AvgFunc` type, but the Symbolics LaTeX pipeline (`_toexpr` in `SymbolicsLatexifyExt`) has no equivalent hook.
  2. **Summation `Σ` prefix on averaged sums** — `average(Σ(b'σ_i, i))` prints as `⟨b† * σ_i₁₂⟩` instead of `Σ(i=1:N)⟨b† * σ_i₁₂⟩`. The summation metadata is preserved (recoverable via `get_sum_indices`) but invisible in both unicode and LaTeX display. The legacy code solved this with dedicated `IndexedAverageSum` types that owned their `show` method; we eliminated those types for cleaner algebra.
  
  **Root cause:** `BasicSymbolic{T <: SymVariant}` constrains `T` to SymbolicUtils' own variants, so downstream packages can't dispatch `show`/`latexify` on the type parameter. The only non-piracy hook is `SymbolicUtils.show_call(io, f, x)` for `Term` nodes, which works for unicode but has no LaTeX equivalent.
  
  **Fix:** Request hooks in SymbolicUtils/Symbolics — e.g. a `_toexpr(f::CustomOp, args)` dispatch point for LaTeX, and a metadata-aware `show` that renders `Σ` prefixes for nodes carrying `SumIndices`. Filed as: (TODO: open issue on JuliaSymbolics/Symbolics.jl or JuliaSymbolics/SymbolicUtils.jl).
