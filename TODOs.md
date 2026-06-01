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
- [ ] migration guide

### Refactors

- [ ] **Hash-cons operator leaves**. Intern each `QSym` value (`Destroy`, `Create`,
      `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`) so equality is `===`.
      Concrete wins: `hash(::QTerm)` becomes O(n) on identity; `==` on `Vector{QSym}`
      is pointer comparison. Implementation: `Dict{Tuple, QSym}` cache, weakly held
      if leak risk materializes. Constructors go through helpers; user-facing API
      unchanged.

- [ ] **Lazy `Commutator` head**. Today `commutator(a, b)` returns `a*b - b*a`
      immediately. A `Commutator(left, right)` head plus `expand_commutator`
      distributor (`[AB,C] = A[B,C] + [A,C]B`) keeps the algebraic structure visible
      until expansion is asked for. Useful for Heisenberg-picture work and for
      keeping intermediate expressions compact.

- [ ] **Specialize `ne` for mergeable operators only**. Structural insight from the
      [QC.jl #230](https://github.com/qojulia/QuantumCumulants.jl/pull/230) analysis:
      `ne` is only meaningful when same-space indexed operators can *merge*, and
      merging only happens for `Transition`, not for `Destroy`/`Create`. Fast-path
      the case where neither term contains a mergeable operator in
      `_accumulate_with_diag!` and skip diagonal-split entirely. Less invasive than
      switching to the full UUID approach below.

- [ ] **`SmallDict` for short `QAdd` sums**. For sums of ≤ 8 terms, a linear-probe
      fixed-buffer dict avoids `Dict`'s constant-factor allocation. Add only after
      profile shows it.

- [ ] **Cache `uses_phys_key` as a `QTerm` field**
      ([src/expressions/qterm.jl:28](src/expressions/qterm.jl#L28)). Today
      `_uses_phys_key(term)` re-runs an O(n²) `phys_ops × phys_ops` scan against
      `ne` on every `isequal` / `hash`, and `QTerm` is the dict key for every
      `QAdd` so this fires on every insertion and lookup. Cache the result as a
      `uses_phys_key::Bool` field set once in `_term_key` (after `_canonical_ne`)
      and have `isequal` / `hash` read the field instead. Tried briefly and
      reverted — re-attempt with a benchmark to confirm the win and to check that
      adding the field doesn't regress construction allocations elsewhere.

- [ ] **Benchmark CI gates**. Suite exists in [benchmark/](benchmark/), regressions
      aren't gated. Add a CI step that fails on regression against tracked baselines
      with a documented "accept regression" override.

### New features

- [ ] **Vacuum expectation value `vev`**. `vev(expr) = scalar_term(normal_order(expr))`.
      First user-facing observable; composes with everything else. Tests:
      `vev(a^n * dag(a)^n) == factorial(n)` for n=1..5; `vev(dag(a)*a) == 0`.

- [ ] **Coherent state expectations**. `expect(expr, CoherentState; alphas)` =
      `normal_order` then substitute `a → α`, `a† → α*`. Multi-mode trivial.

- [ ] **JSON serialization**. Stable schema with top-level `sqa_version`. Round-trip
      `read_json(write_json(expr)) == expr` after `simplify`. Useful for sharing
      Hamiltonians and snapshot tests. Bump version only on schema breakage; reader
      supports prior versions.

- [ ] **Unitary transformation operators**. Opaque heads for `Displacement(α, mode)`,
      `Squeezing(ξ, mode)`, `Rotation(φ, mode)`, `Quadrature(:x/:p, mode)`. Each
      carries its group law and conjugation transformation; no eager BCH expansion.
      Quadrature gets an `expand_quadrature` lowering to `(a±a†)/√2` so
      `latex(x)` can read as `\hat{x}`.

- [ ] **Per-call kwargs on `simplify`** (replaces the LazyOrder workflows). If a
      user needs "keep `|g⟩⟨g|` atomic" or "reductions only" without resurrecting
      `LazyOrder`, add as kwargs: `simplify(expr; expand_completeness=false)` and
      `simplify(expr; commute=false)`. Local control, no global mode. Add only
      when a concrete user need arises.

- [ ] **Bogoliubov transformation**

- [ ] **Index arithmetic / nearest-neighbour shifts**. Today an `Index` is a
      single `name::Symbol`, and structural equality on `name` is what powers
      `_diagonal_split!`, `_site_compare`, and the `NonEqualPair` constraint
      system. Lattice Hamiltonians of the form ``H = -J \sum_i (\sigma_i^+
      \sigma_{i+1}^- + \mathrm{h.c.})`` therefore can't be written natively —
      `i+1` has no representation.
      **MVP:** add a `shift::Int` (or `step`) field to `Index` so
      `Index(name=:i, shift=1)` denotes ``i+1``. Touches:
      - [src/expressions/index_types.jl](src/expressions/index_types.jl) —
        struct, `==`, `hash`, `isless` extended on `(name, shift)`.
      - [src/expressions/index.jl](src/expressions/index.jl) —
        `change_index`, `get_indices`, and crucially
        [`_diagonal_split!`](src/expressions/index.jl#L252): currently
        collides indices by `name`; with shifts, decide when `(i, 0)` and
        `(i, 1)` ever coincide (no, by construction) and when `(i, 1)`
        colliding with free `j` on the same subspace should emit the
        substitution `j → i+1`.
      - Per-operator `_site_compare` (one per `QSym` subtype) — sort key
        becomes `(space_index, name, shift)` so canonical ordering stays
        deterministic.
      - Numeric `Σ` enumeration — enumerate shifted indices respecting
        boundary conditions: open chain (`i+1 ≤ N`) vs periodic
        (`(i mod N) + 1`). Pass the BC via a `boundary=:open|:periodic`
        keyword on `Σ`.
      **Full lattice support** (out of scope for the MVP) would also want:
      sparse-pattern `DoubleIndexedVariable` like `J(i, j)` for arbitrary
      coupling matrices, a generic `j = f(i)` for higher-shell shifts and
      multi-orbital unit cells, and explicit boundary-condition objects.
      **Motivation:** unblocks an indexed Heisenberg/XXZ chain example
      (spin-wave dispersion ``\varepsilon(k) = J(1 - \cos k)`` via
      Holstein-Primakoff). Without it the chain example has to either
      unroll for finite `N` (loses the `Σ` showcase) or be skipped.

### Speculative explorations

- [ ] **UUID-based `i ≠ j` tracking** (per [QC.jl #230](https://github.com/qojulia/QuantumCumulants.jl/pull/230)).
      Replaces `ne::Vector{NonEqualPair}` with `merge_events::Vector{UUID}` on
      indexed `Transition`s; materializes `(i==j)` and `(i!=j)` as CNum factors.
      Deletes `_accumulate_with_diag!` (~200 LOC) and the `phys_ops` shadow field.
      **Conflicts with hash-consing** (each merge mints a fresh UUID, defeating
      leaf interning) and the PR is stale since 2025-01-31 with author-flagged
      open issues around `find_missing` / `change_index` / `complete()`. Revisit
      only if `_accumulate_with_diag!` becomes a measured hot bottleneck *and*
      hash-consing isn't shipping.

- [ ] **SoP (`Vector{OperatorTerm}`) vs Dict storage benchmark**. Consider whether
      a vector primary representation outperforms `Dict{QTerm, CNum}` for QO-sized
      sums (≤ 20 terms). Profile first; switching is a non-trivial rewrite of
      `QAdd` and every `_addto!` site.

- [ ] **`(a†a)^n` compaction during normal ordering**. Eager `*` expands fully;
      for large `n` this is `O(n!)` terms. A Stirling-number identity could keep
      number-operator powers folded. Add only if a benchmark forces it.

## Upstream issues
- [ ] **LaTeX `\langle...\rangle` for averaged expressions.** `latexify(average(a))` renders as `\mathrm{avg}(a)` instead of `\langle a \rangle`, and averaged sums lose their `Σ` prefix in LaTeX output. Unicode rendering is fully wired up (`⟨a⟩` via `SymbolicUtils.show_call` on `AvgFunc`; `Σ(i=1:N)⟨...⟩` via `SymbolicUtils.show_metadata` on `SumIndices`, see `src/printing/printing.jl`). The Symbolics LaTeX pipeline (`_toexpr` in `SymbolicsLatexifyExt`) has no equivalent hook yet. Blocked on JuliaSymbolics/Symbolics.jl#1835 (metadata-aware `_toexpr` dispatch point). Once that lands, mirror the unicode hook with a LaTeX version that emits `\langle ... \rangle` plus the `Σ` prefix from `SumIndices` / `SumNonEqual` metadata.
