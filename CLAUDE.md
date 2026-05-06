# CLAUDE.md — SecondQuantizedAlgebra.jl

## What is this?

A Julia package for symbolic manipulation of second-quantized operators used in quantum many-body theory and quantum optics. Refactored from QuantumCumulants.jl, it provides eagerly-ordered algebraic expressions over bosonic, spin, and N-level operators with automatic commutation relations, configurable ordering conventions, and numeric conversion via QuantumOpticsBase.

## Package layout

```
src/SecondQuantizedAlgebra.jl       # Main module, exports, imports
src/types.jl                        # QField, QSym, OrderingConvention abstract hierarchy
src/precompile.jl                   # PrecompileTools workload
src/average.jl                      # AvgFunc/AvgSym, average(), undo_average()
src/numeric.jl                      # to_numeric(), numeric_average() via QuantumOpticsBase

src/expressions/cnum.jl             # CNum = Complex{Num} arithmetic, fast paths, constants
src/expressions/qterm.jl            # QTerm struct (ops, ne) — dict key for QAdd
src/expressions/qadd.jl             # QAdd — dict-based sum (Dict{QTerm, CNum}); TermInterface
src/expressions/index_types.jl      # Index type, NO_INDEX constant
src/expressions/index.jl            # IndexedVariable, DoubleIndexedVariable, Σ

src/operators/hilbertspace.jl       # HilbertSpace, FockSpace, ProductSpace, ⊗
src/operators/fock.jl               # Destroy, Create (bosonic ladder operators)
src/operators/nlevel.jl             # NLevelSpace, Transition (|i⟩⟨j| operators)
src/operators/pauli.jl              # PauliSpace, Pauli (σx, σy, σz)
src/operators/spin.jl               # SpinSpace, Spin (Sx, Sy, Sz)
src/operators/phase_space.jl        # PhaseSpace, Position, Momentum
src/operators/operators.jl          # fundamental_operators, find_operators, unique_ops

src/arithmetics/ordering.jl         # ScopedValue ORDERING; with_ordering / set_ordering!
src/arithmetics/qadd_arithmetic.jl  # +, -, *, ^ for QAdd; site sorting; worklist commutation
src/arithmetics/simplify.jl         # simplify() — reductions only (no commutation swaps)
src/arithmetics/normal_order.jl     # normal_order(), normal↔symmetric (Weyl) conversion
src/arithmetics/commutator.jl       # commutator(a, b) = a*b - b*a with index collapse
src/arithmetics/substitute.jl       # substitute(expr, Dict) for operators and symbols

src/printing/printing.jl            # Unicode display (†, σ subscripts, ℋ)
src/printing/latexify_recipes.jl    # LaTeX rendering via Latexify.jl
```

## Type hierarchy

```
QField (abstract)
├── QSym (abstract) — leaf operators (each carries space_index, index)
│   ├── Destroy / Create          (FockSpace)
│   ├── Transition                (NLevelSpace)
│   ├── Pauli                     (PauliSpace)
│   ├── Spin                      (SpinSpace)
│   └── Position / Momentum       (PhaseSpace)
└── QAdd                          (dict-based sum: Dict{QTerm, CNum})

QTerm (struct, dict key)          # fields: ops::Vector{QSym}, ne::Vector{NonEqualPair}

HilbertSpace (abstract)
├── FockSpace
├── NLevelSpace
├── PauliSpace
├── SpinSpace
├── PhaseSpace
└── ProductSpace{T<:Tuple}

OrderingConvention (abstract)
├── NormalOrder                   (creation left, annihilation right)
└── LazyOrder                     (no reordering)
```

## Key design decisions

- **Eager ordering produces canonical form**: under `NormalOrder` (default), `*` immediately applies all three transformations — commutation swaps, algebraic reductions, and `NLevelSpace` completeness — so the result is always in the canonical basis. There is no lazy product type.
- **Canonical basis for `NLevelSpace`**: the basis is `{σⁱʲ : (i,j) ≠ (g,g)} ∪ {1}` — ground-state projectors `σᵍᵍ` are *not* part of canonical form. Same-site composition that would produce `σᵍᵍ` (e.g. `σ¹² · σ²¹`) is eagerly rewritten via completeness `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ`. This makes dict-key equality match physical equality.
- **`Transition` carries its own GS info**: `ground_state::Int` and `n_levels::Int` are fields on `Transition`, not arguments to the algebra. The arithmetic does completeness without ever consulting the Hilbert space, keeping operators Hilbert-space-decoupled.
- **`LazyOrder` is the explicit-control opt-out**: under `LazyOrder`, none of the three transformations fire eagerly. Users invoke them piecewise: `simplify(expr)` for reductions only (keeps `σᵍᵍ` atomic), `normal_order(expr)` adds commutation swaps, the `(expr, h)` overloads add completeness. The `h` argument is the LazyOrder opt-in for completeness — it is a no-op under `NormalOrder`.
- **Dict-based term storage**: `QAdd` stores `Dict{QTerm, CNum}` where `QTerm` is a struct bundling `ops::Vector{QSym}` with its `ne::Vector{NonEqualPair}` index-inequality scope. Like terms are collected on construction.
- **Ordering is task-local via ScopedValues**: the active `OrderingConvention` lives in a `ScopedValue` (`ORDERING` in `arithmetics/ordering.jl`). Use `with_ordering(LazyOrder()) do ... end` for transient/task-local scoping; use `set_ordering!(ord)` only to change the process-wide default. Read it via `get_ordering()`.
- **CNum prefactors**: prefactors are `Complex{Num}` (from Symbolics.jl), not parameterized. Dedicated fast paths in `cnum.jl` short-circuit for numeric (non-symbolic) cases.
- **Site-indexed operators**: each `QSym` carries `space_index` and `index::Index`. Operators interact only if `_same_site(a, b)`.
- **Three-layer separation**: `_apply_ordering` handles commutation swaps (used by `normal_order`), `_apply_reductions` handles algebraic identities only (used by `simplify`), and completeness is gated inside `_try_algebraic_reduction!` via the `expand_gs::Bool` parameter (NormalOrder passes `true`, simplify passes `false`) plus `_apply_ground_state` in `normal_order.jl` for the LazyOrder opt-in path.
- **Concrete struct fields**: all struct fields are concretely typed (enforced by CheckConcreteStructs in tests).
- **Index tracking**: `QAdd` carries `indices::Vector{Index}` and `non_equal::Vector{Tuple{Index,Index}}` for summation metadata.

## Git policy

**Never commit or push.** Neither Claude nor any subagent may run `git commit`, `git push`, or any git command that modifies history. All commits are made by the user.

## Development workflow

Common tasks via the Makefile:

```sh
make test       # run all tests (Pkg.test)
make format     # format with Runic (src/ test/ benchmark/)
make docs       # build documentation
make servedocs  # serve docs locally with LiveServer
make bench      # run benchmarks
make benchlocal # run benchmarks, save to data/, print changelog
make all        # setup + format + test + docs
```

### Quick debugging with Julia MCP

Use the `julia-mcp` MCP server (tools: `julia_eval`, `julia_list_sessions`, `julia_restart`) for quick debugging and testing small snippets — e.g., checking a type, evaluating an expression, or verifying a method signature. Prefer this over spinning up a full test run when you just need a quick answer.

### Test structure

Each test file is self-contained (`test/*_test.jl`), run via `test/runtests.jl` (ParallelTestRunner).

**Quality gates:**
- `aqua_test.jl` — Aqua.jl: unbound args, undefined exports, stale deps, piracy
- `jet_test.jl` — JET.jl: type errors and optimization analysis
- `explicit_imports_test.jl` — ExplicitImports.jl: no implicit imports
- `concrete_test.jl` — CheckConcreteStructs: all struct fields concretely typed

**Unit tests:** `fock_test.jl`, `nlevel_test.jl`, `pauli_test.jl`, `spin_test.jl`, `phase_space_test.jl`, `hilbertspace_test.jl`, `qmul_test.jl`, `qadd_test.jl`, `canonical_form_test.jl`, `simplify_test.jl`, `commutator_test.jl`, `normal_order_test.jl`, `macros_test.jl`, `printing_test.jl`, `latexify_test.jl`, `numeric_test.jl`, `operators_test.jl`, `substitute_test.jl`, `average_test.jl`, `indexing_test.jl`

**Integration:** `integration_test.jl` — end-to-end scenarios (Jaynes-Cummings, multi-mode cavities)

### Testing patterns

- `@inferred` for type stability checks
- `@allocations` for zero-allocation verification on hot paths
- Symbolic prefactor tests with `Symbolics.@variables`

### Quality gates

Before merging any PR:
1. `make test` passes (all quality + unit + integration tests)
2. JET reports zero issues
3. No `Any`-typed fields in any struct (CheckConcreteStructs)
4. All imports explicit (ExplicitImports)

## Coding rules

### Function signatures

- **Use the most restrictive signature type possible.** Tight type declarations let JET catch errors. When prototyping it's fine to start loose, but committed code should have specific types.
- **Explicit `;` for keyword arguments.** Always use an explicit semicolon:
  ```julia
  # Good
  FockSpace(; name = :a)
  # Bad
  FockSpace(name = :a)
  ```

### Type system

- **No abstract-typed fields.** Every struct field must be concretely typed.
- **`ProductSpace{T}` uses concrete `Tuple` storage.** The type parameter `T` is a concrete tuple type.

### Performance

- **Type stability first.** All operations should be inferable — verify with `@inferred`.
- **Minimize allocations in ordering.** The worklist algorithm in `ordering.jl` is the hot path.
- **No kwargs in hot paths.** Keyword arguments prevent specialization and can allocate. Expose kwargs at the API boundary, forward to positional-arg inner functions.
- **No kwargs splatting.** Never forward `kwargs...` — it blocks inference. Explicitly name and forward each keyword.

### Imports and style

- **No `using X` without explicit imports.** Use `using X: func1, func2` or `import X`. ExplicitImports.jl enforces this.
- **Format with Runic.** Run `make format` before committing.

## Dependencies

| Package | Purpose |
|---------|---------|
| Combinatorics | Product enumeration for operator generation |
| QuantumOpticsBase | Numeric basis types (FockBasis, NLevelBasis, SpinBasis), `⊗`, LazyTensor |
| SymbolicUtils | Symbolic tree traversal interface |
| Symbolics | Symbolic variables (`@variables`), `Num` type for CNum prefactors |
| TermInterface | `iscall`, `operation`, `arguments` protocol |
| Latexify | LaTeX rendering recipes |
| ScopedValues | Task-local `ORDERING` for `with_ordering`/`get_ordering` |
| PrecompileTools | `@setup_workload`/`@compile_workload` in `precompile.jl` |

**Test dependencies:** Aqua, JET, CheckConcreteStructs, ExplicitImports, LaTeXStrings
