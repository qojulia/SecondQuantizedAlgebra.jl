# CLAUDE.md — SecondQuantizedAlgebra.jl

## What is this?

A Julia package for symbolic manipulation of second-quantized operators used in quantum many-body theory and quantum optics. Refactored from QuantumCumulants.jl, it provides eagerly-ordered algebraic expressions over bosonic, spin, and N-level operators with automatic commutation relations, configurable ordering conventions, and numeric conversion via QuantumOpticsBase.

## Package layout

```
src/SecondQuantizedAlgebra.jl   # Main module, exports, imports
src/types.jl                    # Abstract type hierarchy: QField, QSym, QTerm, OrderingConvention
src/cnum.jl                     # CNum = Complex{Num} arithmetic, fast paths, constants
src/hilbertspace.jl             # HilbertSpace, FockSpace, ProductSpace, ⊗
src/fock.jl                     # Destroy, Create (bosonic ladder operators)
src/nlevel.jl                   # NLevelSpace, Transition (|i⟩⟨j| operators)
src/pauli.jl                    # PauliSpace, Pauli (σx, σy, σz)
src/spin.jl                     # SpinSpace, Spin (Sx, Sy, Sz with rational spin)
src/phase_space.jl              # PhaseSpace, Position, Momentum
src/qadd.jl                     # QAdd — dict-based sum, all arithmetic, site sorting
src/ordering.jl                 # Worklist ordering: commutation rules, algebraic reductions
src/simplify.jl                 # simplify() — reductions only (no commutation swaps)
src/normal_order.jl             # normal_order(), normal↔symmetric (Weyl) conversion
src/commutator.jl               # commutator(a, b) = a*b - b*a with index collapse
src/index_types.jl              # Index type, NO_INDEX constant
src/index.jl                    # IndexedVariable, DoubleIndexedVariable, Σ, expand_sums
src/operators.jl                # fundamental_operators, find_operators, unique_ops
src/average.jl                  # AvgFunc/AvgSym, average(), undo_average()
src/substitute.jl               # substitute(expr, Dict) for operators and symbols
src/interface.jl                # TermInterface/SymbolicUtils integration
src/printing.jl                 # Unicode display (†, σ subscripts, ℋ)
src/latexify_recipes.jl         # LaTeX rendering via Latexify.jl
src/numeric.jl                  # to_numeric(), numeric_average() via QuantumOpticsBase
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
└── QTerm (abstract)
    └── QAdd                      (dict-based sum: Dict{Vector{QSym}, CNum})

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

- **Eager ordering**: multiplication immediately applies the global `ORDERING[]` convention and returns a fully-ordered `QAdd`. There is no lazy product type.
- **Dict-based term storage**: `QAdd` stores `Dict{Vector{QSym}, CNum}` — operator sequences map to prefactors. Like terms are collected on construction.
- **CNum prefactors**: prefactors are `Complex{Num}` (from Symbolics.jl), not parameterized. Dedicated fast paths in `cnum.jl` short-circuit for numeric (non-symbolic) cases.
- **Site-indexed operators**: each `QSym` carries `space_index` and `index::Index`. Operators interact only if `_same_site(a, b)`.
- **Separation of ordering vs reduction**: `_apply_ordering` handles commutation swaps (used by `normal_order`), `_apply_reductions` handles algebraic identities only (used by `simplify`).
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

**Unit tests:** `fock_test.jl`, `nlevel_test.jl`, `pauli_test.jl`, `spin_test.jl`, `phase_space_test.jl`, `hilbertspace_test.jl`, `qmul_test.jl`, `qadd_test.jl`, `simplify_test.jl`, `commutator_test.jl`, `normal_order_test.jl`, `interface_test.jl`, `macros_test.jl`, `printing_test.jl`, `latexify_test.jl`, `numeric_test.jl`, `types_test.jl`, `operators_test.jl`, `substitute_test.jl`, `average_test.jl`, `indexing_test.jl`

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

**Test dependencies:** Aqua, JET, CheckConcreteStructs, ExplicitImports, LaTeXStrings
