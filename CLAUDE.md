# CLAUDE.md — SecondQuantizedAlgebra.jl

## What is this?

A Julia package for symbolic manipulation of second-quantized operators used in quantum many-body theory and quantum optics. Refactored from QuantumCumulants.jl, it provides eagerly-canonicalized algebraic expressions over bosonic, spin, and N-level operators with automatic commutation relations and numeric conversion via QuantumOpticsBase.

## Read the devdocs first

**Always read [`docs/src/devdocs.md`](docs/src/devdocs.md) before making any change to the repo.** It documents the internal architecture and, crucially, the *design rationale*: why eager canonicalization works the way it does, why `σᵍᵍ` stays atomic, the naming policy (the algebra never mints `Index` objects on the user's behalf), diagonal splitting, and the dead-NE invariant. Many "obvious" changes are deliberately not made for reasons recorded there. Check it before proposing or implementing anything.

## Package layout

```
src/SecondQuantizedAlgebra.jl       # Main module, exports, imports
src/types.jl                        # QField, QSym abstract hierarchy
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
src/operators/operators.jl          # fundamental_operators, find_operators, unique_ops, qadjoint

src/algebra/algebra.jl              # normal_order(), simplify(), expand_completeness(), commutator(), substitute()
src/algebra/passes.jl               # _partial_sort!, _reduce_ops, _commute_ops, _expand_gs_ops
src/algebra/pipelines.jl            # _stream!, _canonicalize!, _emit_product!, _accumulate_with_diag!
src/algebra/weyl.jl                 # normal_to_symmetric(), symmetric_to_normal()

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
```

## Key design decisions

- **Eager canonicalization**: every `*` immediately runs the full pipeline — sort operators to canonical order, reduce same-site pairs algebraically (`_reduce_ops`), commute non-canonical same-site pairs (`_commute_ops`), reduce again to catch any composition the commute residual produced. The result of `a * a'` is already `a'*a + 1`. `normal_order` and `simplify` are idempotent on expressions built via `*`.
- **Pipeline order is `reduce → commute → reduce`**: the first reduce handles Transition/Pauli same-site composition (these compose, not commute); commute handles Fock/Spin/PhaseSpace ladder pairs; the trailing reduce catches residuals. The Spin commutator emits a contracted-axis operator (e.g. `[Sy, Sx] = -i Sz`) via a uniform 4-tuple `(swap_b, swap_a, residual_coeff, residual_ops)`.
- **`σᵍᵍ` is a canonical atom**: ground-state projectors stay atomic through `*`, `normal_order`, and `simplify`. The completeness identity `σᵍᵍ = 1 - Σ_{k≠g} σᵏᵏ` does not fire automatically. Call `expand_completeness(expr)` explicitly when you want it. This prevents exponential term blowup in dissipator-style expressions.
- **`Transition` carries its own GS info**: `ground_state::Int` and `n_levels::Int` are fields on `Transition`, keeping operators Hilbert-space-decoupled.
- **`assume_distinct_index(q, pairs)`**: explicit escape hatch when two free indices semantically denote distinct sites but no `Σ` supplies the constraint. Takes a `Vector{Tuple{Index, Index}}` of inequality pairs, augments each term's `ne`, re-canonicalizes, and runs `expand_completeness`.
- **Free indices outside `Σ` stay `Undetermined`**: two operators with different symbolic indices on the same space, neither bound by a sum, are left in physical order. No same-site collapse fires until `assume_distinct_index` or a `Σ`-driven diagonal split resolves the relationship.
- **Dict-based term storage**: `QAdd` stores `Dict{QTerm, CNum}` where `QTerm` bundles `ops::Vector{QSym}` with `ne::Vector{NonEqualPair}` index-inequality scope. Like terms are collected on construction.
- **CNum prefactors**: prefactors are `Complex{Num}` (from Symbolics.jl), not parameterized. Dedicated fast paths in `cnum.jl` short-circuit for numeric (non-symbolic) cases.
- **Site-indexed operators**: each `QSym` carries `space_index` and `index::Index`. Operators interact only if `_same_site(a, b)`.
- **Five operator hooks**: every concrete `QSym` subtype implements `_site_compare`, `_can_commute`, `_commute_pair`, `_reduce_pair`, `_ground_state_expand`. The algebra talks to all operator types exclusively through these. Adding a new operator type only requires defining the struct and filling in the five hooks.
- **Concrete struct fields**: all struct fields are concretely typed (enforced by CheckConcreteStructs in tests).
- **Index tracking**: `QAdd` carries `indices::Vector{Index}` for summation scope; per-term inequality constraints live on `QTerm.ne`, not on `QAdd`.

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

Tests are organized in subdirectories matching `src/`, run via `test/runtests.jl` (ParallelTestRunner).

**Quality gates** (`test/quality/quality_test.jl`): Aqua (unbound args, stale deps, piracy), JET (type errors), ExplicitImports, CheckConcreteStructs.

**Operator tests** (`test/operators/`): `fock_test.jl`, `nlevel_test.jl`, `pauli_test.jl`, `spin_test.jl`, `phase_space_test.jl`, `hilbertspace_test.jl`, `operators_test.jl`

**Expression tests** (`test/expressions/`): `algebra_test.jl`, `indexing_test.jl`, `macros_test.jl`

**Arithmetic tests** (`test/arithmetics/`): `simplify_test.jl`, `commutator_test.jl`, `normal_order_test.jl`, `substitute_test.jl`, `internals_test.jl`

**Other:** `test/printing/printing_test.jl`, `test/average_test.jl`, `test/numeric_test.jl`

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
- **Minimize allocations in the pipeline.** The streaming passes in `passes.jl` are the hot path. Sink type-parameters (`F` in `where {F}`) force specialization so nested `do`-blocks inline with zero closure allocation.
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
| PrecompileTools | `@setup_workload`/`@compile_workload` in `precompile.jl` |

**Test dependencies:** Aqua, JET, CheckConcreteStructs, ExplicitImports, LaTeXStrings
