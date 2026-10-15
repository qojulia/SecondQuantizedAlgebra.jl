# SecondQuantizedAlgebra.jl

## What is this?

A Julia package for symbolic manipulation of second-quantized operators used in quantum many-body theory and quantum optics. Refactored from QuantumCumulants.jl, it provides eagerly-canonicalized algebraic expressions over bosonic, spin, and N-level operators with automatic commutation relations and numeric conversion via QuantumOpticsBase.

## Read the architecture page first

**Always read [`docs/src/dev/architecture.md`](docs/src/dev/architecture.md) before changing anything under `src/` or `ext/`.** It lists the invariants a change must keep, each with its reason and the test that pins it: eager canonical form, the pipeline order, why `σᵍᵍ` stays atomic, the naming policy (the algebra never mints `Index` objects on the user's behalf), diagonal splitting, the dead-NE rule, the `QAdd`/`QExpr` split, and the coefficient tiers. Several of them are not covered by tests, and many "obvious" changes are deliberately not made for reasons recorded there. The developer manual in `docs/src/dev/` also covers the workflow (`contributing.md`), adding an operator kind or a numeric backend (`extending.md`), and the internal docstrings (`internals.md`).

## Package layout

```
src/SecondQuantizedAlgebra.jl       # Main module, exports, imports
src/types.jl                        # QField, QSym abstract hierarchy, SiteCmp, ReduceKind
src/precompile.jl                   # PrecompileTools workload
src/average.jl                      # AvgFunc/SumFunc nodes, average(), undo_average(), Hermitian-Real symtype

src/numeric/backend.jl              # NumericBackend + singletons; open hooks; default backend
src/numeric/coeff.jl                # coefficient folding, parameter/time splitting
src/numeric/core.jl                 # NumericContext and the static/time-dependent conversion core
src/numeric/indexed.jl              # backend-neutral indexed unroll (sites, ne)
src/numeric/api.jl                  # public to_numeric / numeric_average / expect
src/numeric/qexpr.jl                # deliberate numeric refusal for unlowered QExpr
src/numeric/superoperator.jl        # superoperator conversion and to_numeric_master

ext/SecondQuantizedAlgebraQuantumOpticsBaseExt.jl  # QuantumOpticsBase backend (vector LazySum)
ext/SecondQuantizedAlgebraQuantumToolboxExt.jl     # QuantumToolbox backend (VecSum over QobjEvo)

src/expressions/intern.jl           # name and range interning tables, name ranks
src/expressions/gaussian.jl         # GaussianRational{T}: exact (re + im*i)/den numbers, checked Int kernels
src/expressions/exact.jl            # exact scalar tiers, the NativeSlot layout, radical-folding primitives
src/expressions/monomial.jl         # Monomial{E}/Poly structs (radical-canonical by construction)
src/expressions/cnum.jl             # Coeff (CNum) arithmetic, recognition, raw lowering, LinearAlgebra hooks
src/expressions/reduce.jl           # trigonometric and hyperbolic coefficient reduction
src/expressions/qterm.jl            # QTerm struct (ops, ne), dict key for QAdd
src/expressions/qadd.jl             # QAdd, canonical polynomial sum; TermInterface
src/expressions/qexpr.jl            # QExpr, cold-path formal non-polynomial expressions
src/expressions/qexpr_structural.jl # structural rewrites/introspection for QExpr
src/expressions/qexpr_taylor.jl     # explicit exact-coefficient Taylor lowering to QAdd
src/expressions/index_types.jl      # Index type, NO_INDEX constant
src/expressions/index.jl            # IndexedVariable, DoubleIndexedVariable, Σ

src/operators/hilbertspace.jl       # HilbertSpace, FockSpace, ProductSpace, ⊗
src/operators/op.jl                 # OpKind enum, Op struct, hash/isequal, is_* predicates, optype
src/operators/fock.jl               # Destroy, Create constructors
src/operators/nlevel.jl             # NLevelSpace, Transition, CollectiveNLevelSpace, CollectiveTransition
src/operators/pauli.jl              # PauliSpace, Pauli constructor
src/operators/spin.jl               # SpinSpace, Spin constructor
src/operators/phase_space.jl        # PhaseSpace, Position, Momentum constructors
src/operators/operators.jl          # the operator hooks (kind-branched), order_key, adjoint, qadjoint

src/algebra/algebra.jl              # normal_order(), simplify(), expand_completeness(), commutator(), substitute()
src/algebra/passes.jl               # partial_sort! and the reduce/commute/ground-state passes
src/algebra/pipelines.jl            # stream!, canonicalize!, products and diagonal splitting
src/algebra/mutable_arithmetics.jl  # in-place sum/reduce accumulation (QAddBuilder)
src/algebra/weyl.jl                 # normal_to_symmetric(), symmetric_to_normal()
src/algebra/unitary*.jl             # UnitaryTransform: affine actions, constructors, generators, frames, QExpr recursion

src/printing/printing.jl            # Unicode display
src/printing/latexify_recipes.jl    # LaTeX rendering via Latexify.jl
src/printing/qexpr.jl               # precedence-aware text/LaTeX rendering for QExpr
src/printing/superoperator.jl       # rendering of maps and Lindbladians

src/superoperators/                 # SuperTerm, SuperOperator, Channel, Lindbladian, SuperExpr, scope, frames
```

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

Tests are organized by behavioral contract rather than by mirroring `src/`. The detailed taxonomy and placement rule live in [`test/README.md`](test/README.md): put a test where a contributor would look when the behavior breaks, not where the implementing function currently lives.

- `test/fundamentals/`: Hilbert spaces, operator construction, indexing primitives, names, and canonical ordering.
- `test/algebra/`: arithmetic, commutators, normal ordering, simplification, substitution, and unitary transformations.
- `test/symbolic/`: coefficients, phases, averages, indexed constraints, time dependence, and symbolic reductions.
- `test/numeric/`: core conversion, backend interfaces, materialization/performance, expectation values, parameter/operator substitution, indexed conversion, and backend-specific behavior.
- `test/presentation/`: text/Unicode and LaTeX rendering.
- `test/integration/`: workflows and invariants that deliberately cross subsystem boundaries.
- `test/quality/`: static and repository-wide quality gates.

Directories encode concepts, filenames encode contracts, and nested `@testset`s encode scenarios. Tests that intentionally span several concepts belong in `integration/`. The suite is discovered recursively by `test/runtests.jl` using ParallelTestRunner; keep the `quality/JET` path stable because the runner handles it specially.

**Quality gates** (`test/quality/quality_test.jl`): Aqua (unbound args, stale deps, piracy), JET (type errors), ExplicitImports, CheckConcreteStructs.

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

### Imports and style

- **No `using X` without explicit imports.** Use `using X: func1, func2` or `import X`. ExplicitImports.jl enforces this.
- **Use ordinary names for all files and bindings.** Do not create underscore-prefixed filenames, functions, or other bindings; Julia's module system handles visibility.
- **Format with Runic.** Run `make format` before committing.

### Documentation and comments

- **Comments: compact, why not what.** Default to no comment; add one only for a non-obvious *why*. Keep it to a couple of lines.
- **Changelog entries: short.** One to two sentences per entry. State what changed and where it matters, not the mechanism.
- **Developer manual: invariants and recipes, not private-function tours.** State the rule and its reason in a sentence or two and point to the test that pins it. Do not enumerate private functions or struct fields; those belong in internal docstrings.
-
## Performance terminology

Three distinct axes; keep them separate when reporting numbers.

- **Precompilation time**: compiling package code into the on-disk cache. Paid once per version/dependency change, not per session; where `precompile.jl`'s `@compile_workload` runs.
- **TTFX**: cold-process latency to the *first* call: `using` load time plus any JIT not cached during precompilation. A bigger `@compile_workload` trades more precompilation time for less TTFX.
- **Runtime**: steady-state cost of every call after the first, no compilation. What `@btime`/`benchmark/` measure after warmup; target of the zero-allocation passes.

## Dependencies

| Package | Purpose |
|---------|---------|
| Combinatorics | Product enumeration for operator generation |
| QuantumInterface | Lightweight owner of `expect`/`basis` and the `Basis`/`AbstractOperator`/`StateVector` types (hard dep) |
| TensorCore | Owner of the shared `⊗`/`tensor` generic function, imported directly since QuantumInterface 0.4.4 (hard dep) |
| SymbolicUtils | Symbolic tree traversal interface |
| Symbolics | Symbolic variables (`@variables`), `Num` type for scalar symbolic boundaries |
| TermInterface | `iscall`, `operation`, `arguments` protocol |
| Latexify | LaTeX rendering recipes |
| LinearAlgebra | `det` on a `Matrix{Coeff}` via minor expansion; stdlib, already in the load closure |
| PrecompileTools | `@setup_workload`/`@compile_workload` in `precompile.jl` |
| Primes | `isprime` for radicand cofactors above `2^32`; already loaded by Symbolics |

**Weak dependencies (numeric extensions):** `QuantumOpticsBase` (FockBasis/NLevelBasis/SpinBasis, vector `LazySum`/`TimeDependentSum`), `QuantumToolbox` + `SciMLOperators` (QuantumObject builders, the `VecSum` over `QobjEvo`). The numeric API errors until one backend is loaded.

**Test dependencies:** Aqua, JET, CheckConcreteStructs, ExplicitImports, LaTeXStrings, QuantumOptics (master-equation solver checks), QuantumOpticsBase, QuantumToolbox, SciMLOperators, LinearAlgebra