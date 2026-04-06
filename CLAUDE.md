# CLAUDE.md — SecondQuantizedAlgebra.jl

## What is this?

A Julia package for symbolic manipulation of second-quantized operators used in quantum many-body theory and quantum optics. Refactored from QuantumCumulants.jl, it provides eagerly-ordered algebraic expressions over bosonic, spin, and N-level operators with automatic commutation relations, configurable ordering conventions, and numeric conversion via QuantumOpticsBase.

## Package layout

```
src/SecondQuantizedAlgebra.jl   # Main module, exports, imports
src/types.jl                    # Abstract type hierarchy: QField, QSym, QTerm, OrderingConvention
src/hilbertspace.jl             # HilbertSpace, FockSpace, ProductSpace, ⊗
src/fock.jl                     # Destroy, Create (bosonic ladder operators)
src/nlevel.jl                   # NLevelSpace, Transition (|i⟩⟨j| operators)
src/pauli.jl                    # PauliSpace, Pauli (σx, σy, σz)
src/spin.jl                     # SpinSpace, Spin (Sx, Sy, Sz with rational spin)
src/phase_space.jl              # PhaseSpace, Position, Momentum
src/qmul.jl                     # QMul{T} — lazy product with c-number prefactor
src/qadd.jl                     # QAdd{T} — sum of QMul terms
src/macros.jl                   # @qnumbers convenience macro
src/interface.jl                # TermInterface/SymbolicUtils integration
src/simplify.jl                 # Worklist simplifier: commutation rules, like-term collection
src/commutator.jl               # commutator(a, b) = ab - ba with short-circuits
src/normal_order.jl             # normal_order(): creation left, annihilation right
src/printing.jl                 # Unicode display (†, σ subscripts, ℋ)
src/latexify_recipes.jl         # LaTeX rendering via Latexify.jl
src/numeric.jl                  # to_numeric(), numeric_average() via QuantumOpticsBase
```

## Type hierarchy

```
QField (abstract)
├── QSym (abstract) — leaf operators
│   ├── Destroy / Create          (FockSpace)
│   ├── Transition                (NLevelSpace)
│   ├── Pauli                     (PauliSpace)
│   ├── Spin                      (SpinSpace)
│   └── Position / Momentum       (PhaseSpace)
└── QTerm (abstract) — compound expressions
    ├── QMul{T<:Number}           (lazy product)
    └── QAdd{T<:Number}           (sum)

HilbertSpace (abstract)
├── FockSpace
├── NLevelSpace
├── PauliSpace
├── SpinSpace
├── PhaseSpace
└── ProductSpace{T<:Tuple}
```

## Key design decisions

- **Lazy evaluation**: multiplication builds `QMul` trees; commutation rules only apply during `simplify()`
- **Space-indexed operators**: each `QSym` carries a `space_index::Int` for composite Hilbert spaces
- **Type-stable prefactors**: `QMul{T}` and `QAdd{T}` are parameterized by the number type; simplification may widen to `Complex{T}` when commutation rules introduce `im`
- **Concrete struct fields**: all struct fields are concretely typed (enforced by CheckConcreteStructs in tests)

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

Each test file is self-contained (`test/*_test.jl`):

**Quality gates:**
- `aqua_test.jl` — Aqua.jl: unbound args, undefined exports, stale deps, piracy
- `jet_test.jl` — JET.jl: type errors and optimization analysis
- `explicit_imports_test.jl` — ExplicitImports.jl: no implicit imports
- `concrete_test.jl` — CheckConcreteStructs: all struct fields concretely typed

**Unit tests:** `fock_test.jl`, `nlevel_test.jl`, `pauli_test.jl`, `spin_test.jl`, `phase_space_test.jl`, `hilbertspace_test.jl`, `qmul_test.jl`, `qadd_test.jl`, `simplify_test.jl`, `commutator_test.jl`, `normal_order_test.jl`, `interface_test.jl`, `macros_test.jl`, `printing_test.jl`, `latexify_test.jl`, `numeric_test.jl`, `types_test.jl`

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
- **`QMul{T}` / `QAdd{T}` parameterized by number type.** Arithmetic preserves or widens `T` as needed.
- **`ProductSpace{T}` uses concrete `Tuple` storage.** The type parameter `T` is a concrete tuple type.

### Performance

- **Type stability first.** All operations should be inferrable — verify with `@inferred`.
- **Minimize allocations in simplification.** The worklist algorithm in `simplify.jl` is the hot path.
- **No kwargs in hot paths.** Keyword arguments prevent specialization and can allocate. Expose kwargs at the API boundary, forward to positional-arg inner functions.
- **No kwargs splatting.** Never forward `kwargs...` — it blocks inference. Explicitly name and forward each keyword.
- **Column-major access.** First index varies fastest. Inner loops over `i` (rows), outer loops over `j` (columns): `for j in 1:n, i in 1:m`.
- **Fuse broadcasts.** Use `@.` or dot syntax to avoid temporary arrays. Use in-place fused assignment: `y .= @. 3x^2 + 4x`.
- **`@views` for slices.** Array slicing copies; use `@view` or `@views` to avoid allocation.
- **`@inbounds` with `eachindex`.** Use `@inbounds` only when indices are provably valid. Prefer `eachindex(x)` over `1:length(x)`.
- **`abs2(z)` over `abs(z)^2`.** Avoids intermediate allocation for complex numbers.
- **Avoid string interpolation in I/O.** Use `println(file, a, " ", b)` not `println(file, "$a $b")`.

### Imports and style

- **No `using X` without explicit imports.** Use `using X: func1, func2` or `import X`. ExplicitImports.jl enforces this.
- **Format with Runic.** Run `make format` before committing.

## Dependencies

| Package | Purpose |
|---------|---------|
| QuantumOpticsBase | Numeric basis types (FockBasis, NLevelBasis, SpinBasis), `⊗`, LazyTensor |
| SymbolicUtils | Symbolic tree traversal interface |
| Symbolics | Symbolic variables (`@variables`) for parameterized prefactors |
| TermInterface | `iscall`, `operation`, `arguments` protocol |
| Latexify | LaTeX rendering recipes |

**Weak dependency:** ModelingToolkit (extension point, not in core)

**Test dependencies:** Aqua, JET, CheckConcreteStructs, ExplicitImports, LaTeXStrings
