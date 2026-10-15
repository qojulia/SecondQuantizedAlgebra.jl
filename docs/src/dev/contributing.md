# Contributing

Contributions follow the [ColPrac](https://github.com/SciML/ColPrac) guidelines for pull requests, review and releases. Before changing code under `src/`, read the [Architecture](@ref) page: several of its invariants are not covered by tests.

```@contents
Pages = ["contributing.md"]
Depth = 2:3
```

## Workflow

Common tasks are targets of the `Makefile`:

| Command | Effect |
| :-- | :-- |
| `make setup` | install the formatter Runic and the docs tools |
| `make test` | run the whole test suite |
| `make jet` | run the JET analysis on its own |
| `make format` | format all Julia files with Runic |
| `make docs` | build the documentation |
| `make servedocs` | serve the documentation locally with live reload |
| `make bench` | run the benchmarks |

A subset of the tests runs by passing path filters to the test runner, for instance `Pkg.test(; test_args = ["algebra/simplify"])`.

## Tests

The tests are organized by behavior, not by source file: a test goes where a contributor would look when the behavior breaks. The directories and the placement rule are in [`test/README.md`](https://github.com/qojulia/SecondQuantizedAlgebra.jl/blob/main/test/README.md). Every bug fix comes with a test that fails without it.

A pull request passes the quality gates in `test/quality/`:

- **Aqua**: no unbound type parameters, stale dependencies or type piracy;
- **JET**: no type errors or runtime dispatch in the analyzed entry points;
- **ExplicitImports**: every imported name is listed explicitly;
- **CheckConcreteStructs**: every struct field has a concrete type.

The continuous integration also runs the tests of the packages that build on this one, QuantumCumulants.jl, QuantumInputOutput.jl and FloquetExpansions.jl, against the change.

## Code style

- Use the most specific argument types possible, and an explicit `;` before keyword arguments.
- Keep keyword arguments out of hot paths: accept them at the public function and forward positional arguments.
- Check type stability with `@inferred` and allocations with `@allocations` instead of guessing.
- Import names explicitly, `using X: f, g`.
- Name files and functions without a leading underscore; the module controls visibility.
- Add a comment only for a reason that the code cannot show.

Docstrings of public functions are user documentation: a one-sentence summary, every calling form, an example and a See also line. Implementation details belong in the developer documentation or in docstrings of internal functions.

## Changelog

Every change to user-facing code gets an entry in `Changelog.md` under the unreleased version, which the continuous integration checks: one or two sentences on what changed and where it matters, not on how. A breaking change goes under *Changed (breaking)*.
