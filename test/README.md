# Test suite organization

The test tree is organized by behavioral contract rather than by the implementation layout under `src/`.

Use this rule when adding or moving a test:

> Put a test where a contributor would look when the behavior breaks, not where the implementing function currently lives.

The top-level categories are:

- `fundamentals/`: Hilbert spaces, operator construction, indexing primitives, names, and canonical ordering semantics.
- `algebra/`: arithmetic, commutators, normal ordering, simplification, substitution, and unitary transformations.
- `symbolic/`: symbolic coefficients, phases, averages, indexed constraints, and symbolic reductions.
- `numeric/`: symbolic-to-numeric conversion, indexed materialization, backend contracts, and expectation-value conversion.
- `presentation/`: text/Unicode and LaTeX rendering behavior.
- `integration/`: tests that deliberately cross subsystem boundaries, such as canonical model workflows and algebraic invariants.
- `quality/`: static and repository-wide quality gates such as JET, Aqua, ExplicitImports, and CheckConcreteStructs.

Directories encode concepts, filenames encode contracts, and nested `@testset`s encode scenarios. Tests that intentionally span several concepts belong in `integration/` rather than being attached to an implementation module.
