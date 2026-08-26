# Production materializer experiment

Status: **experimental, uncommitted, and not a deliverable**.

Captured: 2026-08-26, Europe/Zurich.

This experiment followed Wave 2A to test one terminal seam in real package code. It does not
implement the closed-adjoint engine, a public generic constructor, or multimode Gaussian-affine
Hamiltonian analysis. The changes are intentionally left in the worktree while the architecture
is discussed.

## Question

Can a private coefficient-row materializer replace a substantial existing matrix-to-rule loop
while preserving the current `UnitaryTransform`, exact public behavior, inference, and acceptable
warm construction cost?

This is an end-to-end pressure test of the Wave 2A recommendation. It is not evidence that a
standalone materializer should be merged.

## Worktree state

- Branch: `displace`
- Starting HEAD: `d23ecb4ced6d8bc7297714a0bd6b809686bb3a09`
- Files modified by the experiment:
  - `src/algebra/unitary.jl`
  - `src/algebra/unitary_constructors.jl`
  - `test/arithmetics/unitary_test.jl`
  - `test/quality/JET.jl`
- No commit was created.
- The user's other work remains untouched.

## Change under test

The experiment adds `_materialize_rule_map`, a private terminal bridge accepting:

```julia
generators::Vector{Op}
coordinates::Vector{QTerm}
coefficients::Matrix{Coeff}
offsets::Union{Nothing,Vector{Coeff}}
```

It immediately constructs the existing `Dict{Op,QAdd}` representation. Coordinate matrices,
analysis state, closures, callbacks, strategy metadata, and provenance do not survive in
`UnitaryTransform`.

The existing dense N-level `_matrix_unit_rules` loop was selected as the real consumer because it
is a substantial matrix-to-rule contraction. It now fills a short-lived concrete coefficient
matrix and calls the materializer. No public signature or terminal transform field changed.

## Behavioral oracle

The focused unitary test gained a dense rational orthogonal three-level matrix. For every one of
the nine matrix units, the public symbolic result

```julia
conjugate(Transition(atom, :tau, i, j), Rotation(transition, W))
```

is compared with the independent dense-matrix oracle `W' * E_ij * W`. The test also checks a
forward/inverse round trip through public transformation operations.

The final focused run completed with **265/265 tests passing**. This number belongs to the current
branch state and should not be compared directly with the frozen Wave 0 count of 277, whose test
contents differed.

## Inference and storage checks

- `@inferred Rotation(dense_transition, dense_basis)` returns `UnitaryTransform`.
- The materializer returns the existing concrete rule dictionary; no new public or retained
  storage type was introduced.
- The dense constructor was added to the package JET call-site list.

## Performance evidence

Two measurement forms were used and must not be conflated.

The package's `_warm_allocated` fixture measured **5,444 allocations** for the experimental dense
three-level constructor in the root environment. The test ceiling is 6,129, derived from the
previous recorded median of 5,108 plus 20% headroom.

A same-session BenchmarkTools comparison against a scratch reconstruction of the old direct loop
produced:

| Implementation | Time | Memory | Allocations |
| --- | ---: | ---: | ---: |
| Old direct loop | 983,093 ns | 215,696 bytes | 5,525 |
| New coefficient-row materializer | 853,677 ns | 239,440 bytes | 5,861 |

The rule dictionaries were exactly equal. In this run the shared path was about 13% faster, used
about 11% more memory, and performed about 6% more allocations. These are single-session prototype
comparisons, not CI performance claims. They support continued use of the seam for a real vertical
slice, but do not justify merging it by itself.

## JET evidence and unresolved check

A temporary quality environment was needed because attaching the nested test environment through
the Julia MCP resolved it against the root project unexpectedly.

With JET 0.12:

- package correctness analysis passed;
- all 41 targeted call-site checks passed, including the new dense N-level constructor;
- optimization analysis reported two findings in existing bounded-response kernels. They did not
  point to the new materializer, but the version-sensitive attribution was not fully closed.

The temporary environment was then constrained to JET 0.11.6 to match the repository's intended
quality line. That rerun was interrupted before a final result was captured. Therefore the
production JET gate remains **open**. It must not be summarized as passing.

## Interpretation

The experiment validates a plausible terminal materialization seam:

- it preserves public behavior and the existing terminal type;
- it provides a concrete consumer for coefficient-row materialization;
- it remains within the provisional allocation gate;
- it does not demonstrate a complete public feature;
- it does not yet have a completed repository-version JET result;
- it does not establish that tiny named constructors should migrate.

The materializer should land only in the same pull request as the first complete closed-adjoint or
Gaussian-affine vertical slice that needs it. If that slice fails representative end-to-end gates,
retain the specialized constructor and discard this experiment.

## Required follow-up before production

1. Select the first public vertical slice and make it the actual consumer.
2. Repeat old/new complete-workflow benchmarks in one controlled session.
3. Complete package-scoped JET correctness and optimization checks using the repository-supported
   JET version.
4. Run Aqua ambiguity and piracy checks.
5. Run the focused and full test-project suites.
6. Run fresh precompile/load and documentation checks.
7. Run `git diff --check`.
8. Review the resulting full diff and either retain or remove the experiment as one coherent
   implementation decision.

