# Exact closed-adjoint prototype evidence

This directory is the durable evidence archive for the agentic prototype workflow described in
[`../exact-closed-adjoint-plan.md`](../exact-closed-adjoint-plan.md). The reports distinguish
prototype evidence from production commitments. All Wave 1 and Wave 2 spikes were developed in
detached worktrees at commit `4a26c868a6a76238e46b24805e2f98fde69004d3` and reproduced by the
coordinator through the Julia MCP before their conclusions were accepted.

The concise decision matrices answer *what was decided*. The per-spike reports and retained Julia
harnesses answer *why*, *against which alternatives*, and *with what evidence*.

## Evidence index

| Wave | Question | Detailed report | Retained harness |
| --- | --- | --- | --- |
| 0 | What behavior and complete-workflow cost must later work preserve? | [`wave0-baseline.md`](wave0-baseline.md) | Baseline is the package test suite and benchmark expressions recorded in the report |
| 1A | Which closure coordinate and Hilbert-context boundary are correct? | [`wave1/closure-and-hilbert.md`](wave1/closure-and-hilbert.md) | [`closure_prototype.jl`](harnesses/wave1/closure_prototype.jl), [`run_closure_prototype.jl`](harnesses/wave1/run_closure_prototype.jl) |
| 1B | Which multimode ordering, ladder/quadrature convention, and timed-gauge decomposition are correct? | [`wave1/gaussian-and-gauge.md`](wave1/gaussian-and-gauge.md) | [`gaussian_conventions.jl`](harnesses/wave1/gaussian_conventions.jl), [`run_gaussian_conventions.jl`](harnesses/wave1/run_gaussian_conventions.jl) |
| 1C | Which exact bounded-response solver is viable over `CNum`? | [`wave1/exact-response-solver.md`](wave1/exact-response-solver.md) | [`exact_response_solver.jl`](harnesses/wave1/exact_response_solver.jl), [`run_exact_response_solver.jl`](harnesses/wave1/run_exact_response_solver.jl) |
| 2A | Should analyzed maps materialize directly or survive in a new wrapper/storage form? | [`wave2/materializer.md`](wave2/materializer.md) | [`materializer_prototype.jl`](harnesses/wave2/materializer_prototype.jl), [`run_materializer_prototype.jl`](harnesses/wave2/run_materializer_prototype.jl) |
| 2B | Can exact strategies use ordinary concrete dispatch without a function registry? | [`wave2/strategy-dispatch.md`](wave2/strategy-dispatch.md) | [`exact_strategy_interface.jl`](harnesses/wave2/exact_strategy_interface.jl), [`run_exact_strategy_interface.jl`](harnesses/wave2/run_exact_strategy_interface.jl) |
| 2C | Do flat coordinates and an explicit-Hilbert API survive product-space pressure? | [`wave2/product-space-api.md`](wave2/product-space-api.md) | [`product_api_spike.jl`](harnesses/wave2/product_api_spike.jl), [`run_product_api_spike.jl`](harnesses/wave2/run_product_api_spike.jl) |
| Production seam | Does the direct materializer remain acceptable when inserted into a real constructor? | [`production-materializer-experiment.md`](production-materializer-experiment.md) | Current uncommitted source and test diff; no separate copied harness |

The aggregate views remain useful:

- [`wave1-decision-matrix.md`](wave1-decision-matrix.md) maps the Wave 1 questions to the
  coordinator reproductions.
- [`wave1-foundation-report.md`](wave1-foundation-report.md) records the Wave 1 checkpoint.
- [`wave2-decision-matrix.md`](wave2-decision-matrix.md) maps the Wave 2 questions to the
  coordinator reproductions.
- [`wave2-architecture-report.md`](wave2-architecture-report.md) records the accepted architecture
  boundary and migration rule.

## Decision traceability

| Plan decision | Decisive evidence |
| --- | --- |
| Reuse `QTerm`; do not add `OperatorAlgebraBasis` or `BasisMonomial` identity | Wave 1A wrapper comparison, equality semantics, and product-valued closure fixtures |
| Treat explicit `HilbertSpace` as authoritative | Wave 1A renamed, duplicate-factor, indexed, ordinary N-level, and collective-transition fixtures |
| Persist site-interleaved bosonic coordinates | Wave 1B permutation equivalence and 33-mode materialization comparison |
| Share only the nonscalar timed-gauge reconstruction | Wave 1B body-velocity checks and the Fock-versus-quadrature scalar-lift counterexample |
| Split structural response blocks before exact solving | Wave 1C two-by-four versus whole-eight solve measurements and certificates |
| Retain scalar and explicit `2x2` solvers | Wave 1C expression-size and warm-cost comparison |
| Keep `UnitaryTransform` as the sole terminal representation | Wave 2A direct, transient-wrapper, and retained-nested storage comparison |
| Select exact strategies through concrete method dispatch | Wave 2B inference, JET, ambiguity, refusal, and branch-safety fixtures |
| Use one flat product-space coordinate vector | Wave 2C coupled four-coordinate Gaussian and cross-site Pauli fixtures |
| Require explicit `hilbert` in the first generic public slice | Wave 2C mismatch/refusal and keyword-adapter checks |
| Defer free indexed families to a binder-aware design | Wave 1A missing-scope evidence and Wave 2C wildcard-equality failure analysis |
| Migrate a named constructor only after an end-to-end gate | Wave 2A prototype comparison and the production materializer experiment |

## Reproduction contract

The retained harnesses are historical experiments, not package tests. Reproduce them against the
frozen commit, not by assuming that private APIs on a later branch are compatible.

1. Create a detached worktree at `4a26c868a6a76238e46b24805e2f98fde69004d3`.
2. Copy one report's two harness files into a `prototype/` directory in that worktree.
3. Instantiate the frozen root project in the Julia MCP.
4. Set Julia's working directory to the harness directory and include the corresponding
   `run_*.jl` file. Running from that directory is required for the Wave 1B and Wave 1C runners,
   which intentionally retain their original relative `include` calls.
5. For JET, use the version and target-module scope stated in the detailed report. Do not compare
   unrestricted dependency reports with package-scoped production JET gates.

The Wave 0 environment was Julia 1.12.6 with 12 threads. Its root and test project hashes are in
[`wave0-baseline.md`](wave0-baseline.md). The Wave 1A isolated JET environment used JET 0.11.4;
its manifest hash was
`7847149245e4b45fe0c7e132121c86f793b36f23b47807380b4a244ae761762f`. The Wave 2A isolated
environment used JET 0.12.1; its manifest hash was
`d52b6e223c7ebdb904a274e436753d065977df45471484d31a4ed24023814f96`.

## Evidence limits

- Prototype timings are warm directional measurements from the sample counts stated in each
  report. They are not CI thresholds.
- Agent conclusions entered the decision matrices only after the coordinator reproduced the
  smallest decisive fixture. Raw terminal transcripts were not retained.
- Generated manifests were not copied into this archive. Their relevant package versions and
  hashes are recorded above; the Wave 0 root and test manifests remain the authoritative package
  environment snapshots.
- Disposable code may use broad or allocation-heavy implementations to isolate a mathematical or
  architectural question. A successful prototype is not permission to copy it into production.
- The production materializer experiment occurred after the disposable Wave 2A spike and has a
  separate, explicitly incomplete verification status.

