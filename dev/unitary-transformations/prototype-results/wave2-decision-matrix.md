# Wave 2 architecture-pressure decision matrix

Status: complete. The coordinator reproduced all decisive harnesses through the Julia MCP on the
frozen Wave 0 commit.

| Decision | Required evidence | Decision |
| --- | --- | --- |
| Map-to-rules materializer | Exact equality with named constructors and complete-workflow time/allocation comparison | Adopt one private direct coefficient-row materializer; retain tiny named builders and route substantial matrix-rule loops only after end-to-end measurement |
| Hot unconditional storage | Conjugation, application, inversion, composition, display, and conversion with direct versus wrapped fields | Keep `UnitaryTransform` direct fields; transient wrapper added 7--10 allocations and retained nesting showed no repeatable benefit |
| Scalar-lift handoff | Complete timed gauge with concrete construction-time value and no persistent callbacks/provenance | Strategy passes positional `CNum`; construction consumes it into gauge `QAdd` and discards lift state |
| Exact-strategy dispatch | Concrete certificates or `nothing`, branch-safe fixtures, inference/JET, ambiguity and invalidation inspection | Ordinary concrete methods; supported blocks return inferred certificate, unsupported blocks return `nothing`; no registry or generic symbolic exponential |
| Certificate lifetime and fields | No duplicated terminal transform state; exact inverse and verification available before materialization | Ephemeral marker, forward, inverse, verification token, and lift token only; consume immediately |
| ProductSpace coordinate support | Cross-mode Gaussian and finite cross-site closed-adjoint projection | Flat `QTerm` coordinates work; derive factor support as analysis scratch rather than persistent basis state |
| Indexed-family compatibility | Concrete sites work and future alpha-equivalence does not require wildcard dictionary equality | Accept concrete indexed sites, refuse free families; #234 gets a separate binder-aware key without weakening `QTerm` equality |
| `hilbert` API boundary | Mismatch refused before analysis; inferred and explicit paths remain concrete and ambiguity-free | Require explicit `hilbert` for the first generic constructor; thin keyword validation forwards to a concrete positional method |

## Acceptance rule

The coordinator reruns each smallest decisive harness through the Julia MCP on the frozen Wave 0
commit. A design is accepted only when it preserves exact behavior, retains concrete storage, and
meets the plan's 10% runtime and 20% allocation replacement gates on representative complete
workflows. Prototype-only microbenchmarks may reject an abstraction but may not establish a CI
threshold.

## Coordinator reproduction summary

- Materializer: named and direct rules/applications/inverse laws matched for displacement,
  one- and two-mode Gaussian chains, and dense three-level rotation. Direct materialization stayed
  below the replacement gates; the wrapper added bytes and allocations. Numeric conversion agreed.
- Strategy interface: 47/47 tests passed, including all exact inverse identities and deterministic
  refusal of unsupported and root-unsafe blocks. Targeted `Matrix{CNum}` JET correctness and
  optimization checks reported zero findings; no ambiguities were detected.
- Product/API pressure: 27/27 tests passed. Exact flat-coordinate actions reproduced the analytic
  4x4 Gaussian and 2x2 cross-site Pauli matrices. Mismatches refused before analysis; targeted JET,
  package Aqua ambiguity/piracy, inference, precompile, and docs-discovery checks passed.
- Generic `LinearAlgebra` kernels remain outside the architecture because `Coeff` intentionally
  lacks `zero(::Type{Coeff})`. Shared matrix kernels use `_CNUM_ZERO` and explicit coefficient
  helpers.
