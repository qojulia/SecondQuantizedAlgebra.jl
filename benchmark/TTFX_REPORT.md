# TTFX investigation: operator-product specialization

This note records a time-to-first-execution (TTFX) investigation of `SecondQuantizedAlgebra`, motivated by the first-call latency of `QuantumCumulants.meanfield` + `complete`. It documents where the compile-time cost goes, the change made in this branch, and the larger architectural option for future work.

## Methodology

TTFX is reported as **wall-clock time of the first call in a fresh `julia --startup-file=no` process, median of 5 runs**, after precompilation is already done. This isolates first-call latency (inference, native codegen, execution) from precompilation. The `@snoop_inference` instrumented time is used only for *ranking* MethodInstances and packages; its absolute numbers carry roughly ±1.5 s of noise on this machine and are not reliable for A/B magnitude. Measurements were taken on Julia 1.12.6 using the Jaynes-Cummings model (one Fock mode, one two-level atom) as the workload.

## Where the first-call time goes

Breakdown of first-call inference for `meanfield` + `complete` (JC model), grouped by the package owning each inferred method:

| Bucket | Share | What it is |
| :-- | --: | :-- |
| Base / stdlib | ~46% | `Dict` hashing/rehash, array growth, `Broadcast`, specialized on symbolic key/value types |
| SymbolicUtils | ~28% | v4 `Add`/`Mul` worker buffers specialized on many argument-tuple types |
| SecondQuantizedAlgebra | ~16-19% | operator algebra (`*` and the canonicalization passes) |
| Moshi | ~2.5% | the sum-type machinery the SymbolicUtils v4 backend uses |

There is no single hotspot. About three quarters of the cost is the symbolic coefficient/expression stack (Base plus SymbolicUtils), which is upstream of this package and independent of how operators are represented. The remaining ~16-19% is this package's own operator algebra, which is the part addressable here.

## Operator specialization inventory

There are 7 concrete `QSym` types: `Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`. The per-pair methods grow as O(n^2):

| Function | Methods | Growth |
| :-- | --: | :-- |
| `_site_compare` | 12 | per operator pair |
| `_can_commute` | 12 | per operator pair |
| `_commute_pair` | 5 | per operator pair |
| `_reduce_pair` | 3 | per operator pair |
| `adjoint`, `_ground_state_expand`, `_type_order` | ~9 | per operator type |

Operator products flow through `Vector{QSym}` (in `QTerm.ops` and the canonicalization pipeline) and `QTermDict = Dict{QTerm, CNum}`. Because the element type is the abstract `QSym`, every `Base` container method specializes on the abstract type and element access dispatches dynamically.

## Change in this branch

Two source-level changes, both verified to preserve results:

1. **Function-barrier return-type assertions** in the canonicalization passes (`src/algebra/passes.jl`, `src/algebra/pipelines.jl`, `src/expressions/index.jl`). The leaf functions `_site_compare`, `_can_commute`, `_commute_pair`, `_reduce_pair`, and the `.index` field access each have a concrete return type, but inference widens them to `Any` because they are called on elements of the abstract `Vector{QSym}`. Pinning the true return type at each call site (`::SiteCmp`, `::Bool`, `::Tuple{...}`, `::Index`) stops `Any` from cascading into the downstream code. This does not remove the dynamic dispatch into those functions (that is inherent to the heterogeneous product), but it removes the secondary dispatch cascade: `@report_opt` on a product goes from 15 to 5 runtime-dispatch sites.

2. **`@nospecialize` / `Base.@nospecializeinfer`** on the three operator-product methods `*(::QSym, ::QSym)`, `*(::QAdd, ::QSym)`, `*(::QSym, ::QAdd)`. These are single methods that Julia specialized once per concrete operator pair (`*(::Destroy, ::QAdd)`, `*(::Create, ::Create)`, and so on), even though their bodies only box the operand into a `Vector{QSym}` and run the pipeline, so the concrete operator type is irrelevant to the body. Collapsing those specializations cuts native codegen for the product entry points.

### Measured effect

Wall-clock first-call (`meanfield` + `complete`, JC model), median of 5 fresh runs:

| State | build | meanfield | complete | total first-exec |
| :-- | --: | --: | --: | --: |
| baseline | 2.49 s | 13.10 s | 1.01 s | **17.13 s** |
| this branch | 1.92 s | 12.29 s | 0.70 s | **15.40 s** |

A reproducible **~1.7 s (~10%)** reduction in first-call latency (run-to-run spread ~0.4 s). Operator-product results and runtime are unchanged (for example `a*a'` still gives `1 + a'*a`, product runtime ~1.5 us). Note that `@nospecialize` reduces codegen rather than inference, so this win shows in wall-clock first-call but not in `@snoop_inference` totals.

## Architectural option (future work)

The largest remaining lever on this package's share is to collapse the 7 `QSym` subtypes into one concrete sum type, for example via Moshi's `@data` (the same mechanism SymbolicUtils v4 uses for `BasicSymbolicImpl`). This would:

- make `Vector{QSym}` a concrete-element vector, so `Base` container methods specialize once and element access stops being dynamic dispatch; and
- collapse the ~32 per-pair O(n^2) methods into a handful of functions that branch on a tag (via `@match`), turning `_site_compare(::QSym, ...)` runtime dispatch into static calls.

This targets the SecondQuantizedAlgebra share fully plus a slice of the Base share. It does not touch the SymbolicUtils coefficient stack, which is the remaining floor. It is a substantial refactor (rewriting the 7 operator types and the per-pair methods as tagged variants and `@match` branches, with a closed-world variant set), and the per-pair methods encode real algebraic correctness (Pauli epsilon, Transition composition, Spin commutators), so each branch must be preserved and guarded by the test suite. A reference datapoint for the magnitude: a 100-variant stress test compiles in ~6.7 s as a sum type versus ~56 s as a `Union`.

## Summary

- The dominant first-call cost is the upstream symbolic stack (~74%), not the operator algebra.
- This branch removes the operator-product codegen explosion for a measured ~10% first-call reduction, with no change to results or runtime.
- A sum-type collapse of the operator hierarchy is the principled next step for this package's remaining share, scoped as its own project.
