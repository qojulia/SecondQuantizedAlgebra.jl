# Wave 2C: flat product coordinates and public API pressure

## Question

Can the planned exact closed-adjoint engine use one flat `QTerm` coordinate space over an
explicit `HilbertSpace`, including cross-factor Gaussian and finite cross-site closures, without
introducing a hierarchical product-space basis? What should the first public `hilbert` keyword
boundary accept, infer, and refuse, and can its adapter remain concrete and type stable?

## Competing designs

1. **Flat exact coordinates (tested):** `Vector{QTerm}` plus `Dict{QTerm,Int32}`, with the supplied
   `HilbertSpace` as authoritative context. Factor support is computed from each `Op.space_index`.
2. **Hierarchical per-factor coordinates (rejected for the first boundary):** nested basis objects
   whose tensor structure duplicates information already present in `QTerm` and `ProductSpace`.
3. **Permanent support metadata (investigated, not recommended):** a parallel sorted
   `Vector{Vector{Int32}}` or `UInt64` mask per coordinate.
4. **Automatic exact Hilbert reconstruction (rejected):** reconstruct a complete `HilbertSpace`
   from operators before dispatch.
5. **Wildcard indexed `QTerm` equality (rejected):** weaken ordinary dictionary identity so free
   index labels compare alpha-equivalent.
6. **Explicit future family descriptor (viable later):** canonicalize declared binders into a
   separate exact key before ordinary dictionary lookup; do not alter `QTerm` equality.

The tested public shape is one thin method
`analyze_frame(generator::QAdd; hilbert=nothing)`. It validates the keyword and expression, then
immediately calls `_analyze_frame(generator, hilbert::H) where H<:HilbertSpace`.

## Fixtures and refusal cases

- Two renamed bosonic modes on a `ProductSpace` containing duplicate `FockSpace` factor types.
  The Hermitian generator contains hopping and anomalous terms:
  `a' b + b' a + a b + a' b'`.
- A finite cross-site Pauli closure under `G = sigma_z1 sigma_z2`, using the flat coordinates
  `{sigma_x1, sigma_y1 sigma_z2}`.
- An ordinary three-level transition with non-default ground level `2`, both unindexed and on a
  concrete indexed site.
- A renamed unindexed collective transition on a three-level collective space.
- Renamed duplicate factor types are distinguished by `space_index`, not operator name.
- Refusals cover an absent factor, wrong factor role/order, wrong transition dimension/ground
  metadata, a mismatched operator/index factor, a missing `hilbert` keyword, and a free abstract
  indexed family.
- A validation counter proves all Hilbert mismatches and missing-context cases fail before private
  analysis starts.

## Independent oracle

The Gaussian action was compared entry-by-entry with the hand-derived matrix for
`L(X) = im [G,X]` in `(a,b,a',b')` order:

```text
[  0  -im   0   im
  -im   0   im   0
    0  -im   0   im
  -im   0   im   0 ]
```

The cross-site oracle follows directly from Pauli multiplication and was compared with
`[0 2; -2 0]`. Hilbert validation was independently checked against the factor type and stored
ordinary-transition `(n_levels, ground_state)` metadata, rather than inferred operator names.

## Correctness result

- Final Julia-MCP run: **27/27 tests passed**.
- Both action matrices projected exactly into flat `QTerm` coordinates; no product hierarchy was
  needed.
- Gaussian support was `[[1], [2], [1], [2]]`; the cross-site support was `[[1], [1,2]]`.
- Renaming either the Hilbert factors or operators did not affect positional validation.
- Duplicate factor types were unambiguous when operators carried explicit positions.
- Concrete indexed sites validated. Free indexed families were deterministically refused.
- Raw `QTerm` keys for renamed index labels remained unequal, while the disposable explicit-binder
  signature mapped the two terms to the same alpha-equivalence signature. This establishes
  feasibility without wildcard dictionary equality; it is not an implementation of #234.

## Inference and JET result

- `@inferred` passed for Gaussian projection, cross-site projection, the Hilbert validator, and the
  keyword-to-positional adapter.
- `FlatCoordinates{H}` was concrete; its keys were `QTerm` and values `Int32`.
- Targeted `JET.report_call` reported no errors for projection or the adapter.
- Targeted `JET.report_opt` reported no optimization errors for either path.
- The private positional method precompiled successfully.
- `Test.detect_ambiguities(ProductAPIPrototype; recursive=false)` returned no ambiguities.
- Aqua cannot inspect a non-toplevel disposable module and reported that limitation. The package
  baseline `Aqua.test_ambiguities([SecondQuantizedAlgebra]; recursive=false)` and
  `Aqua.test_piracies(SecondQuantizedAlgebra)` both passed. The spike adds no methods to external
  generics whose argument types it does not own.
- The concrete public discovery shape was one `QAdd` method with keyword `hilbert`; its docstring
  was discoverable through `Base.Docs.doc`.

## Runtime, allocations, and expression size

Warm medians, 20 samples with one evaluation per sample:

| Complete spike operation | Time | Bytes | Allocations |
|---|---:|---:|---:|
| Project four-coordinate Gaussian action | 11.873 us | 22,256 | 230 |
| Project two-coordinate cross-site action | 3.562 us | 7,232 | 51 |
| Validate and dispatch public adapter | 1.188 us | 1,696 | 20 |

The four-coordinate basis occupied 1,568 bytes and its dense `4x4` coefficient action occupied
432 bytes. Parallel sorted-vector support metadata occupied 320 bytes. The same four supports as
`UInt64` masks occupied 72 bytes, a 4.4x reduction, but impose an unjustified 64-factor limit.
These microbenchmarks validate shape and bookkeeping only; they are not production performance
ceilings.

## Production concepts and lines introduced

- Production files or public API lines changed: **0**.
- Disposable prototype: 180 lines; disposable runner: 156 lines.
- The smallest production basis remains the Wave 1 shape:
  `hilbert::H`, `terms::Vector{QTerm}`, and `lookup::Dict{QTerm,Int32}`.
- Factor support should be analysis scratch derived from `term.ops`, not a fourth persistent basis
  field. A connected-block pass can use masks opportunistically for at most 64 factors and a
  general sorted-vector/bit-set fallback without changing coordinate identity.
- The public keyword adapter should contain validation/refusal only and dispatch immediately to a
  concrete private positional method.

## Failure modes

- Exact `HilbertSpace` reconstruction from operators is impossible in general. Operators do not
  retain Hilbert factor names; untouched factors are invisible; ordinary transitions omit symbolic
  level names; collective transitions do not carry the collective dimension. Therefore a partial
  scan can at most infer compatibility tags, not the authoritative Hilbert object.
- A supplied space with the same factor roles and transition metadata but different factor names is
  indistinguishable from the operator payload. The supplied `hilbert` must be documented as
  authoritative; validation can reject structural mismatches, not provenance mismatches.
- Free indexed families require binder scope, alpha-renaming, inequalities, coefficient indices,
  and substitution semantics. The small signature experiment covers only term-level feasibility
  and must not ship as a family key.
- `Matrix{Coeff}` does not support generic matrix multiplication because `zero(::Coeff)` is absent.
  The spike uses explicit structural comparisons. Production coefficient linear algebra must keep
  using package coefficient primitives/specialized loops rather than generic matrix algorithms.
- Persisting `Vector{Vector{Int32}}` support adds allocations; persisting `UInt64` support silently
  caps product factors. Neither belongs in the core coordinate record without a complete-workflow
  benchmark.
- Keyword arguments do not dispatch independently in Julia. Public overload separation must occur
  in positional argument types; the keyword adapter should validate `hilbert` and forward.

## Recommendation: adopt with a strict initial boundary

Adopt flat `QTerm` coordinates over the supplied `HilbertSpace`. Do not introduce a product-basis
tree or encode support in the persistent coordinate identity. Require `hilbert` for the first
generic exact closed-adjoint constructor; keep existing named constructors free to infer their
already-explicit operator sites. Validate factor count, operator role, ordinary-transition
metadata, and index/factor consistency before closure discovery, then call a concrete private
positional method.

Initially accept unindexed operators and concrete indexed sites. Refuse free indexed families with
a diagnostic pointing to the future indexed-family work. Preserve exact `QTerm` equality. In #234,
add an explicit binder-aware family descriptor and canonical family key if its full scope and
coefficient semantics pass their own spike.

This boundary supports cross-mode Gaussian and finite cross-site closed adjoints now, ordinary and
collective transitions through explicit Hilbert validation, and future family rules without
premature wildcard equality.

## Files safe to discard

- `prototype/product_api_spike.jl`
- `prototype/run_product_api_spike.jl`
- `prototype_report.md` after its evidence is copied into the coordinator decision matrix

All are confined to the isolated detached worktree. No production file was edited.

