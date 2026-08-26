# Wave 2B prototype report: method-based exact strategies

## Question

Can exact closed-adjoint maps be selected with ordinary Julia dispatch, without a registry or
stored `Function`, while returning a concrete short-lived certificate that is sufficient for
rule materialization and does not become a second transform representation?

## Designs compared

1. **Direct method dispatch (implemented):** one `exact_certificate(::ConcreteBlock)` method per
   proven structure. A supported method returns a concrete `ExactCertificate`; a deterministic
   unsupported structure returns `nothing`; malformed claims throw `ArgumentError`.
2. **Function registry (discarded without implementation):** would add mutable global state,
   erase callable types, complicate invalidation, and provide no capability that methods lack.
3. **Generic symbolic matrix exponential (explicitly not implemented):** violates the refusal and
   branch policies and is unnecessary for the exercised exact structures.
4. **Persistent analysis/result object (discarded):** would duplicate the final rule maps and
   retain provenance after it stops being useful.

The prototype certificate is:

```julia
struct ExactCertificate{S,M,V,L}
    strategy::S
    forward::M
    inverse::M
    verification::V
    scalar_lift::L
end
```

`S` is a zero-size strategy marker. `V` is an exact-verification token. `L` is either
`NoScalarLift` or `ScalarLift{CNum}`. Production matrices are `Matrix{CNum}`. The closed-adjoint
basis, diagnostics, action matrix, and coordinate ordering remain in the analysis problem and are
not copied into the certificate.

## Fixtures and oracle

The runner exercises:

- two exact diagonal unit phases and supplied reciprocal phases;
- an exact plane rotation using the rational Pythagorean pair `(3/5, 4/5)`;
- an exact hyperbolic squeeze using `(5/4, 3/4)`;
- normalized involutions with both `A^2 = -I` and `A^2 = I`;
- a third-order nilpotent Jordan block, exponentiated by its finite polynomial;
- a user-supplied exact orthogonal map and exact inverse;
- an unsupported dense block;
- a scaled involution whose normalization would require choosing `sqrt(4)`;
- a four-coordinate Gaussian map converted from Wave 1's canonical site-interleaved order
  `(a1,a1',a2,a2')` to a temporary split-Nambu view and back;
- nonscalar body velocity `M^-1*dM` and static versus timed scalar-lift payloads.

The independent oracle is exact forward/inverse multiplication in both orders. Named strategies
also verify their defining identities before constructing a certificate: circular or hyperbolic
normalization, the declared involution square, or nilpotency at the declared order. The user-map
fixture supplies an exact orthogonal map; production must additionally invoke the Hilbert-space
algebra/canonical-map validator before certification.

## Correctness and refusal results

- `47/47` tests pass.
- Every accepted certificate passes `forward*inverse == inverse*forward == I` exactly.
- Every accepted certificate and every one of its field types is concrete.
- Unsupported dense and branch-unsafe scaled-involution blocks return `nothing`.
- Incorrect rotation/squeeze normalization, false involution claims, false nilpotency claims, and
  incorrect user inverses fail with deterministic `ArgumentError`s.
- Site-interleaved -> split-Nambu -> site-interleaved round trip is exact. This supports keeping
  the Wave 1 site-interleaved basis authoritative and using split order only inside a strategy.

## Branch and root boundary

Adopt these branch-safe initial strategies:

- diagonal canonical phase atoms with an exact reciprocal;
- named rotation and squeezing formulae whose parameter relations are already certified;
- normalized involutions proven to square to `+I` or `-I`;
- nilpotent actions, using the finite exponential polynomial;
- user-supplied exact maps/bases with supplied exact inverse and domain validation.

Do not synthesize a normalization with `sqrt(lambda)` for a general claim `A^2 = lambda I`.
Accept it only if a named strategy supplies a normalized action and verifies that normalization, or
if the caller supplies the complete exact map/inverse. Likewise, do not add a generic symbolic
spectral exponential. Root-producing dressed strategies remain behind the branch-aware work
already deferred in the plan.

## Gauge and scalar lift

`M^-1*dM` (and, for affine maps, `M^-1*dd`) determines the shared nonscalar body-gauge data. It
does not determine the c-number phase lift. A timed strategy must therefore construct a concrete
`ScalarLift{CNum}` under its named convention. An opaque timed user map must supply and validate
that lift; absence is a refusal, not an inferred zero. The lift is consumed into `gauge::QAdd` and
does not survive in `UnitaryTransform`.

The prototype intentionally models static input with `NoScalarLift`; the production call path,
which already knows `StaticTime` versus `DynamicTime`, should enforce the corresponding lift type
before materialization.

## Inference, JET, ambiguities, and invalidation surface

- `@inferred` passes for all seven accepted concrete families and the unsupported `nothing` path.
- `Test.detect_ambiguities(...; recursive=true)` reports zero ambiguities.
- There are eight `exact_certificate` methods and no registry.
- Targeted JET correctness and optimization reports are empty for rotation, nilpotent, user-map,
  and unsupported rational fixtures.
- A realistic `Matrix{CNum}` rotation is also `@inferred`, fully concrete, and has zero targeted
  JET correctness and optimization reports when analysis is scoped to the prototype module, as
  the repository's JET tests scope analysis to package code.

Adding another strategy adds one ordinary method specialized on a concrete analysis record.
Concrete public constructor call sites infer that method directly. Only callers that erase their
argument to an abstract block type expose an open-world return union and the usual method-table
backedge; such an abstract entry must not be the production hot path. This is a substantially
smaller invalidation surface than a mutable registry plus dynamically typed callables.

## Coefficient-kernel finding

The first generic prototype failed on the actual coefficient domain because `Coeff` intentionally
does not define `zero(::Type{Coeff})`; generic `zeros`, matrix multiplication, and several
`LinearAlgebra` paths therefore are not valid production kernels. The repaired prototype uses
private scalar operations backed by `_CNUM_ZERO`, `_CNUM_ONE`, `_add_cnum`, `_mul_cnum`, and
`_neg_cnum`, plus explicit dense matrix loops.

This is a foundation constraint, not merely an optimization: exact-strategy code must not call
generic `zeros(CNum, ...)`, `*`, `exp`, `inv`, or `\`. It should share the same private CNum matrix
kernels as the adjugate solver and materializer.

## Performance, allocation, and size

Exploratory warm medians use 31 samples and the complete certificate call, including exact
verification. They are suitable for architecture comparison, not CI ceilings.

| Fixture | Median time | Median bytes |
|---|---:|---:|
| diagonal 2-coordinate rational | 2.4 us | 624 B |
| rotation 2-coordinate rational | 2.4 us | 720 B |
| squeeze 2-coordinate rational | 2.4 us | 720 B |
| involution 2-coordinate rational | 3.6 us | 1,152 B |
| nilpotent 2-coordinate rational | 3.6 us | 1,824 B |
| user map 2-coordinate rational | 2.4 us | 432 B |
| unsupported dispatch | 1.2 us | 0 B |
| rotation 2-coordinate `CNum` | 0.11 ms | 20,928 B |

Repeated `CNum` medians varied between roughly 0.06 and 0.11 ms in this shared prototype session,
but allocation was stable at about 21 KB. This remains small relative to Wave 0's complete
millisecond-scale construction workflows; the coordinator should measure the integrated vertical
slice before setting a gate.

`Base.summarysize` for accepted rational certificates is 248--408 B. Those figures include both
dense maps. The strategy and static-lift markers are zero-size fields; the verification byte is
effectively padding-level storage.

## Duplication against `UnitaryTransform`

- `forward` and `inverse` logically duplicate the final forward/inverse rule dictionaries, but
  both are required temporarily: generic symbolic inversion is forbidden and the materializer
  must build both rule sets once. They must be dropped immediately after materialization.
- `scalar_lift` becomes the c-number part of `gauge` and must then be dropped.
- `strategy` is useful only for typed dispatch and diagnostics and must not persist.
- `verification` is a construction token, not public provenance, and must not persist.
- basis, coordinate order, action, offset, derivatives, and diagnostic traces do not belong in
  this certificate at all.

Thus the certificate is not a sibling transform type; it is a typed hand-off between exact solving
and one-time rule materialization.

## Concepts and lines

Disposable code is 279 lines; the runner is 201 lines. Much of that is fixture scaffolding:

- 2 abstract dispatch categories;
- 6 zero-size strategy markers;
- 1 certificate, 1 verification token, and 2 lift variants;
- 8 prototype block/input records;
- 1 small private exact matrix kernel;
- 8 strategy methods;
- ordering/body-velocity probes.

Production should not add the eight prototype block records as a parallel public model. Methods
should dispatch on the actual concrete closed-adjoint block/problem records, leaving roughly the
certificate/lift tokens, coefficient kernels, and selected strategy bodies as the reusable core.

## Failures and limitations

1. Generic matrix allocation/multiplication failed for `CNum`; repaired with private CNum kernels.
2. Unscoped JET descends through Symbolics internals and produces irrelevant dependency reports;
   repository-style `target_modules` scoping gives zero reports for the prototype code.
3. The prototype verifies inverse identities and named structural claims, but the final
   Hilbert-space algebra/canonical validator belongs at the certificate boundary and was not
   duplicated here.
4. Timed derivative construction is represented by the body-velocity and scalar-lift contract;
   a complete named gauge is owned by the Wave 1 Gaussian prototype/materializer, not rebuilt in
   this spike.
5. Measurements were made in a shared agent session and are not stable benchmark ceilings.

## Recommendation

Adopt the method-based design with a strict boundary:

1. Dispatch directly on concrete analyzed block types; no registry and no stored callable.
2. Return a concrete short-lived certificate or `nothing`; format refusal only after the caller
   knows no exact strategy accepted the block.
3. Store only strategy marker, forward map, exact inverse, verification token, and static/timed
   lift token in the certificate.
4. Validate the exact inverse, the strategy's defining identity, and the Hilbert-space algebra
   before materializing rules.
5. Use private explicit `CNum` matrix kernels, not generic `LinearAlgebra` arithmetic.
6. Keep site-interleaved order canonical; convert temporary split-Nambu strategy views back before
   certification/materialization.
7. Support branch-safe diagonal, named rotation/squeeze, normalized involution, nilpotent, and
   verified user-map strategies first. Refuse root-producing or generic spectral blocks.
8. Consume the certificate immediately into the existing `UnitaryTransform` and retain none of
   its analysis provenance.

## Discard list

- mutable strategy registries;
- `Function` fields, closures, or `Any` payloads;
- a generic symbolic matrix exponential;
- automatic roots for `A^2 = lambda I`;
- generic `inv`, `exp`, `zeros(CNum, ...)`, or matrix `*` on the coefficient hot path;
- a persistent `ExactAdjointMap` sibling to `UnitaryTransform`;
- persisted basis/order/action/derivative/provenance data in the certificate;
- inference through an abstract erased-block hot path;
- reconstruction of a timed scalar gauge from `(M,d,dM,dd)` alone.

## Reproduction

Julia was run only through the Julia MCP in the isolated detached worktree at `4a26c868`:

```julia
include("prototype/run_exact_strategy_interface.jl")
```

No production source or test file was edited. The isolated environment generated an ignored
`Manifest.toml` during instantiation; disposable tracked output consists only of `prototype/` and
this report.

