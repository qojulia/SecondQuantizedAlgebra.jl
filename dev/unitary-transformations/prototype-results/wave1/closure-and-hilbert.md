# Wave 1A prototype report: closure coordinates and Hilbert context

## Outcome

The prototype supports the proposed closed-adjoint layer, with two changes to the plan:

1. Reuse `QTerm` as the private monomial coordinate instead of introducing
   `BasisMonomial`. `QTerm` already has the required concrete storage, inequality constraints,
   structural equality, and cached hash.
2. Make explicit `HilbertSpace` context the reliable generic path. Operator inspection can be
   a checked convenience path, but it cannot reconstruct the original space and must refuse
   incomplete metadata.

The successful examples all close quickly in a concrete `Matrix{Coeff}`. Kerr and mixed
spin-boson generators demonstrate deterministic bounded refusal rather than accidental runaway.
The production implementation should not copy this spike wholesale: the code intentionally
recomputes no columns, but it does use general `commutator` and `expand_completeness` pipelines
and allocates ordinary dense action matrices.

## Artifacts and execution

- `closure_prototype.jl`: disposable records, discovery loop, metadata spike, measurements.
- `run_closure_prototype.jl`: fixtures and independent algebraic checks.
- `jet_env/`: disposable environment used only for JET 0.11.4.

The runner was executed through the Julia MCP on Julia 1.12.6. All assertions passed. Timings
below are warm medians from 21 samples; they are directional prototype measurements, not
benchmark thresholds.

## Closure algorithm tested

The worklist begins with user-selected generators and their adjoints. Each coordinate is a
`QTerm`. For every unvisited coordinate it computes

```math
[G, X_j],
```

applies transition completeness, and immediately inserts unseen monomials. Identity terms are
stored in a separate affine offset. The final action is a concrete `Matrix{Coeff}`.

Three independent limits are checked in deterministic order:

- commutators computed;
- number of basis monomials;
- maximum operator-product degree.

The refusal result records the limit and first offending degree. Production diagnostics should
also print the generator, source coordinate, and unseen monomial; the prototype omits those to
keep its records small.

## Exact fixture results

| Fixture | Status | Coordinates | Commutators | Max degree | Warm time | Warm bytes |
|---|---:|---:|---:|---:|---:|---:|
| Number rotation | closed | 2 | 2 | 1 | 5.9 us | 12,304 |
| One-mode squeeze | closed | 2 | 2 | 1 | 7.1 us | 14,064 |
| Two-mode rotation | closed | 4 | 4 | 1 | 11.9 us | 26,944 |
| Two-mode squeeze | closed | 4 | 4 | 1 | 11.9 us | 26,944 |
| Pauli rotation | closed | 3 | 3 | 1 | 8.3 us | 24,576 |
| Spin rotation | closed | 3 | 3 | 1 | 8.3 us | 25,408 |
| Three-level mixing, ground 2 | closed | 8 | 8 | 1 | 16.6 us | 53,664 |
| Concrete indexed number rotation | closed | 2 | 2 | 1 | 5.9 us | 12,304 |
| Kerr on linear Fock seed | max degree | 8 at refusal | 7 | 9 | 65 us | 96,320 |
| Spin-boson coupling | max basis | 48 | 30 | 5 | 163 us | 331,384 |

The exact Fock action for `G = a'a` in the ordered basis `(a,a')` was

```text
[-1  0
  0  1].
```

For `G = (a'^2-a^2)/2` it was

```text
[ 0 -1
 -1  0].
```

For Pauli `G = sigma_z` in `(sigma_x,sigma_y,sigma_z)` it was

```text
[  0 -2im  0
  2im  0    0
   0   0    0].
```

For `NLevelSpace(:atom, 3, 2)`, `fundamental_operators` plus adjoints produced
the expected `n^2-1 = 8` independent coordinates, excluding the selected ground projector
`sigma_22`. Completeness projected every commutator back into those coordinates.

## Independent oracles

The runner checks known algebra laws directly, separately from discovery:

- `[a'a,a] = -a` and `[a'a,a'] = a'`;
- `[sigma_z,sigma_x] = 2im sigma_y`;
- `[S_z,S_x] = im S_y`;
- `[sigma_12 + sigma_21, sigma_12] = sigma_22 - sigma_11`, followed by the
  configured ground-projector completeness relation.

These pin sign, normalization, and transition-completeness conventions. Matrix exponentiation
was intentionally left to the exact-strategy spike; this report establishes only exact closure
and projection.

## `QTerm` versus a new `BasisMonomial`

The competing wrapper duplicated `QTerm`'s `Vector{Op}`, inequality vector, cached hash,
equality, and dictionary behavior. On the eight-coordinate N-level fixture, 200 dictionary
builds measured:

| Representation | Construction bytes | Lookup bytes |
|---|---:|---:|
| Reused `QTerm` | 134,400 | 0 |
| Copied `BasisMonomial` | 416,000 | 0 |

Lookup speed was indistinguishable at this scale. The wrapper allocated about 3.1 times as much
during construction because it copied the already-owned vectors. A wrapper also risks drifting
from canonical `QTerm` equality and inequality semantics.

Recommendation: reuse `QTerm`, but add one internal, trusted way to turn a coordinate back into
a one-term `QAdd` while preserving `ne`. Do not reconstruct it by multiplying `term.ops` as this
spike does; that loses scoped inequality information. Summation-scoped terms need a separate
indexed-family prototype before they are accepted.

## Minimal `AdjointBasis` layout

The smallest useful production record is:

```julia
struct AdjointBasis{H<:HilbertSpace}
    space::H
    terms::Vector{QTerm}
    lookup::Dict{QTerm,Int32}
end
```

`max_degree` is discovery telemetry, not permanent basis state. The action, affine offset,
status, limits, and diagnostics belong to a separate analysis result. A successful transform
should materialize into the existing public transformation type and should not retain the
discovery worklist or lookup table.

This layout is concrete for every concrete `H`. The prototype returned
`ClosureResult{NLevelSpace}`, `Matrix{Coeff}`, and `Vector{Coeff}` under `@inferred`.

## Hilbert-space inference boundary

### Metadata that survives in `Op`

- algebra family through `OpKind`;
- `space_index` within a product;
- operator display/family name;
- abstract or concrete index tag;
- N-level dimension and nondefault ground state for ordinary `Transition`;
- transition level pair or spin/Pauli axis.

This was sufficient for renamed Fock operators, duplicate Fock factors at distinct product
positions, a concrete indexed site, and `NLevelSpace(n=3, ground=2)`. A concrete
`Dict{Int32,InferredSlot}` worked without abstract values.

### Metadata that is absent

- the Hilbert-space name;
- untouched factors and therefore total `ProductSpace` shape;
- symbolic N-level labels;
- the index range and originating Hilbert space;
- the level count of `CollectiveTransition` (`nlev == 0`);
- whether inferred occupied slots with gaps correspond to omitted factors.

Operator names also matter to closure: differently named operators at the same space index
commute, so `fundamental_operators(h)` with default names is not a valid replacement for renamed
operators found in `G`.

### Recommended API policy

- Keep `hilbert = h` as the authoritative generic path.
- Inspect actual operators in `G` to preserve names and concrete indices even when `h` is given.
- Permit inference only when every occupied slot has complete algebra metadata and the requested
  transformation does not require unknown product factors.
- Refuse collective transitions without explicit `CollectiveNLevelSpace` because their `Op`
  omits `n`.
- Do not claim that inference recreates the user's original `HilbertSpace`; it produces only the
  algebra context needed for closure.

## Closure seeding and transition completeness

For transformation coverage, seed the complete independent local algebra for every touched
operator family, using names taken from `G`, then add adjoints. This gives:

- two ladder coordinates per touched Fock mode;
- both quadratures per phase-space mode;
- three Pauli or Spin coordinates;
- `n^2-1` matrix units for ordinary N-level spaces, excluding the configured ground projector;
- all `n^2` collective matrix units when explicit space context supplies `n`.

Starting from only `sigma_12` and its adjoint also closed, but completeness introduced four
monomials `(sigma_12,sigma_21,sigma_11,sigma_33)`. That smaller worklist does not provide rules
for every fundamental transition. Complete algebra seeding is therefore the correct default for
a reusable transform, while named constructors may use smaller certified bases when their rule
coverage is known in advance.

## Deterministic refusal semantics

Kerr acting on `(a,a')` produced successively higher odd-degree monomials and hit degree 9 after
seven commutators with `max_degree=7`. The spin-boson generator `(a+a')S_z` generated mixed
products and reached the 48-coordinate bound after 30 commutators.

Recommended default diagnostic fields:

```text
reason, configured limit, observed count or degree,
source basis coordinate, first unseen monomial
```

Limit order must be documented and tested so the same input reports the same reason. A refusal
is not evidence that the physical algebra never closes; it means exact finite closure was not
certified within declared work bounds.

## Type stability and JET

- Every successful and refused fixture passed `@inferred`.
- Result, basis, action, offset, lookup, and inferred-slot storage have concrete element types.
- JET 0.11.4 on Julia 1.12.6 reported zero correctness findings and zero optimization findings
  when restricted to `target_modules=(ClosurePrototype,)` for the Fock closure entry point.
- An unrestricted JET descent produced 603 reports in Base/SymbolicUtils internals. This matches
  the repository's existing practice: production gates must restrict `target_modules` to
  `SecondQuantizedAlgebra` rather than enlarge allowlists for upstream compiler noise.

The production version should add targeted JET fixtures for one Lie-algebra closure and one
Gaussian closure. The refusal path should remain in `@noinline` formatting helpers so it does not
pollute successful inference.

## Competing designs rejected or deferred

### New `BasisMonomial`

Rejected. It duplicates `QTerm`, increases construction allocation, and would need to replicate
index-constraint canonicality.

### `Vector{Op}` coordinates

Rejected. Kerr, mixed spin-boson, higher moment, and other valid finite closures require operator
products.

### Infer and recreate `HilbertSpace` from `Op`

Rejected. The required metadata is absent for collective transitions, symbolic level labels,
index ranges, untouched product factors, and original space names.

### Always seed only operators appearing in the generator

Rejected as the general default. It can certify a small invariant subspace but does not provide
complete transformation coverage. It remains appropriate for narrow named constructors.

### Indexed sums in the first generic slice

Deferred. `QTerm.ne` is adequate for monomial identity, but reconstructing summation scope and
projecting indexed families requires a separate design. Concrete indexed operators already work.

## Production recommendation

Implement the first closed-adjoint discovery slice with:

1. explicit Hilbert context plus names/indices collected from the actual generator;
2. complete independent local coordinates for touched spaces;
3. `QTerm` coordinates and a concrete `Dict{QTerm,Int32}` lookup;
4. separate identity offset;
5. fixed deterministic limits and compact refusal diagnostics;
6. the existing `commutator` and `expand_completeness` canonical pipelines;
7. no public analysis object retained in `UnitaryTransform`;
8. concrete indexed operators accepted, summation-scoped/index-family closure deferred.

The one-mode and named-constructor fast paths should continue to bypass discovery where their
closed actions are already known.

## Safe to discard

Everything under `prototype/` is disposable and no production source file was changed. The
`jet_env/Manifest.toml` and `jet_env/Project.toml` exist only to reproduce the JET check. The root
worktree's ignored `Manifest.toml` was also instantiated locally for the Julia MCP session and is
safe to discard with the isolated worktree.

