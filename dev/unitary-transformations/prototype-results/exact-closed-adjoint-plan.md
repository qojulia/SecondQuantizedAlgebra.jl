# Exact closed-adjoint and Gaussian-affine transformation plan

Status: design draft for iteration. This file records architectural decisions and open
questions; it is not a public API commitment.

Related issues:

- [#148: Bogoliubov transformation](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/148)
- [#231: automatic exact rotating frames and dressed bases](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/231)
- [#234: generalized canonical maps](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/234)

## Decision

Build one exact closed-adjoint transformation architecture for every operator algebra supported
by the package. Treat Gaussian-affine transformations as the specialization for bosonic
Heisenberg algebras, not as the root abstraction.

The architecture should be complete in representation and deliberately limited in exact solution
strategies. Unsupported symbolic exponentials, spectra, or non-closing generators must be refused
instead of converted to floating point.

Keep the existing `UnitaryTransform` as the public result and rule-application type. New analysis
records and matrix representations remain private until a concrete public need justifies exposing
them.

## Goals

- Support exact transformations for Fock, phase-space, Pauli, Spin, transition, and collective
  transition operator bases through one internal architecture.
- Share basis discovery, adjoint-closure extraction, exact map materialization, inversion, gauge
  handling, and validation.
- Support general multimode Gaussian-affine maps, including coupled quadratic Hamiltonians,
  anomalous terms, and linear drives.
- Retain inexpensive closed-form paths for named one- and two-mode transformations.
- Make exactness, conditional validity, resonance surfaces, and refusal behavior explicit.
- Keep all persistent storage concrete and type-stable.
- Reduce maintained concepts and duplicate formulas rather than optimizing only for source-line
  count.

## Non-goals

- A universal symbolic matrix exponential or eigensolver.
- Time-ordered exponentials for noncommuting time-dependent generators.
- Treating approximately canonical or approximately closing maps as exact.
- Silently expanding unbounded or polynomially growing adjoint bases.
- Making every operator-valued generator supported merely because its leaf operator kind is known.
- Folding indexed-family addressing, collective algebra laws, and Gaussian algebra laws into one
  undifferentiated implementation.

## Hilbert-space algebra interface and adjoint bases

Do not introduce a second public or private type hierarchy parallel to `HilbertSpace`. The
existing concrete spaces already identify the supported algebras:

```text
HilbertSpace
├── FockSpace
├── PhaseSpace
├── PauliSpace
├── SpinSpace
├── NLevelSpace
├── CollectiveNLevelSpace
└── ProductSpace
```

Algebra-specific hooks should dispatch on these types. `fundamental_operators(h)` already supplies
the minimal generator seed for every space and should be reused rather than recreated inside the
transformation subsystem.

`FockSpace` and `PhaseSpace` remain two coordinate realizations of a bosonic Heisenberg algebra,
but that relationship is expressed through shared methods rather than another set of Julia types.

A particular transformation still needs a deterministic, independent coordinate set. Represent
that with a short-lived private `AdjointBasis{H}` value parameterized by the existing Hilbert-space
type. This calculation record owns:

- externally addressable rule generators;
- an independent closure basis;
- adjoint pairing;
- commutator structure;
- identity and completeness relations;
- exact map-validation laws;
- projection of an expression onto the closure basis.

Rule coverage and closure coordinates must be separate concepts. For example, all transition
operators need rules, while diagonal completeness can make a naive transition-plus-identity
coordinate set linearly dependent.

An expression does not currently retain its originating `HilbertSpace`: each `Op` stores a role,
site index, name, and packed role-specific metadata. Most local space information can be inferred
from those fields, but a collective transition does not carry the level count needed to reconstruct
its complete algebra. Generic APIs should therefore infer Hilbert-space context when the operator
metadata is sufficient and accept explicit `HilbertSpace` context for the general or ambiguous
case. Do not enlarge every `Op` or `QAdd` merely to persist construction provenance.

## Shared closed-adjoint model

For an independent basis `X = (X₁, …, Xₙ)`, extract the action of a generator `G` as

```math
i[G,X_j] = \sum_k A_{kj}X_k + b_j 1.
```

The affine part is homogenized internally by augmenting the action:

```math
\widehat A =
\begin{pmatrix}
A & b\\
0 & 0
\end{pmatrix}.
```

An exact strategy produces

```math
U^\dagger X U = M X + d.
```

A possible private representation is shown below. Wave 1 narrowed `AdjointBasis`: closure
coordinates reuse `QTerm`, while rule coverage, action matrices, affine offsets, and diagnostics
remain separate short-lived records rather than duplicated basis state.

```julia
struct AdjointBasis{H<:HilbertSpace}
    space::H
    terms::Vector{QTerm}
    lookup::Dict{QTerm,Int32}
end

struct ClosedAdjointProblem{H<:HilbertSpace}
    basis::AdjointBasis{H}
    action::Matrix{CNum}
    offset::Vector{CNum}
end

struct ExactAdjointMap{H<:HilbertSpace,T}
    basis::AdjointBasis{H}
    forward::Matrix{CNum}
    inverse::Matrix{CNum}
    shift::Vector{CNum}
    inverse_shift::Vector{CNum}
    gauge::QAdd
    time::T
    relations::Vector{ParamRelation}
end
```

These declarations are illustrative. Before adding them, compare the required fields against
`UnitaryTransform` and avoid storing the same information twice. If an analysis result is consumed
immediately, prefer a short-lived concrete record and materialize the rule dictionaries once.

### Basis extraction

1. Use an explicit `HilbertSpace` when supplied; otherwise infer the local space only when the
   operator metadata is sufficient.
2. Discover every touched concrete site and choose a deterministic basis ordering.
3. Seed the basis with `fundamental_operators(h)`, mapped to the observed operator names and sites.
4. Obtain an independent closure coordinate set using methods on the concrete `HilbertSpace`.
5. Compute commutators or use a specialized single-pass scanner.
6. Project every result exactly onto the basis plus its allowed central component.
7. Reject the first unsupported term with the generator and commutator in the error.
8. Decompose structurally disconnected blocks before invoking an exact solver.

For bosonic quadratic Hamiltonians, scan the Hamiltonian once rather than constructing one full
symbolic commutator per generator. Given canonical brackets `Cᵢⱼ = [zᵢ,zⱼ]`, use

```math
[z_i z_j,z_k] = z_i C_{jk} + C_{ik} z_j
```

to accumulate the affine action directly.

### Exact map materialization

The common materializer must:

- build complete forward and inverse `QAdd` rules;
- verify exact forward/inverse round trips;
- validate the map using methods on its concrete `HilbertSpace`;
- attach parameter relations or explicit conditional obligations;
- carry a complete gauge for timed transformations;
- return the existing concrete `UnitaryTransform{StaticTime}` or
  `UnitaryTransform{DynamicTime}`.

## Algebra-specific invariants

The exact validation law depends on the algebra:

| Hilbert spaces | Required invariant |
| --- | --- |
| `FockSpace`, `PhaseSpace` | canonical commutator matrix, adjoint/reality pairing |
| `PauliSpace`, `SpinSpace` | `su(2)` bracket, adjoints, and representation-specific identities |
| `NLevelSpace` | matrix-unit product, adjoints, and completeness |
| `CollectiveNLevelSpace` | collective Lie bracket and adjoints without local matrix-unit products |

A general map can carry an augmented central coordinate, but the concrete `HilbertSpace` methods
decide whether a nonzero affine shift is legal. Genuine scalar translations belong to Fock and
phase-space Heisenberg algebras. Pauli, Spin, and transition transformations are normally
homogeneous in their complete operator bases.

Canonicality proves preservation of the algebra; it does not prove that a map equals a requested
exponential. Construction certificates and canonicality diagnostics must remain distinct.

## Gaussian-affine specialization

For a bosonic canonical vector `z`, write

```math
[z_i,z_j] = C_{ij}, \qquad U^\dagger z U = Mz+d.
```

The homogeneous map is canonical when

```math
M C M^T = C
```

and it satisfies the appropriate adjoint or real-coordinate condition. Its inverse is derived
structurally:

```math
z \mapsto M^{-1}z-M^{-1}d.
```

For a quadratic-affine Hamiltonian,

```math
H = \frac12 z^T K z + f^T z + c,
```

the adjoint action closes on `z` plus the identity. The same extracted problem supports several
different operations; these operations must not be conflated into one ambiguous constructor.

### Static centering

`Displace(modes, Hlin)` computes the algebraic equilibrium shift that removes the linear terms.
The general solver should support:

- coupled modes and hopping;
- anomalous quadratic terms;
- arbitrary exact one-mode and multimode quadratic blocks;
- drives on any selected canonical coordinate.

Structurally singular driven systems are refused. Symbolic quotients are valid away from the
reported denominator surfaces. A singular system with zero forcing needs a deterministic identity
result rather than an arbitrary null-space displacement.

### Bounded harmonic response

`Displace(modes, Hlin, t)` scans a finite sum of constants and exact unit phases and solves each
harmonic response with zero homogeneous contribution. The harmonic scanner is shared across all
bosonic coordinate systems and immediately accumulates each solved component.

Arbitrary envelopes, nonlinear phases, transients, initial conditions, and numerical solutions
remain the responsibility of explicit-field displacement constructors.

### Exact Gaussian evolution

`UnitaryTransform(G, θ)` and `RotatingFrame(H0, t)` may use the homogenized affine action when an
exact matrix-exponential strategy is available. Initial strategies should include:

- diagonal phase blocks;
- exact one- and two-mode rotation blocks;
- exact one- and two-mode squeezing blocks;
- nilpotent blocks;
- involutive blocks with a proven scalar square;
- registered structured exact spectral decompositions.

An unsupported block is refused rather than passed to a generic numerical exponential.

### Gaussian diagonalization

`DressedFrame(H0)` is a separate objective from centering or evolution. Exact symplectic or
Bogoliubov diagonalization must define:

- stability and Hermiticity requirements;
- deterministic mode ordering and phase conventions;
- behavior under degeneracy;
- exact square-root and branch requirements;
- verification that the returned map is canonical and diagonalizes `H0`.

Only registered exact structures are supported. A general symbolic Williamson or Bogoliubov
eigensolver is not required.

## Gauge and phase convention

A timed affine operator map `(M(t), d(t))` does not determine the scalar part of
`im*(∂ₜU')*U`. Different unitary phase choices have the same operator action and different scalar
gauges.

Every timed solution strategy must therefore return or derive:

- the quadratic gauge coefficients;
- the linear gauge coefficients;
- the c-number gauge under a documented phase convention.

Examples:

- `RotatingFrame(H0, t)` has the certified gauge `-H0` under its convention;
- Weyl displacement uses the complete half-symplectic scalar term;
- named rotations and squeezes define a particular continuous lift;
- an opaque moving matrix needs a caller-supplied gauge plus compatibility validation.

The common map materializer must never invent the scalar gauge from `M` and `d` alone.
Wave 1 established that the shared nonscalar reconstruction uses body velocities
`M^-1*dot(M)` and `M^-1*dot(d)`. Each exact strategy supplies its scalar lift as a concrete
coefficient during construction; no callback or lift metadata survives in `UnitaryTransform`.

## Public API direction

Retain semantic constructors rather than exposing a single ambiguous `GaussianFrame(H)`:

```julia
# Explicit named maps
Displace(...)
Rotation(...)
Squeeze(...)
Bogoliubov(...)

# Hamiltonian-derived operations
Displace(modes, Hlin)
Displace(modes, Hlin, t)
UnitaryTransform(G, θ)
UnitaryTransform(G, θ; hilbert = h)  # explicit context when inference is insufficient
RotatingFrame(H0, t)
RotatingFrame(H0, t; hilbert = h)
DressedFrame(H0)
DressedFrame(H0; hilbert = h)
```

Normally the Hilbert-space context and adjoint basis are inferred from the concrete operators
touched by the arguments. Explicit `HilbertSpace` context is required where the operator
representation is ambiguous or incomplete, especially for collective transitions whose level
count is not encoded by one operator. The exact positional or keyword spelling of the explicit
context uses a keyword so the ordinary inferred spelling remains unchanged. The public keyword
method must immediately dispatch to positional, concretely typed private methods so inference and
JET can see the supplied Hilbert-space type.

Raw user-supplied Bogoliubov coefficients that are not provably canonical should return a future
conditional transformation carrying visible constraints, as planned in #234. Parametrized
`Rotation` and `Squeeze` remain unconditional because their construction proves the constraints.

## Exactness and refusal policy

- Require exact structural Hermiticity and closure.
- Manufacture no floating-point values from exact input.
- Use structural zeros; never delete small coefficients by tolerance.
- Reject exact resonances and singular driven blocks with component-level diagnostics.
- Preserve symbolic quotients and state their nonzero-domain requirements.
- Reject generators whose adjoint basis grows beyond the supported finite domain.
- Reject exact blocks for which no registered exponential or spectral strategy exists.
- Keep approximate numerical transformations in the distinct API tracked by #232.

## Performance and type stability

- Store only concrete `Op`, `CNum`, `QAdd`, matrix/vector, relation, and marker types.
- Store no `Any`, abstract container elements, `Function` fields, closures, or symbolic callbacks.
- Do not encode arbitrary mode count in the type solely to claim static sizing. `Matrix{CNum}` and
  `Vector{CNum}` are concrete; runtime dimensions avoid compilation blow-up.
- Retain specialized scalar and explicit `2×2` methods for common one-mode blocks.
- Add structured small-block methods only when representative benchmarks justify them.
- Scan each input Hamiltonian once and decompose its structural coupling graph before solving.
- Materialize `QAdd` rules only once after the coefficient-level map is certified.
- Keep refusal-message formatting behind `@noinline` boundaries.
- Use internal `CNum` arithmetic helpers on coefficient hot paths.
- Do not use symbolic-pivot Bareiss elimination with the current `CNum` domain: composite exact
  polynomial divisions are retained as unevaluated quotients. Prefer specialized small-block
  formulas and prototype a symbolic-pivot-free adjugate solve for larger blocks.

Inference and quality requirements:

- `@inferred` public constructor return types;
- JET correctness checks for all new entry points;
- optimization checks on dispatch-free extraction and solver kernels;
- Aqua ambiguity and piracy checks;
- no allowlist expansion merely to hide new runtime dispatch.

Benchmarks should exercise complete millisecond-scale workflows rather than isolated noisy
micro-operations:

- extraction, solve, materialization, and transformation;
- one-mode and coupled multimode Gaussian problems;
- one and 33 drive harmonics;
- representative Spin and transition closed-adjoint flows;
- inversion, composition, and timed gauges.

## Implementation phases

Each phase should leave the tree tested and reviewable. The phases share one architecture but do
not need to land in one pull request.

### Phase 0: lock down existing behavior

- Treat current named transforms and one-mode Hamiltonian-derived displacement tests as behavioral
  oracles.
- Record current warm allocations and end-to-end benchmark baselines.
- Add no abstraction until it replaces duplicated behavior or enables the next phase.

Exit criterion: current results, gauges, inference, and allocation gates are reproducible.

### Phase 1: Hilbert-space algebra hooks and adjoint bases

- Extend the existing `HilbertSpace` interface with private closure, projection, and map-validation
  hooks.
- Introduce the short-lived private `AdjointBasis{H}` calculation record and deterministic basis
  discovery.
- Reuse `fundamental_operators(h)` as the generator seed.
- Separate complete rule coverage from independent closure coordinates.
- Implement hooks for Fock/phase-space, Pauli/Spin, and local transition spaces.
- Implement collective hooks without pretending collective transitions obey local product
  identities, and require explicit space context when their dimension cannot be inferred.

Exit criterion: every currently supported named transform can validate through its Hilbert-space
methods without changing public behavior or introducing a parallel algebra hierarchy.

### Phase 2: shared exact map materializer

- Build forward/inverse rules from coefficient matrices and shifts.
- Centralize adjoint pairing, inverse construction, relations, and algebra validation.
- Migrate named `Rotation`, `Squeeze`, and explicit `Displace` constructors where this removes
  duplicated rule construction.
- Keep a direct constructor path when the generic materializer is measurably worse and sharing
  would not reduce maintained logic.

Exit criterion: named-map laws and performance remain within measured thresholds, with a net
reduction in independent rule-building logic.

### Phase 3: general bosonic quadratic-affine extraction

- Replace separate one-mode Hamiltonian scanners with a shared exact quadratic scanner.
- Share bosonic action and response logic between `FockSpace` and `PhaseSpace` hooks.
- Build structurally connected mode blocks.
- Retain scalar and `2×2` response solvers as specialized methods.

Exit criterion: all existing one-mode displacement behavior is unchanged and ladder/quadrature
forms agree under independently checked equivalent examples.

### Phase 4: multimode static and bounded displacement

- Add deterministic mode-list APIs.
- Solve coupled static equilibrium systems.
- Solve finite-harmonic response block by block.
- Support hopping, anomalous terms, and cross-quadrature couplings.
- Verify exact cancellation, resonance refusal, and Kerr/nonlinear full-model application.

Exit criterion: coupled Fock and quadrature references produce equivalent exact maps, and warm
construction costs have representative allocation gates.

### Phase 5: explicit Bogoliubov maps and conditionality

- Add `Bogoliubov` using the shared Heisenberg map materializer.
- Prove canonicality for numeric or structurally parametrized maps.
- Integrate visible symbolic constraints with the conditional-transform work in #234.
- Validate constraints before numerical backend conversion.

Exit criterion: #148 behavior is covered without a second rule-construction architecture.

### Phase 6: generic exact closed-adjoint exponentiation

- Extract finite affine actions for supported operator bases.
- Add the initial exact exponential strategy registry.
- Implement `UnitaryTransform(G, θ)` and `RotatingFrame(H0, t)` with exact refusal diagnostics.
- Compare generic results with every equivalent named constructor.

Exit criterion: representative bosonic, Pauli/Spin, and transition generators agree exactly with
their named transformations; unsupported closure and exponential cases fail clearly.

### Phase 7: exact dressed strategies

- Implement exact two-level mixing and already diagonal blocks.
- Add verified registered spectral structures.
- Add supported Gaussian symplectic/Bogoliubov diagonalizers only with explicit branch and
  degeneracy rules.

Exit criterion: every returned frame is independently verified to be canonical/unitary and to
diagonalize the requested Hamiltonian exactly.

### Phase 8: indexed and collective integration

- Connect the algebra-basis interface to indexed-family rule representation.
- Implement collective basis rotations using collective bracket laws.
- Preserve alpha-renaming invariance and complete coverage semantics from #234.

Exit criterion: indexed and collective maps reuse the same exact-map and diagnostic layers without
introducing wildcard dictionary semantics.

## Test laws

Prefer behavior and algebra laws over type-only coverage:

- forward followed by inverse is exactly the identity;
- adjoints transform consistently;
- the basis-specific bracket and product identities are preserved;
- displacements eliminate certified linear terms exactly;
- timed transforms include the complete operator and scalar gauge;
- named and generic derivations agree on shared supported cases;
- equivalent ladder and quadrature Gaussian problems agree numerically and symbolically where the
  coefficient domain permits;
- transformations compose correctly on disjoint and overlapping sites;
- a deliberately wrong inverse, noncanonical matrix, or incomplete basis is rejected;
- every exact spectral strategy verifies its own decomposition before returning.

## Documentation plan

- Keep named transformations as the approachable entry point.
- Explain the distinction between an operator map, a moving-frame gauge, and a Hamiltonian-derived
  solving objective.
- Document Gaussian centering, bounded response, evolution, and diagonalization as separate
  operations.
- Include multimode driven-resonator and Kerr examples.
- Show Pauli/Spin and transition examples using the same generic closed-adjoint API.
- State exactness domains, resonance surfaces, conditional constraints, and refusal cases next to
  each constructor.

## Issue ownership

- #231 owns the closed-adjoint engine, exact generator flows, bounded affine Hamiltonian solving,
  rotating frames, and dressed-frame strategies.
- #148 owns the public Bogoliubov spelling and its user-facing examples.
- #234 owns indexed families, collective map semantics, conditional canonicality, and structured
  diagnostics.
- #232 owns approximate and partial transformations.

One foundational implementation may advance several issues, but acceptance and review should stay
split along these semantic boundaries.

## Resolved design decisions

These decisions are the current default for implementation. Items in the following prototype-gate
section can still narrow support without changing the architecture.

### Hilbert-space context is an optional keyword

Keep the ergonomic inferred APIs and add `hilbert = h` only when complete space metadata cannot be
recovered from the operators:

```julia
UnitaryTransform(G, θ; hilbert = h)
RotatingFrame(H0, t; hilbert = h)
DressedFrame(H0; hilbert = h)
```

Keywords avoid a second family of visually competing positional APIs. The public keyword method
is a thin adapter: it calls private positional methods specialized on `Nothing` or the concrete
`HilbertSpace` type. Do not introduce an analyzed-model object until a benchmark demonstrates that
reusing extracted closure data is important.

Wave 1 further fixes the correctness boundary: explicit context is authoritative. Inference is a
checked convenience only when every touched site has complete algebra metadata. It cannot recover
untouched product factors, original space names, index ranges, symbolic level labels, or the level
count of collective transitions.

### Closure coordinates reuse `QTerm`

Use canonical `QTerm` equality, hashing, operator storage, and inequality constraints for private
closure coordinates. The minimal basis record stores only the concrete Hilbert space, a
`Vector{QTerm}`, and a `Dict{QTerm,Int32}`. A separate analysis result owns action coefficients,
affine offsets, work limits, and diagnostics, and is discarded after materialization.

Closure discovery checks commutator count, coordinate count, and maximum monomial degree in a
documented deterministic order. It reports the source coordinate and first unseen monomial. The
limit is never based on elapsed time. Concrete indexed sites are allowed; free indexed-family
closure remains deferred to its dedicated gate.

### Gaussian coordinates are site-interleaved with split solver views

Persist bosonic ladder coordinates as `(a1,a1',a2,a2',...)`, matching current generator ordering,
site coverage, and rule materialization. A Gaussian strategy may construct a temporary
split-Nambu permutation `(a1,a2,...,a1',a2',...)` to expose particle/hole blocks. Do not store a
second basis or change public generator ordering.

The ladder/quadrature conversion uses `a=(x+im*p)/sqrt(2)` and retains normal-order scalar
offsets. These signs and factors are pinned by exact commutator fixtures and an independent
numeric coordinate-change oracle.

### Product-space adjoint coordinates are flat

Store one deterministic flat coordinate vector and flat rule-generator vector. Record contiguous
per-site ranges for the initial local bases, plus the structural coupling graph used to find solve
blocks. A tuple of heterogeneous local basis objects would complicate global matrix indices,
cross-site Gaussian couplings, and materialization into the already flat rule dictionaries.

Closure elements involving several sites cannot belong to one local tuple. If later closed-adjoint
strategies admit such products, append their `QTerm` coordinates and record their site support
explicitly. Do not widen coordinates to `Vector{QField}` or another abstract container.

### Local transition closure uses the non-ground matrix units

For an `NLevelSpace` with ground level `g`, use the independent noncentral coordinates

```math
\{\sigma^{ij}: i \ne j\}\;\cup\;
\{\sigma^{ii}: i \ne g\},
```

with the identity represented separately. This gives `n² - 1` independent coordinates and agrees
with the existing completeness normalization

```math
\sigma^{gg}=1-\sum_{i\ne g}\sigma^{ii}.
```

Seed these coordinates with `fundamental_operators(h)` and their adjoints. Public rule coverage
still includes all `n²` matrix units; construct the ground-projector image from completeness and
verify that the complete rule set agrees with direct matrix conjugation.

Do not switch to Gell-Mann coordinates: they introduce additional normalization coefficients and
make rule materialization less direct without improving closure.

### Rational-function domains are not canonicality conditions

A map with coefficients such as `1 / det(K)` is an exact rational map on the domain where its
coefficients are defined. Keep this domain behavior on `UnitaryTransform`:

- reject structurally zero denominators immediately;
- retain symbolic quotients;
- document that they are valid away from their poles;
- do not turn every symbolic displacement into a `ConditionalTransform`.

Canonicality conditions are different: a raw Bogoliubov matrix with symbolic `u` and `v` is not
known to preserve the algebra until an equality such as `abs2(u)-abs2(v)=1` is established. Those
proof obligations belong on `ConditionalTransform`.

A future `domain_conditions(U)` may expose coefficient poles if the coefficient layer can extract
and simplify them reliably. Do not store solver-provenance conditions that become stale after
coefficient cancellation or composition.

### Generic exact solves avoid symbolic pivot divisions

Retain the existing explicit scalar and `2×2` formulas. For a larger structurally connected block,
prototype an adjugate solve based on the Faddeev-LeVerrier recurrence:

```math
x=\frac{\operatorname{adj}(A)b}{\det(A)}.
```

The recurrence divides only by exact integer step counts, not by symbolic pivots, and the returned
map has the single maximal invertibility condition `det(A) ≠ 0`. It avoids the artificial branch
conditions introduced by choosing symbolic Gaussian-elimination pivots.

Do not implement Bareiss elimination before the coefficient domain supports exact polynomial
division. A Julia MCP probe of the current `CNum` arithmetic showed that even
`((a*d-b*c)*k)/(a*d-b*c)` remains an unevaluated quotient rather than reducing to `k`, so Bareiss's
supposed exact intermediate divisions would accumulate rational expression noise.

Keep this solver private to the transformation work until a second coefficient-level consumer
needs it. Promote it to a shared exact-linear-algebra layer only after that reuse is concrete.

Wave 1 accepts Faddeev--LeVerrier only behind a strict boundary:

- decompose structural connected blocks before solving;
- keep scalar and explicit `2x2` paths;
- validate `A*numerator == det*forcing` before constructing quotients;
- initially cap connected blocks by a measured expression-work budget;
- map an isolated structural zero with zero forcing to zero;
- refuse nonzero forcing on a structural null coordinate and every connected structurally singular
  block, even when one numeric substitution happens to be compatible.

Do not require the completed quotient residual to simplify structurally. Cross multiplication is
the exact rational-equivalence oracle under current coefficient simplification.

### Initial spectral strategies are branch-safe and mostly root-free

The current coefficient domain can store symbolic fractional powers and can verify some squared
radicals, but it does not carry the sign and branch assumptions required to manufacture arbitrary
spectral roots safely. For example, the current symbolic pipeline reduces `sqrt(a^2)` to `a`
without an explicit sign condition.

Initially support:

- already diagonal blocks;
- named rotation and squeeze parametrizations whose identities are certified by relations;
- user-supplied exact unitary eigenbases that pass independent verification;
- blocks whose required root is already supplied with an explicit certificate;
- `2×2` blocks whose discriminant is a structurally recognized perfect square.

Defer general symbolic square-root construction, cyclotomic values, and higher algebraic roots to
the safe-root and exact-number-domain work in #231. Storing a radical expression is not itself a
proof that its branch makes the returned eigenbasis exact and continuous.

### Dressed frames use structural, not energy, ordering

Symbolic eigenvalues cannot generally be ordered. `DressedFrame` therefore uses deterministic
structural conventions:

- order disconnected blocks by the existing canonical site/operator ordering;
- preserve each exact strategy's declared bare-state ordering;
- do not reorder states by symbolic energy magnitude;
- label nondegenerate dressed states by their deterministic bare pivot or strategy position;
- leave a block unchanged when the Hamiltonian is structurally scalar on it;
- require a strategy or user-supplied exact basis for a nontrivial degenerate subspace;
- let each strategy define a deterministic algebraic eigenvector phase and verify it, rather than
  trying to make a symbolic component "positive".

Every strategy verifies `H*W == W*Diagonal(E)` and `W'*W == I`. Possible parameter crossings remain
documented rational or spectral domains; they do not trigger substitution-dependent reordering.

### Timed matrix APIs distinguish unitary matrices from operator-action matrices

- A static exact matrix map needs no gauge.
- A finite-level `W(t)` is the actual Hilbert-space unitary, including its phase, so the complete
  gauge `im*Ẇ'W` is derived exactly as in the existing timed transition rotation.
- Named bosonic rotations, squeezes, and Weyl displacements define their lift and derive their
  complete gauge.
- A frame constructed from an explicit Hermitian generator derives its gauge from that generator.
- An opaque moving symplectic/Bogoliubov operator-action matrix does not specify an arbitrary
  time-dependent scalar phase. Its timed API must either select and document a canonical lift plus
  scalar convention, accept an explicit scalar phase, or accept a complete caller-supplied gauge.

When a caller supplies a gauge, validate its Hermiticity and verify that its nonscalar adjoint
action reproduces the matrix derivative. The scalar part is a phase convention and cannot be
validated from the operator action alone.

### The shared materializer has an end-to-end replacement gate

Capture the current direct constructors before migrating them. Compare grouped, millisecond-scale
workflows so timer noise does not decide the architecture:

1. construct, compose, invert, and apply a one-mode displacement/rotation/squeeze chain to a Kerr
   Hamiltonian;
2. construct and apply one- and 33-harmonic bounded displacements;
3. construct and apply a coupled four-mode Gaussian map;
4. construct static and timed four-level matrix-unit rotations, including their gauge;
5. perform forward/inverse round trips and render the resulting transforms.

Replace a specialized rule builder only when the shared path clearly removes independent logic,
introduces no JET dispatch, stays within 10% of the warm median runtime, and stays within 20% of
warm allocations for its representative workflow. Otherwise keep the direct coefficient builder
as a specialized method feeding the shared validation and materialization boundary.

Wave 2 rejects both transient and retained rule-storage wrappers. The common terminal step scans
concrete coefficient rows directly into forward and inverse `Dict{Op,QAdd}` values and returns the
existing `UnitaryTransform`. Its inputs are short-lived `QTerm` coordinates, `CNum` maps/offsets,
and relations. A timed strategy passes its scalar lift positionally as `CNum`; the materializer
consumes it into the gauge and retains no lift object.

Keep small direct builders such as displacement specialized. The dense N-level matrix-to-rule loop
is the first plausible shared-materializer consumer, but migrate it only when its complete
construction and application workflow meets the same gates.

### Exact strategies use methods and ephemeral certificates

Dispatch on concrete analyzed block types. A supported method returns a concrete short-lived
certificate; an unsupported exact structure returns `nothing`, after which the caller formats a
deterministic refusal. Do not introduce a mutable registry, stored `Function`, or abstract-block
hot path.

The certificate contains only a zero-size strategy marker, forward map, exact inverse,
verification token, and static/timed scalar-lift token. It is consumed immediately by validation
and materialization. It is not a sibling transformation and never retains basis, action,
derivatives, diagnostic traces, or provenance.

Initial branch-safe strategies are diagonal atoms with exact reciprocals, named normalized
rotations and squeezes, normalized `+I`/`-I` involutions, nilpotent finite polynomials, and verified
user maps with exact inverses. Root-producing normalization and generic symbolic matrix
exponentiation remain refused.

All coefficient matrix operations use private explicit loops with `_CNUM_ZERO`, `_CNUM_ONE`, and
the internal arithmetic helpers. `Coeff` deliberately has no generic `zero`, so generic
`LinearAlgebra` allocation, multiplication, inversion, solving, and exponentiation are outside the
implementation contract.

### Conditional transformations are siblings of unitary transformations

Do not put unproven maps inside `UnitaryTransform`, and do not make `ConditionalTransform` a subtype
of `UnitaryTransform`. Introduce a common transformation interface and shared private rule storage:

```julia
abstract type AbstractTransformation end

struct UnitaryTransform{T} <: AbstractTransformation
    # proven exact rule data
end

struct ConditionalTransform{T} <: AbstractTransformation
    # the same formal rule data
    conditions::Vector{SymbolicCondition}
end
```

The final storage layout should factor common rule data without adding an access indirection to
hot unconditional application unless benchmarks show it is negligible. `ParamRelation` remains a
private simplification identity and is not reused as a logical proof obligation.

Formal conjugation, inversion, and composition may operate on conditional maps while preserving
and simplifying their conditions. Conversion to a numerical backend must prove the substituted
conditions first. When all conditions are proven, an explicit operation may promote the result to
`UnitaryTransform`; the public type must never imply that unresolved assumptions are true.

## Prototype gates still open

These gates are not all prerequisites for Phase 1. Run each immediately before the phase whose
architecture it can change. The immediate foundation blockers are closure-coordinate storage,
Hilbert-space inference, multimode ladder conventions, gauge reconstruction, singular-response
semantics, and the generic adjugate-solver experiment. The remaining gates protect later public
APIs and optimizations from premature abstraction.

### Foundation gates

- [x] Prototype closure discovery and projection before fixing the `AdjointBasis` coordinate type.
      Start from canonical monomials, repeatedly evaluate `i[G, X]`, expand transition
      completeness, and stop at exact closure or a deterministic work bound. Cover number
      rotation, one- and two-mode squeezing/rotation, Pauli/Spin rotation, and N-level mixing as
      accepted cases; use Kerr acting on a linear Fock basis and a spin-boson coupling as refusal
      cases. The prototype must determine whether a small concrete `BasisMonomial` record can reuse
      `QTerm` storage or needs its own identity/site-support fields. Do not settle on
      `Vector{Op}`—which cannot represent product closure elements—or `Vector{QField}`—which is
      abstract—before this spike.
- [x] Define and measure the closure work bound. Record basis size, maximum monomial degree,
      commutator count, and elapsed symbolic work. Refusal must report the first expanding
      commutator and the bound reached; it must not depend on wall-clock time.
- [x] Prototype Hilbert-space inference and name reconciliation for renamed operators, duplicate
      subspace types in a `ProductSpace`, concrete indexed sites, ordinary transitions with a
      non-default ground level, and collective transitions. Confirm exactly which inputs require
      `hilbert = h` and that a mismatched supplied space is rejected before analysis.
- [x] Implement the Faddeev-LeVerrier adjugate solve in an isolated prototype and measure symbolic
      expression size, warm allocations, and residual verification for coupled `4×4`, `8×8`, and
      structurally block-diagonal systems.
- [x] Decide whether the adjugate solver remains transformation-private or moves to the coefficient
      layer after identifying a second real consumer.
- [x] Prototype the flat product-space basis with a cross-mode Gaussian block and one finite
      cross-site closed-adjoint example to validate coordinate-support bookkeeping.

### Gaussian-affine gates

- [x] Fix the multimode ladder-coordinate convention with a coefficient-level prototype. Compare
      site-local ordering `(a₁,a₁†,a₂,a₂†,…)` with split Nambu ordering
      `(a₁,a₂,…,a₁†,a₂†,…)` for commutator matrices, adjoint permutation, block discovery, and rule
      materialization. Scan hopping, single- and two-mode anomalous terms, and complex drives once;
      verify the resulting action against direct commutators.
- [x] Establish ladder/quadrature equivalence on the same one- and two-mode quadratic models,
      including normal-order scalar offsets. This prototype owns the signs and factors of `im`,
      `1//2`, and `sqrt(2)` used by every later Gaussian solver.
- [x] Prototype gauge reconstruction from `(M(t), d(t), ∂M, ∂d)` under a named Weyl/metaplectic
      convention. Compare its complete quadratic, linear, and scalar gauge with the existing timed
      displacement, rotation, one- and two-mode squeeze, and beamsplitter constructors, then test
      inversion and composition. If one formula cannot reproduce the named phase conventions
      without extra lift data, keep gauge construction strategy-specific and share only its
      materialization and verification.
- [x] Prototype structurally singular but compatible affine systems. Distinguish: unique response;
      nonunique response with the null-space contribution fixed to zero; inconsistent forcing; and
      a resonant timed forcing with no bounded solution. Do not add symbolic rank branching: if
      compatibility is not structurally decidable, return the ordinary quotient for a structurally
      nonsingular determinant or refuse with the unresolved block identified.
- [ ] Measure dense coefficient storage against structural block storage before adding sparse
      matrix machinery. Include uncoupled 8-mode, nearest-neighbor 8-mode, and dense 4-mode
      references. The result decides whether analysis matrices stay dense per connected block;
      avoid a new sparse abstraction unless it wins a complete construction workflow.
- [ ] Benchmark application, not only construction, on displaced multimode Kerr and cross-Kerr
      Hamiltonians. Record transformed term count and allocations to ensure a general affine map
      does not introduce avoidable expression expansion in `_substitute_op_rules`.

### Exact-strategy and API gates

- [x] Prototype the exact-strategy interface without a registry of stored functions. Candidate
      strategies should be ordinary methods returning `nothing` or a concrete certificate that
      contains the map and verification data. Exercise diagonal, rotation/squeeze, involutive, and
      unsupported blocks; inspect method ambiguities, invalidations, and JET dispatch before
      declaring the interface extensible.
- [x] Define a minimal construction certificate and prove it is not redundant with
      `UnitaryTransform`. It should retain only information needed before materialization: the
      selected strategy, coefficient map, exact inverse, and verification result. Do not persist
      analysis provenance after rule construction unless diagnostics need it publicly.
- [x] Run the materializer replacement workflows and retain specialized builders wherever the
      10% runtime or 20% allocation gates fail.
- [ ] Specify `SymbolicCondition` only as part of #234, including equality/nonzero condition kinds,
      substitution, simplification, composition, and numerical proof behavior.
- [ ] Prototype shared rule storage for sibling `UnitaryTransform` and `ConditionalTransform`
      types. Compare direct fields with a private `RuleTransformData` wrapper on conjugation,
      composition, inversion, display, and numeric conversion; do not add an access indirection to
      the unconditional hot path without measurements.
- [x] Run API-shape checks for the `hilbert` keyword and new generic constructors: `@inferred`, JET,
      Aqua ambiguity/piracy checks, method invalidation inspection, fresh precompile/load, and
      Documenter `checkdocs = :exports` entries for every concrete exported signature.
- [x] Perform an indexed-family compatibility spike before freezing basis keys. The closed-adjoint
      engine may continue to reject free indexed families initially, but its site and coordinate
      identity must not prevent the alpha-equivalent family-rule representation planned in #234.
- [ ] Revisit root-producing dressed strategies only after a branch-aware root or assumption
      representation is available.

## Agentic prototype workflow

Use agents to falsify design choices and collect evidence, not to have several workers implement
competing production architectures. The coordinator owns the design file and integration branch.
Each specialist receives one bounded question, isolated files, common fixtures, and an explicit
report contract.

### Isolation and artifact rules

- Give every code-writing specialist its own git worktree. Agents sharing one working directory
  may perform read-only analysis concurrently, but must not edit overlapping production files.
- Keep exploratory code under an ignored prototype directory or on disposable spike branches.
- Never merge a prototype merely because it runs. Port only the smallest proven idea into
  production after the decision checkpoint.
- Preserve the user's pre-existing working tree and do not fold unrelated changes into a spike.
- Run Julia only through the Julia MCP so every result has a captured environment and completion
  status.
- Give every agent the same exact symbolic fixtures and numeric oracle points. This makes reports
  comparable and lets the coordinator reproduce disagreements.

Every prototype report contains:

```text
Question
Competing designs
Fixtures and refusal cases
Independent oracle
Correctness result
Inference/JET result
Runtime, allocations, and expression size
Production concepts and lines introduced
Failure modes
Recommendation: adopt / reject / investigate
Files safe to discard
```

### Wave 0: coordinator establishes the baseline

Before delegation, the coordinator:

1. records the branch, complete working-tree state, Julia version, and project manifests;
2. freezes behavior fixtures for every current named transform and automatic displacement;
3. records the existing complete-workflow benchmark and warm-allocation baselines;
4. prepares shared symbolic and numeric oracle cases;
5. gives each agent explicit file ownership and a refusal policy;
6. creates a decision matrix with one row per prototype gate.

The fixtures include creator-operator inputs, complete timed gauges, transition completeness,
multimode Kerr application, symbolic possible resonances, and exact structural resonances.

### Wave 1: three independent foundation spikes

Run three specialists concurrently while the coordinator maintains fixtures and reviews their
intermediate evidence.

#### Agent A: closure and Hilbert-space context

- Prototype `BasisMonomial` versus reused `QTerm` closure coordinates.
- Exercise exact closure, basis growth, identity residuals, and deterministic refusal bounds.
- Reconcile inferred operator metadata with explicit `hilbert = h` for ProductSpace, renamed,
  indexed, transition, and collective cases.
- Return a proposed minimal `AdjointBasis` layout without editing `UnitaryTransform`.

#### Agent B: Gaussian conventions and gauge

- Implement both candidate Nambu orderings in scratch code.
- Scan the same one- and two-mode Hamiltonians in ladder and quadrature coordinates.
- Compare actions with direct commutators and independent numeric matrices.
- Attempt common gauge reconstruction and compare it with all existing named timed constructors,
  including inversion and composition.
- Return the chosen ordering, sign/factor convention table, and whether lift-specific gauge data is
  unavoidable.

#### Agent C: exact response solvers

- Implement the Faddeev-LeVerrier adjugate prototype and specialized singular-compatible solves.
- Compare `2×2`, `4×4`, `8×8`, block-diagonal, and structurally singular fixtures.
- Record determinant size, result term counts, residual-verification cost, allocations, and warm
  runtime.
- Return a strict support boundary; a failed generic solver is an acceptable and useful result.

### Foundation decision checkpoint

Do not start production implementation until the coordinator can fill these rows with reproduced
evidence:

| Decision | Required evidence |
| --- | --- |
| Closure coordinate storage | accepted/refused closure fixtures and concrete-storage/JET results |
| Hilbert inference boundary | exact matrix of inferred versus explicit-context cases |
| Nambu ordering | direct-commutator equivalence and simpler materialization/block behavior |
| Gaussian gauge sharing | equality with every named gauge or a demonstrated counterexample |
| Singular semantics | unique/nonunique/inconsistent/resonant fixtures |
| Generic solver boundary | expression-size and end-to-end timing results |

If two reports disagree, assign one specialist a short adversarial reproduction using the other
agent's fixture before deciding. Update this design file with the measured result, then discard
superseded prototype code.

### Wave 2: architecture-pressure spikes

After Wave 1 fixes the foundation contracts, run three more specialists concurrently.

#### Agent A: materializer and hot-path storage

- Build the smallest coefficient-map-to-rules materializer compatible with the chosen basis.
- Compare direct fields with shared rule storage for unconditional and conditional transforms.
- Run the complete replacement workflows and identify which named builders remain specialized.

#### Agent B: exact-strategy interface

- Prototype method-based strategy selection returning concrete certificates or `nothing`.
- Cover diagonal, rotation/squeeze, involutive, user-basis, and unsupported blocks.
- Verify exactness, refusal messages, ambiguity, invalidation, and branch-safe root boundaries.

#### Agent C: product, indexed, and API pressure

- Exercise flat ProductSpace coordinates with cross-mode coupling and one finite cross-site closure.
- Check that basis keys can later represent alpha-equivalent indexed families without wildcard
  `Dict` semantics.
- Test the `hilbert` keyword API with inference, JET, Aqua, precompile, and documentation discovery.

### Wave 2 decision checkpoint

The coordinator selects the smallest combination that:

- preserves all exact behavior fixtures;
- carries no unresolved semantic assumption inside `UnitaryTransform`;
- introduces no abstract or closure payloads;
- meets the materializer runtime/allocation gates;
- has clear refusal boundaries;
- removes or centralizes real duplicated logic.

Record rejected alternatives and the fixture that rejected them. This prevents a later agent from
reintroducing an already disproven abstraction.

### First production vertical slice

Implement only enough shared architecture to migrate one representative path end to end:

1. infer or accept Hilbert-space context;
2. construct the chosen concrete adjoint basis;
3. extract one supported exact action;
4. solve or exponentiate it with one certified strategy;
5. materialize a complete static or timed `UnitaryTransform`;
6. compare it exactly with the existing named constructor;
7. run behavior, inference, allocation, JET, and documentation gates.

The first slice should reuse current one-mode displacement or a named rotation as its oracle. Do
not add multimode, conditional, indexed, and dressed public APIs simultaneously merely because the
internal representation can describe them.

### Three-agent review wave

Once the vertical slice is complete, give the full diff concurrently to three independent
reviewers:

1. **Physics and exactness reviewer:** commutator conventions, gauges, resonances, branch safety,
   inverse maps, and refusal semantics.
2. **Julia and performance reviewer:** inference, concrete storage, dispatch, allocations,
   invalidations, expression growth, and benchmark representativeness.
3. **API and maintainability reviewer:** reuse of Hilbert-space/operator infrastructure, duplicate
   state, constructor semantics, issue boundaries, documentation, and source complexity.

The coordinator aggregates every finding, fixes it or records a concrete reason to reject it, and
reruns the affected evidence. Reviewers do not merge their own suggested refactors.

### Verification and promotion

A production slice is promoted only after:

- focused unitary tests and all independent oracle comparisons pass;
- the complete test-project suite passes through the Julia MCP;
- JET and Aqua ambiguity/piracy checks pass without new allowlist entries;
- complete-workflow benchmarks meet the recorded gates;
- documentation builds with exported concrete signatures covered;
- fresh precompile/load succeeds;
- `git diff --check` passes;
- the design file records which prototype decisions are now verified rather than proposed.

The agentic loop stops when the current phase's acceptance criteria are met. It does not continue
adding strategies or abstractions simply because unused budget or available agents remain.

## Definition of architectural success

The work is successful when adding a new exact transformation normally requires:

1. implementing the required algebra hooks on an existing `HilbertSpace` subtype once;
2. supplying an exact action or solve strategy;
3. returning a certified coefficient-level map;
4. using the shared materializer, gauge, composition, inversion, and diagnostics;
5. adding behavior tests and a representative benchmark.

It should not require another independent implementation of rule coverage, adjoint pairing,
inverse rules, timed gauges, canonicality checks, or refusal policy.
