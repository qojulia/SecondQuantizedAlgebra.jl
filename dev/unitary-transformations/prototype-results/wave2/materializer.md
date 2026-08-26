# Wave 2A: coefficient-map materialization

## Question

What is the smallest type-stable bridge from the Wave 1 closure/Gaussian results to the
existing `UnitaryTransform`, and does a retained `RuleTransformData` wrapper reduce enough
code or cost to justify changing the hot representation?

The tested input contract was a pair of concrete affine maps:

- site-interleaved `Vector{Op}` generator rows;
- `Vector{QTerm}` closure coordinates (including nonlinear monomials if a future strategy
  needs them);
- dense concrete `Matrix{Coeff}` forward/inverse coefficients;
- concrete `Vector{Coeff}` scalar offsets;
- existing `Vector{ParamRelation}` certificates;
- for moving frames, an operator `QAdd` gauge plus a separately supplied scalar `CNum`.

The terminal value is always the existing `UnitaryTransform{StaticTime}` or
`UnitaryTransform{DynamicTime}`. No closure, solver, callback, provenance, coefficient
matrix, or scalar-lift object survives construction.

## Competing designs

1. **Direct terminal materializer (adopt).** Scan each coefficient row directly into a
   `QTermDict`, add its scalar offset, build `Dict{Op,QAdd}` forward/inverse rules, and pass
   those dictionaries to `_static_transform`/`_timed_transform`.
2. **Transient `RuleTransformData` wrapper (reject).** Materialize the same dictionaries,
   package rules/generators/sites/relations in a private immutable record, then unpack it
   into the existing transform.
3. **Retained nested wrapper (reject).** Store rules/generators/sites/relations behind a
   `data` field and reproduce conjugation, inversion, composition, and display through the
   extra field access.
4. **Named builders only (retain as front ends, not as terminal storage paths).** Existing
   `Displace`, `Rotation`, and `Squeeze` remain useful validated strategies. They can emit
   the common coefficient-map contract where that removes duplicate rule-building code.

## Fixtures and refusals

Accepted fixtures:

- one-mode displacement (2 generators, 4 forward/inverse terms including offsets);
- one-mode rotation/squeeze chain (2 generators, 4 terms in each direction);
- coupled two-mode rotation/squeeze chain (site-interleaved `a,a',b,b'`, 16 terms in each
  direction);
- dense exact three-level rotation using the rational orthogonal matrix
  `[1 2 2; 2 1 -2; -2 2 -1]/3` (9 generators and 81 terms in each direction).

Refusal probes were deterministic:

- a coefficient matrix with the wrong number of rows throws
  `DimensionMismatch("generator rows")` before materialization;
- a rule map containing only `x` for a phase-space site delegates to the existing validator
  and throws `ArgumentError("incomplete rule set: `x` has no rule for its conjugate variable")`;
- the existing validator continues to reject empty maps, mismatched forward/inverse keys,
  incomplete Fock/spin/transition sites, and free indexed generator families.

## Independent oracle and correctness

The existing named constructors were the independent rule oracle. For all four fixtures,
the generic materializer reproduced exact forward dictionaries, inverse dictionaries,
generator order, and gauge. Applying named and generic transforms to nontrivial quadratic
Hamiltonians produced exactly equal `QAdd`s.

Every generic fixture also passed the independent inverse law
`simplify(conjugate(conjugate(H,U),inv(U))) == simplify(H)`. The dense N-level fixture
exercises all 81 coordinate coefficients rather than a permutation-only shortcut.

Numeric conversion was checked after substituting real values into the coupled two-mode
Hamiltonian. Direct and nested paths produced equal 16-by-16 sparse QuantumOpticsBase
operators. This confirms that the generic terminal storage does not perturb the existing
numeric pipeline.

The timed probe supplied the scalar lift as one `CNum`; after construction it existed only
as the empty-operator entry of `U.gauge::QAdd`. No scalar-lift metadata or callback was
retained.

## Inference and JET

`@inferred` succeeded for:

- `materialize_direct(::AffineRulePair, ::Vector{ParamRelation})` returning
  `UnitaryTransform{StaticTime}`;
- `_materialize_rules(::CoefficientRuleMap)` returning `Dict{Op,QAdd}`;
- nested conjugation returning `QAdd`.

Targeted JET correctness reported zero errors. JET optimization reported nine SymbolicUtils
dispatches in the materializer. The existing `_rule_qadd` reports the same nine dispatches
for a symbolic coefficient, so the prototype adds no new optimization-report class and
would require no allowlist expansion. For comparison, construction of the named coupled
rotation/squeeze chain reported 65 optimization reports because it also performs the
symbolic trigonometric/hyperbolic coefficient work.

## Runtime, allocations, and expression size

Warm medians used 21 single-evaluation BenchmarkTools samples. Absolute timing varied with
concurrent prototype load, so allocation counts are the stronger evidence. The final run
was:

| fixture | named ns / bytes / allocs | direct map ns / bytes / allocs | wrapper ns / bytes / allocs |
|---|---:|---:|---:|
| displacement | 15,435 / 12,720 / 120 | 7,124 / 10,432 / 66 | 7,124 / 11,056 / 73 |
| one-mode chain | 160,286 / 54,400 / 815 | 7,124 / 10,624 / 66 | 7,124 / 11,248 / 73 |
| coupled chain | 84,299 / 85,840 / 920 | 13,061 / 20,320 / 148 | 14,248 / 21,360 / 158 |
| dense 3-level | 1,163,556 / 202,352 / 5,108 | 55,803 / 45,344 / 434 | 56,991 / 47,040 / 442 |

These numbers isolate final materialization: the coefficient maps were already solved, while
named builders also compute their analytic coefficients and, for chains, compose transforms.
They therefore prove that a shared generic terminal step does not miss the 10% runtime or
20% allocation thresholds; they do not claim that generic analysis is faster end to end.
Named builders should remain as strategy/front-end fast paths.

The transient wrapper saved nothing and added 624--1,696 bytes and 7--10 allocations. The
retained nested candidate had identical allocations to direct storage for conjugation,
inversion, display, and numeric conversion; coupled composition differed by only two
allocations. Runtime ratios fluctuated around parity, with no repeatable advantage. The
numeric workflow was 4.89 ms / 1,318,080 bytes / 32,802 allocations direct versus 4.87 ms /
the same bytes and allocations nested.

Final forward/inverse expression sizes were:

| fixture | generators | terms each direction | compact-show characters |
|---|---:|---:|---:|
| displacement | 2 | 4 | 52 |
| one-mode chain | 2 | 4 | 145 |
| coupled chain | 4 | 16 | 399 |
| dense 3-level | 9 | 81 | 1,586 |

## Production concepts and lines

- Keep `UnitaryTransform`'s direct concrete fields unchanged
  (`src/algebra/unitary.jl:23-40`). They already serve conjugation, inversion,
  composition, generators, gauge access, and printing without retained analysis state.
- Put one private coefficient-row-to-`QAdd` helper beside the current rule primitives
  (`src/algebra/unitary.jl:514-529`). It should accept only concrete `QTerm`/`Coeff`
  containers and use `_addto_key!`/`_addto!` directly.
- Feed the resulting dictionaries through the existing central validation and canonical
  site ordering (`src/algebra/unitary.jl:112-143`), then use the existing timed conversion
  (`src/algebra/unitary.jl:153-159`). The coupled fixture confirmed canonical persistent
  order `a,a',b,b'`.
- Leave application, inversion, and composition unchanged
  (`src/algebra/unitary.jl:195-229,249-284,449-511`). They consume exactly the terminal
  representation produced by the generic materializer.
- Refactor named builders only where they duplicate substantial matrix-to-rule loops. Tiny
  builders such as displacement (`src/algebra/unitary_constructors.jl:14-18`) should keep
  their direct fast path unless an end-to-end benchmark shows that routing through a dense
  coefficient matrix pays for itself. The N-level rule loop
  (`src/algebra/unitary_constructors.jl:828-848`) is the clearest first consumer.

## Failure modes and boundaries

- Dense `Matrix{Coeff}` scanning is appropriate for the tested closed bases, but storage and
  scan cost are quadratic in closure dimension even when an action is sparse. Do not add a
  sparse alternative until a representative strategy crosses the agreed thresholds.
- `QTerm` contains vectors, so every inserted dictionary key must be copied. Sharing closure
  keys directly would make hashes vulnerable to later mutation.
- Dense N-level basis rotations intrinsically produce `n^4` rule terms. The 81-term 3-level
  fixture already makes materialization visible, but the analytic matrix action—not the
  terminal transform wrapper—is the scaling concern.
- A generic materializer must not infer or repair missing adjoints/sites. Strategy-specific
  validation produces the map; `_validated_transform` remains the final safety boundary.
- The scalar gauge is not derivable from the affine operator map. Require the strategy to
  pass it explicitly as `CNum`, consume it immediately into the gauge `QAdd`, and reject any
  design that retains a callback or provenance object.
- The nested prototype's composition intentionally exercised a non-diagonal coupled map.
  Production must continue using the existing diagonal/invariant-gauge fast paths.

## Decision

**Adopt** a small private direct coefficient-map materializer returning the existing
`UnitaryTransform`. Use `QTerm` coordinates, concrete `Coeff` matrices/vectors, current
canonical generator/site validation, and a positional scalar `CNum` gauge handoff.

**Reject** both transient and retained `RuleTransformData` wrappers. The transient form adds
allocations; the retained form changes every hot consumer while showing no allocation or
runtime benefit. A short-lived strategy-specific coefficient-map record is sufficient.

**Retain** named public builders as validated front ends and measured specializations. Route
only substantial common matrix-rule materialization through the generic helper, beginning
with the future exact closed-adjoint/N-level path; preserve tiny direct builders until an
end-to-end benchmark justifies changing them.

## Files safe to discard

- `prototype/materializer_prototype.jl`
- `prototype/run_materializer_prototype.jl`
- `prototype/jetenv/Project.toml`
- `prototype/jetenv/Manifest.toml`
- the ignored benchmark `Manifest.toml` created solely to run this isolated checkout

No production source, test, documentation, or benchmark file was edited.

