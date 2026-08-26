# Wave 2 architecture-pressure report

Captured: 2026-08-26, Europe/Zurich.

All three disposable harnesses were reproduced by the coordinator through the Julia MCP. No
production files were copied from the spikes.

## Accepted boundary

The construction pipeline is short-lived analysis followed by the existing terminal type:

```text
explicit Hilbert context + flat QTerm basis
    -> concrete analyzed block
    -> method-selected exact certificate or nothing
    -> direct coefficient-row materialization
    -> existing UnitaryTransform
```

The certificate retains only its zero-size strategy marker, exact forward and inverse maps,
verification token, and static/timed scalar-lift token. It is discarded immediately. Basis,
action, derivatives, diagnostics, provenance, and lift metadata never enter `UnitaryTransform`.

## Rejected alternatives

- A transient `RuleTransformData` wrapper added 624--1696 bytes and 7--10 allocations across the
  tested materializations.
- Retained nested rule storage gave no allocation advantage in conjugation, inversion, display,
  composition, or numeric conversion and would change every hot consumer.
- A strategy registry or stored `Function` adds open-world state without capability beyond method
  dispatch.
- Generic `LinearAlgebra` multiplication, `zeros`, `inv`, `exp`, and `backslash` are invalid for
  `CNum`; private kernels must use package coefficient operations.
- Wildcard `QTerm` equality would corrupt ordinary coordinate identity. Future indexed families
  need an explicit binder-aware key in #234.

## Support boundary

- First branch-safe exact strategies: diagonal atoms with supplied reciprocals; named normalized
  rotations/squeezes; normalized plus/minus involutions; nilpotent finite polynomials; and verified
  user maps with exact inverses.
- Root-producing normalization and generic symbolic spectral exponentiation remain refusals.
- Flat ProductSpace coordinates support cross-mode Gaussian and finite cross-site closures.
- Explicit Hilbert context is required initially. Unindexed operators and concrete indexed sites
  are accepted; free indexed families are refused.
- Strategy-supplied scalar `CNum` is consumed into the complete gauge. Body velocity supplies only
  the nonscalar gauge.

## Production migration rule

Add one private direct coefficient-row-to-`QAdd` materializer beside the current rule primitives.
Keep `UnitaryTransform` fields and hot operations unchanged. Keep tiny named builders specialized;
only migrate a substantial matrix-to-rule loop after its complete constructor/application workflow
meets the 10% runtime and 20% allocation gates.
