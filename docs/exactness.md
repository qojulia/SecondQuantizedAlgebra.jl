# Exactness

What is exact today, where a floating-point value gets manufactured, and what each change
would buy. The goal is that a coefficient the user wrote exactly stays exact through the
algebra, and that the package never invents a float of its own.

## The policy today

Coefficients live on three tiers: `Native` (a `ComplexF64`), `Poly` (a sparse sum of monomials
over parameter atoms, with `ComplexF64` scalars), and `Complex{Num}` (the symbolic fallback).
`_to_cnum(x::Real)` keeps a value native only when `ComplexF64(x) == x`, so a value that cannot
round trip stays symbolic instead of silently truncating.

Two consequences catch people out.

**Dyadic rationals print as decimals.** `1//4` round trips losslessly, so it is stored native
and displays as `0.25`; `1//3` does not, so it stays symbolic and displays as `1//3`. Both are
exact. The spelling tracks which tier the value landed on, not its precision.

**The native tier is not exactness preserving.** Round tripping the *inputs* says nothing about
the *result*. Three exactly representable values can still cancel:

```julia
julia> a = 2.0^-40; b = 2.0^40;      # both exact in Float64
julia> (a + b) - b                   # exact answer is 2.0^-40
0.0
```

So "native" means "representable", not "safe". Any exactness guarantee has to come from an
exact scalar, not from the guard alone.

## Where numerics enter

| Source | What happens | Class |
|:--|:--|:--|
| `Monomial.scalar::ComplexF64` | every polynomial tier scalar is a float, so a non dyadic rational rounds on arrival | manufactured |
| native tier arithmetic | cancellation and overflow, as above | manufactured |
| `_rational_power` on a numeric base | `sqrt(Num(2))` becomes `1.4142135623730951` instead of an exact radical | manufactured, localized |
| `_common_factor` | returns a `Matrix{ComplexF64}`, and requires `_is_native` ratios, so it refuses exact ratios that do not round trip | manufactured, localized |
| `eigen` on an unstructured block | float eigenvectors and eigenvalues | manufactured, needs structure to avoid |
| `_sinpi_exact` table cap | declines outside denominators 1, 2, 3, 4 and 6, so the block falls back to `eigen` | coverage |
| `_snap_spectrum!`, `_scaled_rate` | repair float artefacts (`±√2` returned as two non negative floats, a unit magnitude scale minting a second atom) | consequence |
| `_eigen_weight` float method | drops weights under `1e-12` because rounding noise is not a structural zero | consequence |
| `is_canonical`'s `atol`, `_sup_norm`, `_bounded_by` | wave through residuals that cancel only to rounding | consequence |
| `_eigen_or_throw` Hermiticity test | compares against `1e-10 * scale` | consequence |
| `to_numeric` and the numeric backends | floats on purpose, at the boundary | by design, not a target |
| casus irreducibilis, a float the user wrote | no exact form exists to produce | fundamental |

Three routines are already exact and are the model to copy: `_involution_scale` tests
`C² = κI` by exact cancellation rather than a tolerance, `_exact_root` declines rather than
rooting numerically, and `_mixing_angle` leaves a numeric two level angle as an unevaluated
`atan(0.6, 2.0) / 2`.

The rows marked *consequence* are worth reading together. `_snap_spectrum!`, `_scaled_rate`,
the `1e-12` cut, `_sup_norm`, `_bounded_by` and `atol` exist for no reason other than that
floats are present. None of them is a feature. All of them are deletable once the values they
guard are exact, so the work below pays for itself partly in removed code.

## The keystone: an exact polynomial scalar

`Monomial.scalar` being a `ComplexF64` is the single fact that caps everything else. Until it
changes, no upstream exactness survives: perfectly exact algebraic eigenvector entries are
rounded the moment they become a coefficient.

It is also what makes exactness and phase cancellation mutually exclusive for a radical
coefficient. A phase cancels against its conjugate only on the polynomial tier, where it is an
atom carrying an exponent, while an exact non dyadic rational survives only off that tier.
Admitting a radical as a polynomial atom restores the cancellation but lowers `g√2` to a float
and respells `expim` as `cos`/`sin`; leaving it symbolic keeps the value exact but the residual
stalls at `-1/4 + (1/4)·exp(ix)·exp(-ix)`. Nothing bridges that from outside the coefficient
core, which is why the frame paths currently take the numeric spectrum whenever their weights
carry a radical.

The scalar cannot simply become a `Rational`. Every `Float64` is a dyadic rational, but adding
two of them overflows `Rational{Int}` at denominators around `2^110`, and a user who writes
`0.3` should keep a float rather than acquire a 17 digit fraction. The shape that works mirrors
the tier policy already one level up: a concrete struct holding an exact
`Complex{Rational{Int}}`, a float, and an `exact` flag, with contagion (exact combined with
exact stays exact, anything touching a float becomes float) and demotion to float on rational
overflow. The field stays concretely typed, so `CheckConcreteStructs` is satisfied.

Call sites to convert: `_canonical_terms!`, the polynomial arms of `_add_cnum` and `_mul_cnum`,
`_conj_poly`, `to_num`, `_phase_poly_bound`, `_to_complex`, `_exact_root`, `_split_scalar`, and
`Monomial` equality and hashing. The risk is not correctness but speed: `Monomial` is the hot
path that the `@allocations` tests pin, so this wants benchmarks beside it rather than a
follow up.

## What each stage buys

| Stage | Change | Unlocks | Removes |
|:--|:--|:--|:--|
| 0 | `_rational_power` emits an exact radical for a numeric base | `sqrt(Num(2))` stops becoming a float | |
| 0 | structure and Hermiticity detected on the exact `CNum` matrix instead of the float reduction | blocks whose ratios do not round trip | the `1e-10` Hermiticity tolerance |
| 1 | exact polynomial scalar | non dyadic rationals stay exact; a radical coefficient can be exact *and* cancel phases, so the frame gate goes away and frames become exact wherever the dressed basis is | `_snap_spectrum!`, `_scaled_rate`, the `1e-12` cut, `_sup_norm`, `_bounded_by`, `_bounded_after_rewrite`, `atol` |
| 2 | cyclotomic values in place of the `sinpi` table | every chain and ring size, not the current jagged set | `_RootVal`, `_sinpi_exact`, `_cospi_exact` |
| 3 | more structure detectors (Kronecker products, Hadamard, permutation by cycle, the 2×2 to 4×4 formulas) | further blocks that currently reach `eigen` | |

Stage 1 must come before stage 2. Cyclotomic arithmetic would compute exact values that the
polynomial tier then rounds, so doing it first is wasted.

Stage 2 is the reason the current table is a stopgap rather than something to extend. Widening
`_RootVal` from one root to a *sum* of roots adds only denominators 10 and 12, because
`sinpi(1//5)` is a nested radical rather than a sum, which in sizes means one extra chain
(`n = 11`) and one extra ring (`n = 24`). The right object is a cyclotomic number: every
`sin(kπ/m)` and `cos(kπ/m)` lies in `ℚ(ζ₂ₘ)`, which closes under multiplication, so the field
replaces the table instead of lengthening it. Cyclotomics.jl and CyclotomicNumbers.jl are both
MIT licensed and usable, or the same arithmetic is a sparse exponent vector reduced modulo the
cyclotomic polynomial.

## What stays out of reach

**Casus irreducibilis.** A 3×3 real symmetric matrix with three real eigenvalues has no
expression in *real* radicals, and the complex cube roots that a general formula produces would
poison every residual. Abel Ruffini is a separate and weaker obstruction that only starts at
`n = 5`. For an unstructured block like this the correct behaviour is to refuse, not to return
floats.

**Genuinely transcendental values.** A two level frame's dressed energies need a `sqrt` that no
`ParamRelation` closes, and a mixing angle is an `atan`. These can be carried unevaluated, which
is what the package does, but they will not reduce.

**A float the user wrote.** `0.3` is a float and stays one. Write `3//10` for an exact
coefficient. Scaling a Hamiltonian to integers is often enough to collapse a frame onto a
single phase atom.

## The invariant worth adopting now

> The algebra never manufactures a floating point value the user did not write.

This is testable, and it is the discipline that separates declining from guessing. The
structured eigendecomposition already honours it: it emits exact values or returns `nothing` so
the caller can fall back visibly. `_rational_power` on a numeric base and `_common_factor`'s
conversion both violate it today, and both are stage 0 items. Stating it as a rule stops further
numerics from leaking in while the keystone is decided, and it is what argues for an explicit
fallback over the silent `Complex(cos(θ), sin(θ))` that comparable packages return once their
own closed form table runs out.
