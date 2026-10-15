```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Parameters and coefficients

The coefficients of an operator expression are scalars: exact numbers or symbolic expressions from [Symbolics.jl](https://symbolics.juliasymbolics.org/stable/). They multiply operators without touching the operator algebra, and [`simplify`](@ref) applies symbolic identities to them.

```@contents
Pages = ["coefficients.md"]
Depth = 2:3
```

## Symbolic parameters

Parameters are declared with `@variables`, which is re-exported from Symbolics. A plain declaration is real, so it is its own complex conjugate.

```@example params
using SecondQuantizedAlgebra
a = Destroy(FockSpace(:cavity), :a)
@variables ω η::Number ζ::Complex

(ω * a)'
```

A parameter declared `::Number` is complex but stays a single symbol.

```@example params
(η * a)'
```

A parameter declared `::Complex` is split into independent real and imaginary parts.

```@example params
(ζ * a)'
```

In summary:

| Declaration | `conj` | Representation |
| --- | --- | --- |
| `@variables ω` (or `ω::Real`) | `ω` | one real symbol |
| `@variables η::Number` | `conj(η)` | one complex symbol |
| `@variables ζ::Complex` | `real(ζ) - im*imag(ζ)` | independent real and imaginary parts |

A complex amplitude such as a drive is best declared `::Number`, which keeps every product to a single symbol. `::Complex` carries both parts through every product, and suits a parameter whose real and imaginary parts are wanted as independent unknowns.

## Unit phases

Time-dependent drives and rotating frames produce phase factors ``e^{i\phi}``. These are written with [`expim`](@ref), which keeps them exact. The argument has to be real, since the identities ``\overline{e^{ix}} = e^{-ix}`` and ``|e^{ix}| = 1`` hold only for real ``x``.

```@docs
expim(::Real)
```

The same phase is recognized however its argument is written.

```@example phases
using SecondQuantizedAlgebra
using SecondQuantizedAlgebra: exponential_form, trigonometric_form
import Symbolics
@variables ω t θ r
a = Destroy(FockSpace(:cavity), :a)

expim((ω + 2θ) * t) == expim(ω * t + 2θ * t)
```

Integer powers, negative ones included, and division also stay a single phase.

```@example phases
expim(θ)^3 / expim(ω * t), expim(θ)^-2
```

A fractional power would have to pick a branch, so `sqrt(expim(θ))` and `expim(θ)^(1//2)` have no method rather than being rewritten as `expim(θ/2)`.

Substitution rebuilds the phase from its new argument.

```@example phases
substitute(expim(ω * t), Dict(ω => 2ω))
```

Arguments are compared algebraically, not modulo ``2\pi``, so `expim(θ + 2π)` is not recognized as `expim(θ)`.

Differentiation follows the chain rule through the argument.

```@example phases
Symbolics.derivative(expim(ω * t), t)
```

A coefficient that is exactly one unit phase has the familiar trigonometric projections.

```@example phases
real(expim(ω * t)), imag(expim(ω * t)), abs(expim(ω * t)), abs2(expim(ω * t))
```

A phase keeps the form it was written in, so positive and negative frequencies are not merged into a cosine.

```@example phases
(expim(ω * t) + expim(-ω * t)) * a
```

Switching between the exponential and the trigonometric form is explicit, in both directions. Neither conversion is applied by display or by [`simplify`](@ref), and `expim(cos(θ))` stands for ``e^{i\cos\theta}``, not for a conversion of the cosine.

```@docs
exponential_form
trigonometric_form
```

Besides cosines and sines, `exponential_form` also recognizes `exp(im * θ)` and `cis(θ)`.

```@example phases
exponential_form(cis(θ)) * a
```

A coefficient built from phases splits into its Fourier components.

```@docs
phase_terms
PhaseTerm
```

## Coefficient identities

Symbolic identities are applied only by [`simplify`](@ref). It covers trigonometric and hyperbolic identities, also with composite arguments, and leaves a coefficient unchanged when an identity cannot be applied without blowing up the expression.

```@example phases
simplify((cosh(r)^2 - sinh(r)^2) * (a + a'))
```

## Substitution

Parameters and operators are replaced in a single call, with all rules applied at once. Shifting a mode by a coherent amplitude while fixing a frequency gives the displaced Hamiltonian.

```@example params
@variables α::Number
substitute(ω * a' * a, Dict(a => a + α, ω => 2))
```

```@docs
substitute(::QSym, ::AbstractDict)
```

## Taking expressions apart

To hand an expression to other code, for instance to build a numerical function from its coefficients, its parts can be read off directly. An expression first splits into its terms, and each term into a coefficient and operators.

```@example params
using SecondQuantizedAlgebra: sorted_arguments
terms = sorted_arguments(ω * a' * a + η * a)
```

```@example params
get_prefactor(terms[2]), get_operators(terms[2])
```

```@docs
sorted_arguments
get_prefactor
get_operators
get_variables(::QAdd)
```

Iterating an expression gives its terms as pairs of an operator product and a coefficient, kept in exact form. `to_num` converts a coefficient to Symbolics.

```@example params
using SecondQuantizedAlgebra: to_num
[(term.ops, to_num(c)) for (term, c) in 3ω * a' * a]
```

```@docs
QTerm
Coeff
to_num
```

[`coefficient_terms`](@ref) splits a coefficient into its terms, each an exact or floating scalar times powers of parameters, without going through Symbolics.

```@example params
using SecondQuantizedAlgebra: coefficient_terms
coefficient_terms(3ω^2 - 1 // 2)
```

```@docs
has_coefficient_terms
coefficient_terms
CoefficientTerm
```


A complex amplitude declared as `@variables β::Number` stays a single atom next to a phase,
and its conjugate is the symbolic `conj(β)`. Prefer it over `::Complex`, which splits the
parameter into `real(β) + im*imag(β)` in every coefficient:

```@example pr305
using SecondQuantizedAlgebra
h = FockSpace(:harmonic_example)
@qnumbers a::Destroy(h)
@variables ω::Real t::Real
@variables β::Number
β * expim(ω * t) * a + conj(β) * expim(-ω * t) * a'
```

## Harmonic projection

[`harmonics`](@ref) splits an expression into its Fourier components in a base frequency
`ω`, the step behind a rotating-wave or secular approximation. [`harmonic`](@ref) returns one
component, the coefficient of `expim(n * ω * t)`, and [`secular`](@ref) the static one:

```@example pr311
using SecondQuantizedAlgebra
h = FockSpace(:harmonic_example)
@qnumbers a::Destroy(h)
@variables ω::Real t::Real
import SecondQuantizedAlgebra: harmonic, harmonics, secular
@variables ω0::Real ωd::Real F::Real τ::Real
H = ω0 * a' * a - F * (a + a') * (expim(ωd * τ) + expim(-ωd * τ)) / 2
Hrot = transform(H, Rotation(a, ωd * τ, τ))
harmonics(Hrot, ωd, τ)
```

```@example pr311
secular(Hrot, ωd, τ)
```

The decomposition is exact: summing `harmonic(x, ω, t, n) * expim(n * ω * t)` over the
returned orders gives `exponential_form(x)`. A phase that is not an exact integer multiple of
`ω * t` plus a constant, or a time dependence outside a phase, throws an `ArgumentError`
rather than being dropped, because assigning it would assume a relation between frequencies.
The same functions act on [`SuperOperator`](@ref)s and [`Lindbladian`](@ref)s.



```@docs
harmonics
harmonic
secular
```
