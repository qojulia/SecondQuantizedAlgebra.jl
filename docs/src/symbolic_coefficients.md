# Symbolic Coefficients and Phases

The coefficients of a quantum expression can contain symbolic parameters, complex
amplitudes, trigonometric functions, and exact unit phases. They participate in operator
arithmetic without changing the operator algebra.

Symbolic parameters are created with [`@variables`](@ref), which is re-exported from
Symbolics.jl:

```@example coefficient-basics
using SecondQuantizedAlgebra

h = FockSpace(:resonator)
@qnumbers a::Destroy(h)
@variables ω g t

H = ω * a' * a + g * (a + a')
```

Variables are real unless another symbolic type is specified. This means `ω` is unchanged
by conjugation. Complex-valued amplitudes can be declared with `::Number`:

```@example coefficient-basics
@variables η::Number

drive = η * a + conj(η) * a'
```

See [Symbolic parameters](implementation.md#Symbolic-parameters) for the distinction between
real, atomic complex, and split complex variables.

## Exact unit phases

[`expim(x)`](@ref expim) represents the unit phase ``e^{ix}`` for a provably real argument.
It stays compact under multiplication, integer powers, inversion, conjugation,
substitution, and differentiation.

```@example phase-representations
using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: expim, exponential_form, trigonometric_form

h = FockSpace(:phase_example)
@qnumbers a::Destroy(h)
@variables θ ω t

p = expim(ω * t)
p * conj(p)
```

Integer powers remain exact, and opposite phases cancel exactly:

```@example phase-representations
p^3 * inv(p)^3

expim(ω * t) * expim(-ω * t) * a
```

Only real arguments are accepted because the identities
``\overline{e^{ix}}=e^{-ix}`` and ``|e^{ix}|=1`` require real ``x``. A complex symbolic
argument should instead be represented with the ordinary symbolic `exp` function.

## Choosing an exponential or trigonometric form

Unit phases are printed as exponentials, including when positive- and negative-frequency
terms occur together:

```@example phase-representations
(expim(ω * t) + expim(-ω * t)) * a
```

This preserves the representation selected by the user. Conversion between phase and
trigonometric forms is available explicitly.

[`exponential_form`](@ref) rewrites algebraic occurrences of `cos` and `sin` using
`expim`:

```@example phase-representations
exponential_form(cos(ω * t) * a)
```

[`trigonometric_form`](@ref) performs the reverse conversion:

```@example phase-representations
trigonometric_form((expim(ω * t) + expim(-ω * t)) * a)
```

Both functions act termwise on quantum expressions and leave their operator products
unchanged. They are representation changes: ordinary display and [`simplify`](@ref) do not
invoke them automatically. In particular, `expim(cos(θ))` means ``e^{i\cos\theta}``; it is
not a request to convert a cosine.

## Automatic coefficient identities

[`simplify`](@ref) reduces exact trigonometric and hyperbolic identities in operator
coefficients:

```@example coefficient-simplification
using SecondQuantizedAlgebra

h = FockSpace(:simplification_example)
@qnumbers a::Destroy(h)
@variables ω t r

simplify((cos(ω * t)^2 + sin(ω * t)^2) * a' * a)
```

```@example coefficient-simplification
simplify((cosh(r)^2 - sinh(r)^2) * (a + a'))
```

Composite arguments such as `ω*t` are supported. Reduction is exact and conservative: if
an identity cannot be applied safely without excessive expression growth, the original
coefficient is returned unchanged.
