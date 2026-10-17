```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Algebra and ordering

The Heisenberg equation of a Kerr oscillator, ``\dot a = i[H, a]`` with ``H = K a^{\dagger 2} a^2``, asks to move an annihilation operator past two creation operators. Using ``[a, a^\dagger] = 1`` gives ``\dot a = -2iK a^\dagger a^2``; the same result comes out in normal order directly:

```@example kerr
using SecondQuantizedAlgebra
a = Destroy(FockSpace(:cavity), :a)
@variables K
1im * commutator(K * a' * a' * a * a, a)
```

Operator expressions are kept in a canonical form at all times. Every product is reordered and reduced the moment it is formed, so like terms combine without further work.

```@contents
Pages = ["algebra.md"]
Depth = 2:3
```

## Canonical form

Each algebra brings its own rules, and all of them apply as soon as a product is formed.

```@example canonical
using SecondQuantizedAlgebra
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 3) ⊗ PauliSpace(:qubit) ⊗ SpinSpace(:S) ⊗ PhaseSpace(:osc)
@qnumbers a::Destroy(h, 1) b::Destroy(h, 1) x::Position(h) p::Momentum(h)
σ(i, j) = Transition(h, :σ, i, j, 2)
τ(k) = Pauli(h, :τ, k)
S(k) = Spin(h, :S, k)
nothing # hide
```

Operators of different subsystems are sorted by their position in the product space.

```@example canonical
x * a * σ(2, 1)
```

Each bosonic mode is normal ordered, and modes with different names are grouped by name.

```@example canonical
b' * a' * b * a
```

Transitions on the same system compose.

```@example canonical
σ(1, 2) * σ(2, 3)
```

Products of Pauli matrices reduce to a single Pauli matrix or a number.

```@example canonical
τ(1) * τ(2)
```

Spin components are sorted by axis.

```@example canonical
S(2) * S(1)
```

Position stands to the left of momentum.

```@example canonical
p * x
```

[`normal_order`](@ref) applies the same rules explicitly, for code that expects a separate ordering step.

```@docs
normal_order
```

The ground-state projector is the one deliberate exception. It stays as it is, since rewriting it through completeness would multiply the number of terms in dissipators and equations of motion.

```@example canonical
σ(1, 2) * σ(2, 1)
```

Completeness is applied only on request.

```@example canonical
expand_completeness(σ(1, 2) * σ(2, 1))
```

```@docs
expand_completeness
```

An expression in canonical form has its own type, which is what products, sums and commutators of operators return.

```@docs
QAdd
```

Expressions also have a total order, used to sort reproducibly and to pick one of two expressions as the canonical one, such as an operator or its adjoint. It is a named comparison rather than `<`, since it says nothing about the operators as physical quantities.

```@example canonical
using SecondQuantizedAlgebra: canonical_isless
canonical_isless(σ(1, 2), σ(2, 1)), sort([σ(2, 1), σ(1, 2)]; lt = canonical_isless)
```

```@docs
canonical_isless
```

## Commutators

Since both orderings of a product are already canonical, a commutator costs no more than two products and a subtraction. The Pauli matrices show both brackets at once: different components anticommute, and their commutator is the third component.

```@example canonical
commutator(τ(1), τ(2)), anticommutator(τ(1), τ(2))
```

For a mode, the anticommutator of ``a`` and ``a^\dagger`` is the symmetric combination that enters vacuum fluctuations.

```@example canonical
anticommutator(a, a')
```

```@docs
commutator
anticommutator
```

## Finishing coefficients

Coefficients are exact numbers or Symbolics expressions, as described in [Parameters and coefficients](@ref). Polynomial coefficients cancel as terms combine, while identities between functions are applied only by [`simplify`](@ref simplify(::QSym)).

```@example kerr
@variables θ
simplify((sin(θ)^2 + cos(θ)^2) * a' * a)
```

[`expand`](@ref expand(::QAdd)) multiplies out products that sit inside functions of the coefficients, which are otherwise kept as written.

```@example kerr
@variables x y
expand(exp(x * (x + y)) * a)
```

```@docs
simplify(::QSym)
expand(::QAdd)
```

## Symmetric ordering

Phase-space methods such as the truncated Wigner approximation work with symmetrically ordered operators rather than normal-ordered ones. In symmetric order the photon number ``a^\dagger a`` reads as ``(a^\dagger a + a a^\dagger)/2 - 1/2``, so the half quantum of vacuum noise becomes explicit.

```@example kerr
n_sym = normal_to_symmetric(a' * a)
```

The conversion is exact in both directions.

```@example kerr
symmetric_to_normal(n_sym) == a' * a
```

```@docs
normal_to_symmetric
symmetric_to_normal
```
