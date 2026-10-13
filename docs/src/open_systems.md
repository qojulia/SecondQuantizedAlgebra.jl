```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Open Quantum Systems

Dissipative dynamics are first-class objects. A Lindblad generator

```math
\mathcal{L} = -i[H,\cdot] + \sum \gamma\, \mathcal{D}[J, K],
\qquad \mathcal{D}[J, K]\rho = J \rho K^\dagger - \tfrac12\{K^\dagger J, \rho\},
```

is built, transformed, summed over sites and converted to numerics with the same operations as an operator. Every map acts on operators by a call: `L(ρ)` is the Schrödinger generator applied to `ρ`, and `transpose(L)(A)` is the Heisenberg drift of `A`.

## Three layers

| Type | Represents |
|:--|:--|
| [`Lindbladian`](@ref) | `-i[H,·] + Σ γ 𝒟[J, K]` with its structure retained |
| [`SuperOperator`](@ref) | any polynomial map `Σ c · L ρ R` |
| [`SuperExpr`](@ref) | maps whose sides contain formal operator functions |

A `Lindbladian` prints as `𝒟[J]`, transforms jump by jump and converts to solver inputs. A `SuperOperator` is closed under composition, `transpose` and `adjoint`, which leave the GKSL form. The return type of an operation depends only on the argument types: composing two Lindbladians gives a `SuperOperator`, and anything involving a formal function gives a `SuperExpr`.

```@example open
using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: SuperOperator, sprepost, spre, spost, exponential_form
using Symbolics: @variables

hc = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
a = Destroy(hc, :a, 1)
σ(i, j) = Transition(hc, :σ, i, j, 2)
@variables Δ g κ γ

H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
L = hamiltonian_action(H) + κ * dissipator(a) + γ * dissipator(σ(1, 2))
```

`SuperOperator(L)` lowers the structure to explicit sandwiches, and `sprepost(A, B)`, `spre(A)` and `spost(B)` build `AρB`, `Aρ` and `ρB` directly:

```@example open
@assert SuperOperator(dissipator(a)) ==
    sprepost(a, a') - (1 // 2) * spre(a' * a) - (1 // 2) * spost(a' * a)
nothing # hide
```

`SuperOperator`, `spre`, `spost`, `sprepost` and `liouvillian` are public but not exported, because QuantumOpticsBase and QuantumToolbox export the same names. With a backend loaded, its `spre(a)` forwards to the symbolic one.

## Physical contracts

A `Lindbladian` is a form, not a physicality claim: `H` and the rates may be non-Hermitian or negative. The construction normalizes each channel so that equal generators compare equal:

- a constant in a jump moves into the Hamiltonian, since `𝒟[J + α] = 𝒟[J] - i[(i/2)(ᾱJ - αJ†), ·]`;
- a numeric factor or a frame phase in a jump moves into the rate, so `𝒟[2a] == 4𝒟[a]`;
- an identity term of `H` is dropped.

```@example open
@assert dissipator(2a) == 4 * dissipator(a)
dissipator(a + 3)
```

Every Lindbladian preserves the trace, `transpose(L)(1) == 0`. The identity map does not, so `one(Lindbladian)` throws, `L^0 == one(SuperOperator)`, and `c + L` is a `SuperOperator`.

## Equations of motion

The Heisenberg drift is `transpose(L)(A)`, with no separate helper. Averaging it gives the mean-field equation:

```@example open
transpose(L)(a' * a)
```

```@example open
average(transpose(L)(σ(2, 2)))
```

## Correlated and collective channels

`dissipator(J, K)` is the cross channel `𝒟[J, K]`, and `liouvillian(H, J; rates)` follows the QuantumOptics convention, where `rates[i, j]` multiplies `J[i] ρ J[j]'`.

Over indexed sites, `Σ(L, i)` sums a generator, and `dissipator(J, Γ, i => j)` builds the collective channel `Σᵢ Σⱼ Γ 𝒟[Jᵢ, Jⱼ]`. The diagonal `i = j` is split off, so the result prints as separate local and correlated parts:

```@example open
ha = NLevelSpace(:atom, 2)
@variables N
i = Index(ha, :i, N, ha)
j = Index(ha, :j, N, ha)
s(x, y, k) = IndexedOperator(Transition(ha, :s, x, y), k)
Lc = dissipator(s(1, 2, i), DoubleIndexedVariable(:Γ, i, j), i => j)
```

## GKSL coordinates

A `SuperOperator` that came out of an algebraic computation can be read back as a `Lindbladian`. The extraction is exact in the operator algebra: the two-sided sandwiches give the Kossakowski matrix in an ordered `JumpBasis`, and the one-sided remainder must be a Hamiltonian commutator. A constant in a jump comes back as a Hamiltonian change, as in the normalization above:

```@example open
S = SuperOperator(dissipator(a + 3))
Lindbladian(S)
```

The coordinates depend on the basis. In the basis `(a, a + a†a)`, independent dephasing and loss become a correlated Kossakowski matrix:

```@example open
Ld2 = κ * dissipator(a) + γ * dissipator(a' * a)
kossakowski(Ld2, JumpBasis(a, a + a' * a))
```

Basis operators are taken modulo the identity and with `expand_completeness` applied, so `σ¹¹` on a two-level atom becomes `-σ²²`. A map with a jump outside the basis throws an `ArgumentError`, and a map that is not of GKSL form throws a `GKSLFormError`. The extraction is exact, so write rates as rationals such as `1//10` or as symbolic parameters. Floating-point rates carry roundoff from building the map, and a float map is usually rejected; a positive `rtol` opts into treating small float leftovers as roundoff, at the price of dropping genuine terms below that size. Positivity of the rates is not checked.

## Index rules

A map shares one summation scope between its two sides. The rules that keep it well defined throw instead of guessing:

- composing two maps that bind the same index throws; rename one with `change_index`;
- a bound index of one operand that is free in the other throws (capture);
- a side that carries its own `Σ` throws, since `J ρ J†` with a summed `J` has a double sum;
- powers of a map with bound indices throw for the same reason.

## Frames

`transform(L, U)` moves the generator into the frame of a [`UnitaryTransform`](@ref): the Hamiltonian is transformed, and every jump is conjugated. A rotating frame removes the drive phase from both:

```@example open
@variables t::Real ω::Real ωd::Real F::Real
hf = FockSpace(:c)
b = Destroy(hf, :b)
Ld = hamiltonian_action(ω * b' * b + F * (b' * exp(-im * ωd * t) + b * exp(im * ωd * t))) +
    κ * dissipator(b)
exponential_form(transform(Ld, RotatingFrame(ωd * b' * b, t)))
```

`DisplacementFrame(a, L)` removes a coherent drive from a damped linear generator. The displacement is the fixed point of the exact drift `transpose(L)(a) = λa + c`, so thermal heating and loss both enter:

```@example open
Lth = hamiltonian_action(2 * b' * b + 3 * (b + b')) + (1 // 2) * dissipator(b) +
    (1 // 10) * dissipator(b')
U = DisplacementFrame(b, Lth)
transpose(transform(Lth, U))(b)
```

The generator may be any `SuperOperator` whose drift is affine in the mode, `d/dt (a, a') = M (a, a') + c(t)`, such as the Caldeira–Leggett master equation of quantum Brownian motion. Its damping acts on the momentum alone, so the drift of `a` mixes in `a'`. `DisplacementFrame(x, p, L)` works on a quadrature pair, and with a time argument each harmonic of the drive gets its exact bounded response:

```@example open
m = PhaseSpace(:m)
xq, pq = Position(m, :x), Momentum(m, :p)
@variables ω0::Real η::Real Dp::Real
ad(A) = spre(A) - spost(A)
anti(A) = spre(A) + spost(A)
Lqbm = SuperOperator(hamiltonian_action(pq^2 / 2 + ω0^2 * xq^2 / 2 - F * xq)) -
    Dp / 2 * ad(xq) ∘ ad(xq) - im * η / 2 * ad(xq) ∘ anti(pq)
conjugate(xq, DisplacementFrame(xq, pq, Lqbm))
```

## Formal operator functions

A Hamiltonian or a jump that contains `cos`, `sin` or `expim` gives a `Lindbladian{QExpr}`; `taylor(L, 0:n)` lowers it to a polynomial one.

```@example open
φ = Position(PhaseSpace(:transmon), :φ)
Lt = hamiltonian_action(-cos(φ)) + γ * dissipator(φ)
taylor(Lt, 0:2)
```

## Numerics

`to_numeric` converts any polynomial map to the backend's superoperator, with the same `parameter`, `time_parameter` and indexed `sites` arguments as for operators. For solvers that take a Hamiltonian and collapse operators, `to_numeric_master` diagonalizes each correlated rate block:

```julia
using QuantumOptics
import SecondQuantizedAlgebra: to_numeric_master
H, c_ops = to_numeric_master(L, hc, (10, nothing); parameter = Dict(Δ => 1, g => 0.1, κ => 0.3, γ => 0.05))
timeevolution.master(tspan, ρ0, H, c_ops)
```

A rate matrix that is not positive semidefinite, or a non-Hermitian `H`, throws, since a solver needs a physical generator. A time-dependent `H(t)` is checked at a few generic sample times, so it must be Hermitian at every `t` for the conversion to be meaningful; the check catches a broken pairing, not a coincidence at the samples. Time-dependent maps and time-dependent diagonal rates need QuantumToolbox; a correlated block with time-dependent rates is refused, and `to_numeric(L, …; time_parameter)` gives the full generator for `mesolve` instead.
