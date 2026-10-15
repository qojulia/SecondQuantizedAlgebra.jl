```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Unitary transformations

A change of frame replaces every operator ``X`` by ``U^\dagger X U`` for a unitary ``U``. Rotating frames remove fast oscillations, displacements move a driven mode to its steady state, and squeezing or Bogoliubov transformations diagonalize quadratic Hamiltonians. For the transformations that keep the operator algebra linear, the image of every operator is known in closed form, so the new Hamiltonian follows exactly, without a Baker–Campbell–Hausdorff series. Consider a Kerr oscillator driven at frequency ``\omega_d``,

```math
H = \omega a^\dagger a + K a^{\dagger 2} a^2 + \eta e^{-i\omega_d t} a^\dagger + \eta^* e^{i\omega_d t} a.
```

```@contents
Pages = ["unitary_transformations.md"]
Depth = 2:3
```

## Rotating frames

The drive oscillates at ``\omega_d``. In the frame rotating at the drive frequency, `Rotation(a, ω_d * t, t)`, the drive becomes static and the mode frequency is replaced by the detuning ``\omega - \omega_d``.

```@example frames
using SecondQuantizedAlgebra
a = Destroy(FockSpace(:cavity), :a)
@variables ω ω_d K t η::Number

H = ω * a' * a + K * a' * a' * a * a + η * expim(-ω_d * t) * a' + conj(η) * expim(ω_d * t) * a
H_rot = transform(H, Rotation(a, ω_d * t, t))
```

An observable is transformed with [`conjugate`](@ref), which gives ``U^\dagger X U``. A Hamiltonian in a frame that moves in time also picks up the term ``i(\partial_t U^\dagger)U``, so it is transformed with [`transform`](@ref). Here that term subtracts ``\omega_d a^\dagger a``.

```@example frames
gauge_term(Rotation(a, ω_d * t, t))
```

```@docs
conjugate
transform
gauge_term
```

[`RotatingFrame`](@ref) builds the same frame from the Hamiltonian it rotates with, which is convenient when that Hamiltonian involves several subsystems.

```@example frames
transform(H, RotatingFrame(ω_d * a' * a, t)) == H_rot
```

```@docs
RotatingFrame
Rotation
```

[`Rotation`](@ref) is one of several named constructors, which together cover the standard frames of quantum optics:

| System | Constructor | Action |
| --- | --- | --- |
| Mode | `Displace(a, α)` | ``a \mapsto a + \alpha`` |
| Mode | `Rotation(a, θ)` | ``a \mapsto e^{-i\theta} a`` |
| Mode | `Squeeze(a, r, ϕ)` | ``a \mapsto \cosh r\, a + e^{i\phi} \sinh r\, a^\dagger`` |
| Two modes | `Rotation(a, b, θ)` | beam splitter |
| Two modes | `Squeeze(a, b, r)` | two-mode squeezing |
| Several modes | `Bogoliubov(modes, S)` | general linear mixing of ``a`` and ``a^\dagger`` |
| Quadratures | `Displace(x, p, dx, dp)` | ``x \mapsto x + dx``, ``p \mapsto p + dp`` |
| Quadratures | `Rotation(x, p, θ)` | rotation in phase space |
| Quadratures | `Squeeze(x, p, r)` | ``x \mapsto e^{r} x``, ``p \mapsto e^{-r} p`` |
| Spin or Pauli | `Rotation(S, axis, θ)` | rotation about an axis |
| Finite-level system | `Rotation(σ, W)` | change of basis by a unitary matrix ``W`` |
| Any | `RotatingFrame(H0, t)` | interaction picture of ``H_0`` |
| Driven mode | `DisplacementFrame(a, H)` | displacement to the steady state |

`Displace`, `Rotation`, `Squeeze` and `DisplacementFrame` take a final time variable `t` when their parameters depend on time. The time variable is a real symbolic variable, and without it no time-derivative term is produced, so `transform` then coincides with `conjugate`.

## Displacements

A strong drive pushes the mode to a large coherent amplitude ``\alpha``, and the fluctuations around it are described in a displaced frame. Written with the detuning ``\Delta = \omega - \omega_d``, the displaced Hamiltonian shows how the Kerr term acts on the fluctuations: a frequency shift ``4K|\alpha|^2``, squeezing terms ``K\alpha^2 a^{\dagger 2}``, and a linear term that vanishes when ``\alpha`` solves the classical steady-state equation ``\Delta\alpha + 2K|\alpha|^2\alpha + \eta = 0``.

```@example frames
@variables Δ α::Number
H_det = substitute(H_rot, Dict(ω => ω_d + Δ))
transform(H_det, Displace(a, α))
```

[`DisplacementFrame`](@ref) finds the displacement from the linear part of a Hamiltonian, here ``\alpha = -\eta/\Delta``. It removes the drive exactly when ``K = 0``, and leaves the Kerr correction ``2K|\alpha|^2\alpha`` as a residual drive.

```@example frames
D = DisplacementFrame(a, Δ * a' * a + η * a' + conj(η) * a)
conjugate(a, D)
```

```@docs
DisplacementFrame
Displace
```

## Squeezing and Bogoliubov transformations

Around the steady state, the quadratic part of the displaced Hamiltonian has the form ``\delta a^\dagger a + \frac{\kappa}{2}(a^{\dagger 2} + a^2)``, with ``\delta = \Delta + 4K|\alpha|^2`` and ``\kappa = 2K\alpha^2`` for real ``\alpha``. A transformation that mixes annihilation and creation operators diagonalizes it. Under [`Squeeze`](@ref), the coefficient of ``a^2`` becomes ``\frac{1}{2}(\kappa\cosh 2r + \delta\sinh 2r)``, which vanishes for ``\tanh 2r = -\kappa/\delta``.

```@example frames
@variables δ κ r
conjugate(δ * a' * a + κ / 2 * (a' * a' + a * a), Squeeze(a, r))
```

[`Bogoliubov`](@ref) does the same for any number of modes with a given matrix.

```@docs
Squeeze
Bogoliubov
```

## Composition and inversion

Transformations compose with `*`, the left one applied first, and are inverted with `inv`, which needs no inverse rules from the caller. The time-derivative terms compose in the same order as the operator maps, so a product of frames gives the same Hamiltonian as the two frame changes applied one after the other. The rotating and displaced frames combine into one.

```@example frames
U = Rotation(a, ω_d * t, t) * D
transform(substitute(H, Dict(ω => ω_d + Δ)), U) == transform(H_det, D)
```

Applying a transformation and then its inverse gives back the original expression.

```@example frames
X = a' * a + K * a' * a' * a * a
conjugate(conjugate(X, D), inv(D)) == X
```

Parameters of a transformation are set with [`substitute`](@ref), so a symbolic frame is built once and evaluated later. Time-dependent transformations in one product must share the same time variable.

```@example frames
conjugate(a, substitute(Squeeze(a, r), Dict(r => 0))) == a
```

[`generators`](@ref) lists the operators a transformation acts on.

```@example frames
generators(U)
```

```@docs
generators
```

## Spins and finite-level systems

A spin or a set of Pauli matrices is rotated about one of its axes. A rotation about the ``y`` axis mixes ``S_z`` with ``S_x``.

```@example spins
using SecondQuantizedAlgebra
Sz = Spin(SpinSpace(:S), :S, 3)
@variables θ
conjugate(Sz, Rotation(Sz, 2, θ))
```

The levels of a finite-level system change to a new basis by a unitary matrix. A real rotation by ``\theta`` mixes the two levels of an atom.

```@example spins
h = NLevelSpace(:atom, (:g, :e))
σ = Transition(h, :σ, :g, :e)

W = [cos(θ) -sin(θ); sin(θ) cos(θ)]
conjugate(σ, Rotation(σ, W))
```

The matrix must be unitary, which is not checked. With a time-dependent matrix, `Rotation(σ, W, t)` carries the gauge term ``i \dot W^\dagger W``.

## General generators

A transformation can also be given by its Hermitian generator ``G`` as ``U = e^{-i\theta G}``. The image of each operator is found by closing the commutators with ``G`` and exponentiating the resulting matrix exactly. The tunnelling term ``a^\dagger b + b^\dagger a`` between two modes generates a beam splitter with a phase.

```@example generators
using SecondQuantizedAlgebra
h = FockSpace(:a) ⊗ FockSpace(:b)
a = Destroy(h, :a, 1)
b = Destroy(h, :b, 2)
@variables θ

conjugate(a, UnitaryTransform(a' * b + b' * a, θ))
```

The exponential is exact only for numeric rates, so a generator such as ``g(a^\dagger b + b^\dagger a)`` with a symbolic ``g`` throws an `ArgumentError`; it is written as `UnitaryTransform(a' * b + b' * a, g * θ)`, or with a named constructor.

```@docs
UnitaryTransform
```

Formal operator functions are transformed through their arguments, so the result stays exact.

```@example generators
conjugate(cos(a + a'), Rotation(a, θ))
```


## Ladder bases at different frequencies

One oscillator with quadratures ``x`` and ``p`` has a ladder operator for
every reference frequency, ``a_\omega = (\omega x + i p)/\sqrt{2\omega}``
(unit mass, ``\hbar = 1``). A driven oscillator is often best written in the
basis of the drive frequency ``\omega_r`` rather than its bare frequency
``\omega_0``. The two bases are related by a single-mode squeeze,

```math
a_{\omega_0} = u\, a_{\omega_r} + v\, a_{\omega_r}^\dagger, \qquad
u = \frac{\omega_0 + \omega_r}{2\sqrt{\omega_0\omega_r}}, \quad
v = \frac{\omega_0 - \omega_r}{2\sqrt{\omega_0\omega_r}},
```

which is `Squeeze(a, log(ω0/ωr)/2)` with algebraic coefficients.
[`quadratures(a, ω)`](@ref quadratures) returns ``(x, p)`` in terms of a Fock
mode taken at frequency `ω`, and [`ladder_operator(x, p, ω)`](@ref ladder_operator)
is its inverse. Writing the Hamiltonian in ``x`` and ``p`` selects the basis
through the frequency passed to `quadratures`:

```julia
import SecondQuantizedAlgebra: quadratures, ladder_operator, FrequencyChange

@qnumbers b::Destroy(FockSpace(:b))
@variables ω0::Real ωr::Real

x, p = quadratures(b, ωr)
H = simplify(p * p / 2 + ω0^2 * x * x / 2)
# (ω0²+ωr²)/(2ωr) (b'b + 1/2) + (ω0²-ωr²)/(4ωr) (b² + b'²)
```

An expression already written in the ``\omega_0`` basis is moved to the
``\omega_r`` basis by the unitary transformation
[`FrequencyChange(a, ω0, ωr)`](@ref FrequencyChange), and
`ladder_operator(quadratures(a, ω0)..., ωr)` gives ``a_{\omega_r}`` expressed
through ``a_{\omega_0}``. Integer and rational frequencies give exact radicals.
Symbolic frequencies are taken to be positive; since ``\sqrt{\omega^2}`` is not
folded to ``\omega`` for an unsigned symbol, the coefficients are built from
``\sqrt{2\omega_0}`` and ``\sqrt{2\omega_r}`` so that their squares cancel.



```@docs
FrequencyChange
quadratures
ladder_operator
```
