```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Operator functions

A trapped ion driven by a laser picks up the phase of the light at its position.
Its motion enters through the displacement ``e^{i\eta(a + a^\dagger)}``, with
Lamb–Dicke parameter ``\eta``, giving the Hamiltonian

```math
H = \nu a^\dagger a + \frac{\Omega}{2} \left( \sigma^{21} e^{i\eta(a + a^\dagger)} + \sigma^{12} e^{-i\eta(a + a^\dagger)} \right).
```

We will construct this Hamiltonian, change its parameters and frame, then compare
a polynomial approximation with evaluation on a finite numerical basis.

```@contents
Pages = ["operator_functions.md"]
Depth = 2:3
```

## Construct the Hamiltonian

We write the laser phase with `expim` and pair the excitation term with its adjoint:

```@example functions
using SecondQuantizedAlgebra
h = FockSpace(:motion) ⊗ NLevelSpace(:ion, 2)
@qnumbers a::Destroy(h, 1)
σ(k, l) = Transition(h, :σ, k, l)
@variables η Ω ν

H = ν * a' * a + Ω / 2 * (σ(2, 1) * expim(η * (a + a')) + σ(1, 2) * expim(-η * (a + a')))
```

```@docs
expim(::QField)
sin(::QField)
cos(::QField)
QExpr
```

## Change parameters and frame

The two drive terms give a Hermitian Hamiltonian:

```@example functions
H' == H
```

Setting the Lamb–Dicke parameter to zero recovers the bare internal-state drive:

```@example functions
substitute(H, Dict(η => 0))
```

To see how the motional quadrature changes in a rotated frame, consider its cosine:

```@example functions
@variables θ
R = Rotation(a, θ)
conjugate(cos(a + a'), R) == cos(conjugate(a + a', R))
```

## Obtain a polynomial approximation

For small ``\eta``, the first-order drive contains the sidebands
``a \sigma^{21}`` and ``a^\dagger \sigma^{21}``, which couple the internal state to motion:

```@example functions
taylor(H, 0:1)
```

```@docs
taylor
```

## Evaluate on a finite basis

For a numerical calculation, choose a motional cutoff and parameter values, then
convert the original Hamiltonian with [`to_numeric`](@ref):

```@example functions
using QuantumOpticsBase
backend = QuantumOpticsBackend()
parameter = Dict(η => 0.1, Ω => 1.0, ν => 1.0)
Hn = to_numeric(H, h, (10, 2); backend, parameter)
size(Hn)
```

The same setup lets us compare the first-order approximation with that result:

```@example functions
using LinearAlgebra: norm
Hfirst = to_numeric(taylor(H, 0:1), h, (10, 2); backend, parameter)
norm(dense(Hn - Hfirst).data)
```

Varying the Taylor order examines the symbolic approximation; varying the motional
cutoff examines the numerical truncation. See [`to_numeric`](@ref) for conversion
forms, parameter requirements, and limitations.
