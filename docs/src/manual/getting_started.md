```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Getting started

Consider a cavity mode on resonance with a two-level atom, described by the Jaynes–Cummings Hamiltonian

```math
H = \omega a^\dagger a + \omega \sigma^{22} + g \left( a^\dagger \sigma^{12} + a \sigma^{21} \right).
```

The cavity lives on a [`FockSpace`](@ref) and the atom on an [`NLevelSpace`](@ref); their tensor product [`⊗`](@ref) is the space of the full system. Operators refer to a factor of that product by its position, so the annihilation operator ``a`` acts on the first and the transitions ``\sigma^{ij} = |i\rangle\langle j|`` on the second. The couplings are ordinary Symbolics variables.

```@example getting-started
using SecondQuantizedAlgebra

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)

@qnumbers a::Destroy(h, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
@variables ω g

H = ω * a' * a + ω * σ(2, 2) + g * (a' * σ(1, 2) + a * σ(2, 1))
```

Products are simplified as they are formed. Ladder operators come out normal ordered and transitions on the same atom compose, so there is never an unsimplified expression to clean up afterwards.

```@example getting-started
a * a', σ(1, 2) * σ(2, 1), σ(2, 1) * σ(2, 1)
```

The ground-state projector ``\sigma^{11}`` is kept as is rather than rewritten as ``1 - \sigma^{22}``; [`expand_completeness`](@ref) does that on request. Because every product is already in canonical form, a commutator reduces to a subtraction, and the Heisenberg equation ``\dot{a} = i[H, a]`` is one line.

```@example getting-started
1im * commutator(H, a)
```

Taking expectation values with [`average`](@ref) turns the operator equation for the atomic coherence into an equation for c-numbers, the starting point of any mean-field or cumulant treatment.

```@example getting-started
average(1im * commutator(H, σ(1, 2)))
```

To compute with the model, [`to_numeric`](@ref) builds its matrix on a truncated basis, here with ten photons and through QuantumOpticsBase.

```@example getting-started
using QuantumOpticsBase
using LinearAlgebra: eigvals, Hermitian

backend = QuantumOpticsBackend()
parameter = Dict(ω => 1.0, g => 0.1)
Hn = to_numeric(H, h, (10, 2); backend, parameter)
eigvals(Hermitian(Matrix(dense(Hn).data)))[1:5]
```

Above the ground state the spectrum splits into doublets at ``\omega \pm g`` and ``2\omega \pm \sqrt{2} g``, the vacuum Rabi splitting.
