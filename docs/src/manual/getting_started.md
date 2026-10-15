```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Getting started

We use units with ``\hbar = 1``. Consider a cavity mode on resonance with a two-level atom, described by the Jaynes–Cummings Hamiltonian

```math
H = \omega a^\dagger a + \omega \sigma^{22} + g \left( a^\dagger \sigma^{12} + a \sigma^{21} \right).
```

The cavity lives on a [`FockSpace`](@ref) and the atom on an [`NLevelSpace`](@ref); their tensor product [`⊗`](@ref) is the space of the full system. Each space in the tensor product represents one subsystem. In the constructors below, the final argument selects that subsystem: `1` for the cavity and `2` for the atom. Thus the annihilation operator ``a`` acts on the cavity, while the transitions ``\sigma^{ij} = |i\rangle\langle j|`` act on the atom. The couplings are ordinary Symbolics variables.

```@example getting-started
using SecondQuantizedAlgebra

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)

@qnumbers a::Destroy(h, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
@variables ω g

H = ω * a' * a + ω * σ(2, 2) + g * (a' * σ(1, 2) + a * σ(2, 1))
```

Products apply the local algebraic identities as they are formed: ladder operators come out normal ordered, and transitions on the same atom compose.

```@example getting-started
a * a', σ(1, 2) * σ(2, 1), σ(2, 1) * σ(2, 1)
```

Here ``\sigma^{11} = |1\rangle\langle 1|`` projects onto the atomic ground state, and ``\sigma^{22} = |2\rangle\langle 2|`` onto the excited state. Since these are the only two states, their projectors satisfy ``\sigma^{11} + \sigma^{22} = 1``. The algebra keeps ``\sigma^{11}`` as a single operator; call [`expand_completeness`](@ref) to replace it using this identity.

```@example getting-started
expand_completeness(σ(1, 1))
```

Because every product is already in canonical form, a commutator reduces to a subtraction, and the Heisenberg equation ``\dot{a} = i[H, a]`` is one line.

```@example getting-started
1im * commutator(H, a)
```

Taking expectation values with [`average`](@ref) turns the operator equation for the atomic coherence into an equation for c-numbers, the starting point of any mean-field or cumulant treatment.

```@example getting-started
average(1im * commutator(H, σ(1, 2)))
```

To compute with the model, [`to_numeric`](@ref) builds its matrix through QuantumOptics on a truncated basis, here with cavity Fock states from zero to ten photons.

```@example getting-started
using QuantumOpticsBase
using LinearAlgebra: eigvals, Hermitian

backend = QuantumOpticsBackend()
parameter = Dict(ω => 1.0, g => 0.1)
Hn = to_numeric(H, h, (10, 2); backend, parameter)
eigvals(Hermitian(Matrix(dense(Hn).data)))[1:5]
```

Above the ground state the spectrum splits into doublets at ``\omega \pm g`` and ``2\omega \pm \sqrt{2} g``, the vacuum Rabi splitting.
