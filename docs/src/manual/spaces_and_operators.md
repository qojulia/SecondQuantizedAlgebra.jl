```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Hilbert spaces and operators

A model starts from the Hilbert space of each subsystem. The type of the space fixes which operators live on it and which algebraic relations they obey:

| Subsystem | Space | Operators | Defining relation |
| --- | --- | --- | --- |
| [Bosonic mode](@ref "Bosonic modes") | [`FockSpace`](@ref) | [`Destroy`](@ref), [`Create`](@ref) | ``[a, a^\dagger] = 1`` |
| [Finite-level system](@ref "Finite-level systems") | [`NLevelSpace`](@ref) | [`Transition`](@ref) | ``\sigma^{ij}\sigma^{kl} = \delta_{jk}\sigma^{il}`` |
| [Qubit in Pauli form](@ref "Pauli and spin operators") | [`PauliSpace`](@ref) | [`Pauli`](@ref) | ``\sigma_j\sigma_k = \delta_{jk} + i\epsilon_{jkl}\sigma_l`` |
| [Angular momentum](@ref "Pauli and spin operators") | [`SpinSpace`](@ref) | [`Spin`](@ref) | ``[S_j, S_k] = i\epsilon_{jkl}S_l`` |
| [Canonical quadratures](@ref "Quadratures") | [`PhaseSpace`](@ref) | [`Position`](@ref), [`Momentum`](@ref) | ``[x, p] = i`` |

We build a hybrid quantum system step by step: a cavity mode coupled to an atom with a ground and an excited state, a collective spin, and a mechanical mirror described by its position and momentum.

```@contents
Pages = ["spaces_and_operators.md"]
Depth = 2:3
```

## Hilbert spaces

Each subsystem gets its own space, labelled by a name. The atom lists its levels by name, which then serve to address its transitions.

```@example hybrid
using SecondQuantizedAlgebra
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g, :e))
hs = SpinSpace(:spin)
hm = PhaseSpace(:mirror)
hc, ha, hs, hm
```

All of them share the supertype [`HilbertSpace`](@ref).

```@docs
HilbertSpace
FockSpace
NLevelSpace
PauliSpace
SpinSpace
PhaseSpace
```

### Composite systems

Subsystems are joined with the tensor product. Its factors are numbered from left to right, and this number is how an operator refers to the subsystem it acts on: here the cavity is factor 1 and the mirror factor 4.

```@example hybrid
h = hc ⊗ ha ⊗ hs ⊗ hm
```

```@docs
⊗(::HilbertSpace, ::HilbertSpace)
tensor
ProductSpace
```

## Operators

An operator is built from its space, a name, and whatever else fixes it: two levels for a transition, an axis for Pauli and spin operators. The name belongs to the identity of the operator, so two modes with different names are independent. Writing every name twice, once for the variable and once for the operator, is avoided with [`@qnumbers`](@ref).

```@example hybrid
@qnumbers a::Destroy(h, 1) x::Position(h, 4) p::Momentum(h, 4)
a, x, p
```

```@docs
@qnumbers
```

### Bosonic modes

A mode is described by its annihilation operator. The creation operator follows as its adjoint, so it rarely needs to be constructed on its own.

```@example hybrid
a' == Create(h, :a, 1)
```

```@docs
Destroy
Create
```

### Finite-level systems

Atoms, qudits and other discrete systems are written in terms of transitions between their levels, addressed by the level names. Projectors are the special case of equal levels.

```@example hybrid
σ(i, j) = Transition(h, :σ, i, j, 2)
σ(:g, :e), σ(:e, :e)
```

```@docs
Transition
```

A transition stores its levels as numbers, in the order the space names them, together with the level count and the ground level.

```@example hybrid
using SecondQuantizedAlgebra: n_levels, ground_state, transition_levels
transition_levels(σ(:g, :e)), n_levels(σ(:g, :e)), ground_state(σ(:g, :e))
```

```@docs
transition_levels
n_levels
ground_state
```

### Pauli and spin operators

A two-level system can equally be written with Pauli matrices, which reduce every product to a single matrix or a number.

```@example hybrid
τ(k) = Pauli(PauliSpace(:qubit), :τ, k)
τ(1) * τ(1), τ(1) * τ(2)
```

Spin operators extend this to angular momenta of unspecified size, such as the collective spin of the model.

```@example hybrid
S(k) = Spin(h, :S, k, 3)
S(1), S(2), S(3)
```

```@docs
Pauli
Spin
```

The axis of a Pauli or spin operator is read back as `1`, `2` or `3`.

```@example hybrid
using SecondQuantizedAlgebra: operator_axis
operator_axis(τ(2)), operator_axis(S(3))
```

```@docs
operator_axis
```

### Quadratures

A mode can also be described by its position and momentum quadratures instead of ladder operators, as the mirror is here.

```@example hybrid
p * x
```

```@docs
Position
Momentum
```

## Inspecting operators

Every operator has the same type, so its kind is read with a function rather than a type check. Methods for any operator expression are written for [`QField`](@ref), and for single operators for [`QSym`](@ref).

```@example hybrid
optype(a'), is_create(a'), is_transition(σ(:g, :e))
```

```@docs
QField
QSym
Op
optype
OpKind
is_destroy
is_create
is_transition
is_collective_transition
is_pauli
is_spin
is_position
is_momentum
```

The name, the site index and the subsystems of an operator or an expression are read off directly. The coupling between the cavity and the atom acts on factors 1 and 2.

```@example hybrid
operator_name(σ(:g, :e)), acts_on(a' * σ(:g, :e))
```

```@docs
operator_name
operator_index
acts_on
set_acts_on
```

### Adjoints

The adjoint `'` of an expression conjugates its coefficients along with the operators. [`dagger`](@ref dagger(::QField)) is the same under the name used by QuantumOptics.jl.

```@example hybrid
@variables g::Number
dagger(g * a' * σ(:g, :e))
```

[`qadjoint`](@ref) applies to coefficients as well. Unlike `conj`, it takes the conjugate through to the individual parameters.

```@example hybrid
@variables θ
qadjoint(g^2 * cos(θ)), conj(g^2 * cos(θ))
```

```@docs
dagger(::QField)
qadjoint
```

## Printing and LaTeX

Expressions print with Unicode subscripts in the REPL. With Latexify.jl loaded, `latexify` turns any expression into LaTeX.

```@example latex
using SecondQuantizedAlgebra, Latexify
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
@qnumbers a::Destroy(h, 1) σ::Transition(h, 2, 1, 2)
@variables g

latexify(g * (a' * σ + σ' * a))
```

```@docs
transition_superscript
```
