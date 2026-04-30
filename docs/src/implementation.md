```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Implementation

This page walks through the main concepts in **SecondQuantizedAlgebra.jl**: defining Hilbert spaces, creating operators, building algebraic expressions, and converting to numerics.

## Hilbert spaces

Every system starts with a Hilbert space. The package provides five concrete spaces:

| Space | Description |
|-------|-------------|
| [`FockSpace`](@ref) | Bosonic (quantum harmonic oscillator) |
| [`NLevelSpace`](@ref) | Finite discrete levels (atoms, qubits) |
| [`PauliSpace`](@ref) | Two-level Pauli systems |
| [`SpinSpace`](@ref) | Collective spin angular momentum |
| [`PhaseSpace`](@ref) | Quadrature (position/momentum) |

```@example hilbert-space
using SecondQuantizedAlgebra # hide
hf = FockSpace(:cavity)
nothing # hide
```

[`NLevelSpace`](@ref) requires the number of levels. Levels can be integers or symbolic names, and a ground state can be specified (defaults to the first level):

```@example hilbert-space
ha = NLevelSpace(:atom, 2)                  # levels 1, 2; ground state = 1
ha_sym = NLevelSpace(:atom, (:g, :e))       # symbolic levels; ground state = :g
ha_gs = NLevelSpace(:atom, 4, 2)            # 4 levels; ground state = 2
nothing # hide
```

The ground state projector is eliminated during simplification using the completeness relation ``\sum_j |j\rangle\langle j| = 1``.

Composite systems are built with [`⊗`](@ref) (`\otimes<tab>`) or [`tensor`](@ref):

```@example hilbert-space
h = hf ⊗ ha
h3 = hf ⊗ ha ⊗ NLevelSpace(:spin, 3)
nothing # hide
```


## Operators

Operators are the symbolic building blocks. Each operator type lives on a specific Hilbert space:

| Operator | Space | Description |
|----------|-------|-------------|
| [`Destroy`](@ref) / [`Create`](@ref) | [`FockSpace`](@ref) | Bosonic ladder operators (``a``, ``a^\dagger``) |
| [`Transition`](@ref) | [`NLevelSpace`](@ref) | Level transition ``\|i\rangle\langle j\|`` |
| [`Pauli`](@ref) | [`PauliSpace`](@ref) | Pauli matrices ``\sigma_x, \sigma_y, \sigma_z`` |
| [`Spin`](@ref) | [`SpinSpace`](@ref) | Angular momentum ``S_x, S_y, S_z`` |
| [`Position`](@ref) / [`Momentum`](@ref) | [`PhaseSpace`](@ref) | Quadrature operators ``x``, ``p`` |

All operators are subtypes of `QSym`. Here are some examples:

```@example operators
using SecondQuantizedAlgebra # hide
hf = FockSpace(:cavity)
a = Destroy(hf, :a)
a'   # Create — the adjoint of Destroy
nothing # hide
```

```@example operators
ha = NLevelSpace(:atom, (:g, :e))
σge = Transition(ha, :σ, :g, :e)   # |g⟩⟨e|
nothing # hide
```

```@example operators
hp = PauliSpace(:qubit)
σx = Pauli(hp, :σ, 1)   # axis: 1=x, 2=y, 3=z
nothing # hide
```

In composite systems, specify which subspace the operator acts on by index:

```@example operators
h = FockSpace(:c1) ⊗ FockSpace(:c2) ⊗ NLevelSpace(:atom, 2)
a = Destroy(h, :a, 1)             # acts on first FockSpace
b = Destroy(h, :b, 2)             # acts on second FockSpace
σ12 = Transition(h, :σ, 1, 2, 3)  # acts on the NLevelSpace
nothing # hide
```

For convenience, the [`@qnumbers`](@ref) macro creates operators in one line:

```@example operators
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
@qnumbers a::Destroy(h, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
nothing # hide
```


## Symbolic parameters

Symbolic parameters (c-numbers) are created using `@variables` from [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), which is re-exported by SecondQuantizedAlgebra. Variables are real by default. Annotate with `::Complex` to declare a complex parameter:

```@example c-numbers
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
@variables ω η::Complex     # ω is real, η is complex

H = ω * a' * a + η * (a + a')
nothing # hide
```

Real variables satisfy `conj(ω) == ω`, which simplifies adjoint expressions:


```@example c-numbers
conj(ω)
```

```@example c-numbers
conj(η)
```

## Algebraic expressions and commutation relations

All operator expressions are stored as [`QAdd`](@ref), a dictionary mapping operator sequences to prefactors. The global ordering convention controls whether commutation rules are applied at construction time (i.e. when operators are multiplied with `*`). Use [`set_ordering!`](@ref) and [`get_ordering`](@ref) to switch between them.

**[`NormalOrder`](@ref)** (default): Commutation rules are applied eagerly at construction. Every `*` call immediately normal-orders the result.

```@example ordering
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity) # hide
@qnumbers a::Destroy(h) # hide
set_ordering!(NormalOrder())
a * a'   # immediately gives a†a + 1
```

**[`LazyOrder`](@ref)**: Operators are concatenated as-is at construction — no reordering happens during `*`.

```@example ordering
set_ordering!(LazyOrder())
expr = a * a'               # keeps a·a† as-is
```

```@example ordering
set_ordering!(NormalOrder()) # hide
```

Regardless of the global ordering, [`normal_order`](@ref) can always be called explicitly to normal-order any expression:

```@example ordering
normal_order(expr)           # applies commutation rules: a†a + 1
```

Note that [`simplify`](@ref) applies ordering-independent algebraic identities (Transition composition, Pauli products) but never commutation-based reordering.

## Averaging

In many-body theory, equations of motion are typically written for expectation values rather than operators. The [`average`](@ref) function converts an operator expression into a symbolic scalar representing its expectation value ``\langle \hat{O} \rangle``:

```@example averaging
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity) # hide
@qnumbers a::Destroy(h) # hide
@variables ω # hide
average(a)
```

Averaging is linear — it distributes over sums and pulls out c-number prefactors:

```@example averaging
average(ω * a' * a + a + a')
```

The result is a symbolic expression (a `Num` from Symbolics.jl) that can participate in standard symbolic arithmetic. To recover the underlying operator expression, use [`undo_average`](@ref):

```@example averaging
avg = average(a' * a)
undo_average(avg)
```

Averaged expressions also preserve summation metadata when working with indexed operators, so that symbolic sums carry through the averaging process (see [Symbolic Sums and Indices](@ref)).


## Numeric conversion

Operators can be converted to numeric representations using [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl).

### Direct conversion

Use [`to_numeric`](@ref) to convert an operator expression to a numeric operator:

```@example numeric
using SecondQuantizedAlgebra, QuantumOpticsBase
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)

b = FockBasis(7)
a_num = to_numeric(a, b)
@assert a_num == destroy(b)
nothing # hide
```

### Numeric averages

Use [`numeric_average`](@ref) to compute expectation values for a given state:

```@example numeric
α = 0.3 + 0.1im
ψ = coherentstate(b, α)
numeric_average(a' * a, ψ)
```

### Level mapping for NLevelSpace

When using symbolic levels, provide a `level_map` to specify the mapping from symbolic names to basis indices:

```@example numeric-nlevel
using SecondQuantizedAlgebra, QuantumOpticsBase
h = NLevelSpace(:atom, (:g, :e))
bn = NLevelBasis(2)
s = Transition(h, :s, :g, :e)
level_map = Dict(:g => 1, :e => 2)
@assert to_numeric(s, bn; level_map = level_map) == transition(bn, 1, 2)
nothing # hide
```

### Composite systems

For product spaces, [`to_numeric`](@ref) returns `LazyTensor` operators, which efficiently handle large tensor products:

```@example numeric-composite
using SecondQuantizedAlgebra, QuantumOpticsBase
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:a, :b, :c))
h = hc ⊗ ha

@qnumbers a::Destroy(h, 1)
s(i, j) = Transition(h, :s, i, j, 2)

bf = FockBasis(10)
bn = NLevelBasis(3)
bc = bf ⊗ bn

level_map = Dict(:a => 3, :b => 2, :c => 1)
a_num = to_numeric(a, bc)
s_num = to_numeric(s(:a, :c), bc; level_map = level_map)
nothing # hide
```

For very large systems, use `LazyKet` to avoid materializing the full state vector:

```@example numeric-composite
ψ = LazyKet(bc, (coherentstate(bf, 0.3), (nlevelstate(bn, 1) + nlevelstate(bn, 3)) / sqrt(2)))
avg = average(a' * s(:a, :c))
numeric_average(avg, ψ; level_map = level_map)
```
