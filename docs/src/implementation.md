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

The ground-state projector ``|g\rangle\langle g|`` is eliminated via completeness ``\sum_j |j\rangle\langle j| = 1``. Same-site composition that produces one (e.g. ``\sigma^{12} \cdot \sigma^{21}``) triggers this eagerly inside `*`. User-constructed ground-state projectors (e.g. `Transition(h, :σ, 1, 1)` typed directly, never multiplied) stay atomic until they enter a `*` or until you call [`expand_completeness`](@ref) explicitly.

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

Symbolic parameters (c-numbers) are created using `@variables` from [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), which is re-exported by SecondQuantizedAlgebra. Variables are real by default. The `symtype` annotation controls how conjugation is handled:

| Declaration | `conj` / `adjoint` | Representation |
|-------------|--------------------|----------------|
| `@variables ω` (or `ω::Real`) | identity (`conj(ω) == ω`) | one real symbol |
| `@variables η::Number` | symbolic `conj(η)` | one symbol, complex-valued |
| `@variables ζ::Complex` | `conj(ζ) == real(ζ) - im*imag(ζ)` | split into `real`/`imag` parts |

```@example c-numbers
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
@variables ω η::Number      # ω is real, η is complex

H = ω * a' * a + η * a + conj(η) * a'
nothing # hide
```

Real variables satisfy `conj(ω) == ω`, which simplifies adjoint expressions:

```@example c-numbers
conj(ω)
```

A `::Number` parameter is complex-valued but stays a single symbol, so its conjugate is the symbolic `conj(η)`:

```@example c-numbers
conj(η)
```

Use `::Number` for a complex parameter you want kept atomic (e.g. a coupling amplitude); use `::Complex` when you want the parameter decomposed into independent real and imaginary unknowns. `::Number` keeps coefficient arithmetic on a single symbol, while `::Complex` carries both parts through every product.

## Algebraic expressions and commutation relations

All operator expressions are stored as [`QAdd`](@ref), a dictionary mapping operator sequences to prefactors. Commutation rules are applied **eagerly** at construction time: every `*` immediately normal-orders the result, applies algebraic identities (Transition composition, Pauli products), and expands ground-state completeness. The dict key always reflects the canonical form.

```@example ordering
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity) # hide
@qnumbers a::Destroy(h) # hide
a * a'   # immediately gives a†a + 1
```

[`normal_order`](@ref) can be called explicitly as an idempotent finalizer (useful in tests and for hand-constructed expressions that bypass `*`):

```@example ordering
normal_order(a * a')   # a†a + 1
```

[`simplify`](@ref) runs [`normal_order`](@ref) first, then walks the resulting terms applying `Symbolics.simplify` to each coefficient and dropping any summation indices no surviving term depends on. Since eager `*` already produces canonical form, calling `simplify` on the output of `*` is typically idempotent at the operator level, but it can still reduce symbolic coefficients (e.g. collapsing `sin²(ω) + cos²(ω)` to `1`).

### Substitution

Use [`substitute`](@ref) for both scalar parameters and operator leaves. Mixed
rule dictionaries are split by key: `QSym` keys replace operators, while all
other keys are applied to coefficients through SymbolicUtils.

Operator substitutions are simultaneous and single-pass. The replacement
expression is not searched again for operator keys, so transformations such as
`a => g*a + h*b` are well-defined even though the right-hand side contains `a`.
Missing adjoint rules are generated by default; pass `replace_adjoint=false` to
match only the keys supplied explicitly.

```@example substitution
using SecondQuantizedAlgebra # hide
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
@variables Δ g

substitute(Δ * a, Dict(a => g * a + 1, Δ => 2.0))
```

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

The result is a `BasicSymbolic` (from SymbolicUtils.jl) that participates in standard symbolic arithmetic. To recover the underlying operator expression, use [`undo_average`](@ref):

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

Scalar parameters and custom numeric operators can be supplied directly at the
numeric boundary:

```@example numeric
@variables Δ::Real
H = Δ * a' * a
H_num = to_numeric(H, b; parameter = Dict(Δ => 2.0))
@assert H_num == 2.0 * create(b) * destroy(b)
nothing # hide
```

For time-dependent parameters, pass numbers or functions in `time_parameter`.
The result is a callable `t -> op(t)`.

### Numeric averages

Use [`numeric_average`](@ref) to compute expectation values for a given state:

```@example numeric
α = 0.3 + 0.1im
ψ = coherentstate(b, α)
numeric_average(a' * a, ψ)
```

For symbolic scalar expressions such as `average(a)`, call [`numeric_average`](@ref)
directly. The `expect` alias is intentionally kept for operator expressions (`QField`)
only.

### Symbolic levels for NLevelSpace

Symbolic level names are resolved to integer basis indices at `Transition`
construction time, using the order of the `levels` tuple passed to
[`NLevelSpace`](@ref). The first level maps to basis index `1`, the second to
`2`, and so on:

```@example numeric-nlevel
using SecondQuantizedAlgebra, QuantumOpticsBase
h = NLevelSpace(:atom, (:g, :e))
bn = NLevelBasis(2)
s = Transition(h, :s, :g, :e)
@assert to_numeric(s, bn) == transition(bn, 1, 2)
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

a_num = to_numeric(a, bc)
s_num = to_numeric(s(:a, :c), bc)
nothing # hide
```

For very large systems, use `LazyKet` to avoid materializing the full state vector:

```@example numeric-composite
ψ = LazyKet(bc, (coherentstate(bf, 0.3), (nlevelstate(bn, 1) + nlevelstate(bn, 3)) / sqrt(2)))
avg = average(a' * s(:a, :c))
numeric_average(avg, ψ)
```
