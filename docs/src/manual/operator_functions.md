```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Operator functions

Some Hamiltonians are not polynomials in the ladder operators. A trapped ion driven by a laser picks up the phase of the light at its position, and its motion enters through the displacement ``e^{i\eta(a + a^\dagger)}`` with the Lamb–Dicke parameter ``\eta``,

```math
H = \nu a^\dagger a + \frac{\Omega}{2} \left( \sigma^{21} e^{i\eta(a + a^\dagger)} + \sigma^{12} e^{-i\eta(a + a^\dagger)} \right).
```

Functions of operators such as this exponential, or the cosine potential of a Josephson junction, are kept as formal functions. They take part in the algebra without being expanded, and a Taylor expansion to a chosen order turns them into polynomials when needed.

```@contents
Pages = ["operator_functions.md"]
Depth = 2:3
```

## Formal functions

The exponential ``e^{iA}`` of a Hermitian operator is written [`expim`](@ref expim(::QField)), which also writes unit phases of numbers. With it the trapped-ion Hamiltonian reads as on paper.

```@example functions
using SecondQuantizedAlgebra
h = FockSpace(:motion) ⊗ NLevelSpace(:ion, 2)
@qnumbers a::Destroy(h, 1)
σ(k, l) = Transition(h, :σ, k, l)
@variables η Ω ν

H = ν * a' * a + Ω / 2 * (σ(2, 1) * expim(η * (a + a')) + σ(1, 2) * expim(-η * (a + a')))
```

Sine and cosine of operators are formal in the same way. There is no `exp` of an operator; ``e^{iA}`` is written `expim(A)`.

```@docs
expim(::QField)
sin(::QField)
cos(::QField)
```

A formal function keeps its place in a product. It commutes with operators on other subsystems but not with operators it shares a subsystem with.

```@example functions
a * cos(a + a') == cos(a + a') * a, σ(2, 1) * cos(a + a') == cos(a + a') * σ(2, 1)
```

No identities between functions are applied, not even ``\sin^2 + \cos^2 = 1`` or ``e^{iA} e^{-iA} = 1``. Where such an identity is needed, it is applied by hand.

```@example functions
expim(η * (a + a')) * expim(η * (a + a'))'
```

An expression that contains a formal function has its own type, which holds it as written rather than in canonical form.

```@docs
QExpr
```

## Exact operations

Operations that are exact without an expansion pass through the functions to their arguments. These are the adjoint, [`substitute`](@ref), [`change_index`](@ref), [`normal_order`](@ref), [`simplify`](@ref), [`expand`](@ref), [`expand_completeness`](@ref), [`assume_distinct_index`](@ref) and the [`conjugate`](@ref) of a unitary transformation. The Hamiltonian above is Hermitian as a whole.

```@example functions
H' == H
```

Substitution reaches inside the functions, so the Lamb–Dicke parameter can be set to zero to recover the bare two-level drive.

```@example functions
substitute(H, Dict(η => 0))
```

A unitary transformation acts on the argument of the function, so the result is exact and needs no Baker–Campbell–Hausdorff series.

```@example functions
@variables θ
R = Rotation(a, θ)
conjugate(cos(a + a'), R) == cos(conjugate(a + a', R))
```

A [`commutator`](@ref) with a formal function stays formal, since evaluating it would need the series.

## Taylor expansion

For small ``\eta`` the exponentials are expanded. To first order the drive acquires the sidebands ``a \sigma^{21}`` and ``a^\dagger \sigma^{21}``, which couple the internal state to the motion.

```@example functions
taylor(H, 0:1)
```

```@docs
taylor
```

The operator method of `taylor` takes an operator expression and an order range. It is distinct from the scalar method `taylor(f, x, 0:n)` of Symbolics, which expands in a variable `x`.

## Numerical conversion

Static formal functions can be converted directly on a chosen finite basis. [`to_numeric`](@ref)
materializes the complete argument of each function, then applies the corresponding matrix
function. Outer sums, ordered products, scalar prefactors, and nested functions are preserved.
There is no implicit symbolic Taylor expansion. All free parameters, including those in
arguments and outer prefactors, must be replaced by static values via `parameter`.

```@example functions
using QuantumOpticsBase
backend = QuantumOpticsBackend()
parameter = Dict(η => 0.1, Ω => 1.0, ν => 1.0)
Hn = to_numeric(H, h, (10, 2); backend, parameter)
size(Hn)
```

The finite basis is a numerical truncation; increasing a bosonic cutoff can change the
result. Use [`taylor`](@ref) when a symbolic polynomial approximation is desired.

Formal conversion rejects nonempty `time_parameter`, including time-dependent scalar
prefactors. [`average`](@ref), [`numeric_average`](@ref), `expect`, and formal superoperator
conversion still require explicit polynomial lowering. The positional operator-dictionary
`QExpr` overloads also retain their existing refusals; use the keyword `operators` form.

A symbolic [`Σ`](@ref) over a formal function and conversion between normal and symmetric
order with [`normal_to_symmetric`](@ref) and [`symmetric_to_normal`](@ref) remain unsupported;
both require polynomial lowering. Perturbative frame changes such as Schrieffer–Wolff or
van Vleck transformations are separate effective-Hamiltonian calculations.
