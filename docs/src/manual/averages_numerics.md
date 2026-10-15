```@meta
CurrentModule = SecondQuantizedAlgebra
```

# Averages and numerics

Measured quantities are expectation values, and the equations of motion of an open or closed system are usually written for them. The Heisenberg equation turns into ``\frac{d}{dt}\langle X \rangle = i \langle [H, X] \rangle``, a relation between averages of operator products. Consider a driven Jaynes–Cummings model,

```math
H = \Delta a^\dagger a + g \left(a^\dagger \sigma^{ge} + \sigma^{eg} a\right) + \eta \left(a + a^\dagger\right).
```

```@contents
Pages = ["averages_numerics.md"]
Depth = 2:3
```

## Averages

[`average`](@ref) maps an operator expression to its expectation value. Averaging is linear: it distributes over sums and pulls out the coefficients, so the right-hand side of the Heisenberg equation for ``a`` becomes a sum of averages.

```@example avg
using SecondQuantizedAlgebra
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
@variables Δ g η t

H = Δ * a' * a + g * (a' * σ(:g, :e) + σ(:e, :g) * a) + η * (a + a')
da = average(1im * commutator(H, a))
```

An average is an ordinary symbolic scalar, so it takes part in Symbolics arithmetic. [`undo_average`](@ref) recovers the operator expression.

```@example avg
undo_average(da)
```

```@docs
average
undo_average
is_average
```

### Conjugation

The equation for ``\langle a^\dagger \rangle`` is the complex conjugate of the one for ``\langle a \rangle``. Taken with `conj`, the conjugate of an average stays formal. [`inner_adjoint`](@ref) moves it inside, where it becomes the average of the adjoint operator, ``\langle X \rangle^* = \langle X^\dagger \rangle``, so the conjugate equation is again written in averages.

```@example avg
inner_adjoint(da)
```

```@docs
inner_adjoint
```

### Time dependence

To integrate the equations, every average becomes an unknown function of time.

```@example avg
make_time_dependent(da, t)
```

```@docs
make_time_dependent
```

## Operator bases

The equation for ``\langle a \rangle`` involves ``\langle \sigma^{ge} \rangle``, whose own equation brings in products such as ``\langle a \sigma^{ee} \rangle``. Closing this hierarchy starts from a set of operators to track. [`fundamental_operators`](@ref) gives the generators of the algebra on a space, and [`find_operators`](@ref) all products up to a given order, each counted once together with its adjoint.

```@example avg
fundamental_operators(h)
```

```@example avg
find_operators(h, 2)
```

[`unique_up_to_adjoint`](@ref) applies the same reduction to a list of your own.

```@example avg
unique_up_to_adjoint([a, a', a' * a, σ(:g, :e), σ(:e, :g)])
```

```@docs
fundamental_operators
find_operators
unique_up_to_adjoint
```

## Many emitters

With ``N`` atoms in the cavity, each with its own coupling ``g_i``, the interaction becomes a [`Σ`](@ref) over the atoms. The average of a sum is a sum of averages over the same index, with the same excluded sites, so the equations keep ``N`` symbolic.

```@example sums
using SecondQuantizedAlgebra
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
a = Destroy(h, :a, 1)
@variables N
i = Index(h, :i, N, 2)
σ(k, l) = IndexedOperator(Transition(h, :σ, k, l, 2), i)
g = IndexedVariable(:g, i)

H = Σ(g * (a' * σ(1, 2) + a * σ(2, 1)), i)
average(H)
```

## Numerical conversion

[`to_numeric`](@ref) turns an operator expression into a matrix of [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl) or [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl), ready for their solvers. The numerical side is an extension: it becomes available once one of the two packages is loaded, and the backend is chosen with [`QuantumOpticsBackend`](@ref) or [`QuantumToolboxBackend`](@ref). The portable form takes the symbolic Hilbert space and one truncation per factor, a photon cutoff for a mode. Every parameter needs a value.

```@example avg
using QuantumOpticsBase
backend = QuantumOpticsBackend()
parameter = Dict(Δ => 0.5, g => 0.1, η => 0.2)
Hn = to_numeric(H, h, (10, 2); backend, parameter)
size(Hn)
```

```@docs
to_numeric
QuantumOpticsBackend
QuantumToolboxBackend
```

### Bases and levels

A basis of the backend can be passed in place of the space and its truncation, and then also fixes the backend. Symbolic level names map to basis indices in the order of the levels of the [`NLevelSpace`](@ref), so ``g`` is the first basis state and ``e`` the second.

```@example avg
bf = FockBasis(10)
bn = NLevelBasis(2)
b = bf ⊗ bn
to_numeric(σ(:g, :e), b) == embed(b, 2, transition(bn, 1, 2))
```

### Time-dependent parameters

Parameters that change in time go into `time_parameter`, as functions of time. The result is then the backend's time-dependent operator, a `TimeDependentSum` for QuantumOptics.

```@example avg
Ht = to_numeric(H, h, (10, 2); backend, parameter = Dict(Δ => 0.5, g => 0.1), time_parameter = Dict(η => t -> 0.2 * sin(t)))
Ht isa TimeDependentSum
```

An explicit time dependence of the coefficients, such as a phase from a rotating frame, needs the time variable itself, given as `t => identity`. A constant may stand in `time_parameter` next to it.

```@example avg
Hd = η * (expim(Δ * t) * a + expim(-Δ * t) * a')
to_numeric(Hd, h, (10, 2); backend, parameter = Dict(Δ => 0.5), time_parameter = Dict(t => identity, η => 0.2)) isa TimeDependentSum
```

### Custom operators and lazy output

`operators` replaces the matrix of an operator by one of your choice, given on the full basis, and the adjoint is replaced along with it. Replacing ``a`` by ``a + \alpha`` displaces the mode numerically.

```@example avg
α = 0.5
a_displaced = embed(b, 1, destroy(bf)) + α * identityoperator(b)
to_numeric(a' * a, b; operators = Dict(a => a_displaced)) ≈ to_numeric((a' + α) * (a + α), b)
```

The default result is a sparse matrix. On a large product space, `op_type = identity` keeps the backend's lazy sum of tensor products and never forms the full matrix.

```@example avg
to_numeric(H, b; parameter, op_type = identity) isa LazySum
```

### Indexed systems

A sum over an index is written out over explicit copies of the subsystem, so the basis fixes the number of atoms. It has one factor per copy, and `sites` maps the position of each subsystem in the symbolic space to the positions of its copies. Parameters that depend on the site are keyed by site, obtained by substituting a number for the index.

```@example sums
using QuantumOpticsBase
b = FockBasis(3) ⊗ NLevelBasis(2) ⊗ NLevelBasis(2)
couplings = Dict(substitute(g, Dict(index_sym(i) => k)) => 0.1k for k in 1:2)
Hn = to_numeric(H, b, Dict(1 => [1], 2 => [2, 3]); parameter = couplings)
size(Hn)
```

```@docs
to_numeric(::QAdd, ::Basis, ::AbstractDict{Int, Vector{Int}})
```

## Expectation values

[`numeric_average`](@ref) evaluates an operator expression, or an expression of averages, in a state. For a coherent state with amplitude ``\alpha`` and the atom in its ground state, ``\langle \sigma^{ge} \rangle`` vanishes and the right-hand side of the equation for ``\langle a \rangle`` is ``-i(\Delta\alpha + \eta) = -0.45i``.

```@example avg
ψ = coherentstate(bf, α) ⊗ nlevelstate(bn, 1)
numeric_average(substitute(da, parameter), ψ)
```

A vector of states, such as the output of a time evolution, gives a vector of values. `expect` is the same for operator expressions, under the name used by QuantumOptics.jl and QuantumToolbox.jl.

```@example avg
ψs = [coherentstate(bf, β) ⊗ nlevelstate(bn, 1) for β in (0.1, 0.5)]
numeric_average(a' * a, ψs)
```

```@example avg
expect(a' * a, ψ)
```

On a large product space, a `LazyKet` holds the state of each factor separately and avoids the full state vector.

```@example avg
ψl = LazyKet(b, (coherentstate(bf, 0.5), nlevelstate(bn, 2)))
numeric_average(average(a' * a + σ(:e, :e)), ψl)
```

```@docs
numeric_average
expect(::QField, ::Any, ::AbstractDict{<:QSym})
```


## Symbolic vacuum expectations

`expect(op, Vacuum())` evaluates a vacuum expectation exactly, without a Hilbert-space cutoff or a numeric backend. Matrix elements and overlaps follow by placing the corresponding operator products between the vacuum states.

```@example vacuum
using SecondQuantizedAlgebra
@qnumbers v::Destroy(FockSpace(:vacuum_mode))
expect(v * v', Vacuum())
```

```@docs
Vacuum
expect(::QAdd, ::Vacuum)
kronecker_delta
is_kronecker_delta
```
