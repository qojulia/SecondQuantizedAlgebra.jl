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

A summed average keeps its index, its excluded sites and its summand, and each can be read back.

```@example sums
using SecondQuantizedAlgebra: is_indexed_sum, get_sum_indices, get_sum_body
x = average(Σ(a' * σ(1, 2), i))
is_indexed_sum(x), get_sum_indices(x), get_sum_body(x)
```

```@docs
is_indexed_sum
has_sum_metadata
get_sum_indices(::SecondQuantizedAlgebra.SymbolicUtils.BasicSymbolic)
get_sum_non_equal
get_sum_body
indexed_sum
```

## Numerical conversion

[`to_numeric`](@ref) turns an operator expression into a matrix of [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl) or [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl), ready for their solvers. The numerical side is an extension: it becomes available once one of the two packages is loaded, and the backend is chosen with [`QuantumOpticsBackend`](@ref) or [`QuantumToolboxBackend`](@ref). Both backends accept `to_numeric(op, h; backend)` for N-level and Pauli spaces and their tensor products. When a truncation is needed, use `to_numeric(op, h, dims; backend)` with one entry per factor, such as a photon cutoff for a mode. Parameter values can be supplied at conversion time or, with QuantumToolbox, through solver parameters.

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

For an `NLevelSpace`, a `PauliSpace`, or a product of these, dimensions can be omitted
with either `QuantumOpticsBackend()` or `QuantumToolboxBackend()`. The N-level state
count comes from the space, and a Pauli space has two states.

```julia
import QuantumOpticsBase
import QuantumToolbox

atom = NLevelSpace(:atom, 2)
P = Transition(atom, :P, 2, 2)
Poptics = to_numeric(P, atom; backend = QuantumOpticsBackend())
Ptoolbox = to_numeric(P, atom; backend = QuantumToolboxBackend())
# Equivalent explicit-dimension forms:
Poptics = to_numeric(P, atom, 2; backend = QuantumOpticsBackend())
Ptoolbox = to_numeric(P, atom, 2; backend = QuantumToolboxBackend())

hfinite = atom ⊗ PauliSpace(:qubit)
Pfinite = Transition(hfinite, :P, 2, 2, 1)
for be in (QuantumOpticsBackend(), QuantumToolboxBackend())
    @assert to_numeric(Pfinite, hfinite; backend = be) ≈
        to_numeric(Pfinite, hfinite, (2, 2); backend = be)
end

@variables λ
Hp = to_numeric(λ * cos(P), atom;
    backend = QuantumToolboxBackend(), default_time_dependence = false)
Hp(Dict(λ => 0.5), 0.7)
```

Both Hilbert-space forms accept `default_time_dependence`. This keyword affects only
QuantumToolbox: the default `true` uses `p[λ](t)`, and `false` uses `p[λ]`.
QuantumOptics continues to use explicit `parameter` and `time_parameter` mappings;
`default_time_dependence` has no effect on that backend. The supplied Hilbert space
retains every tensor factor, including spectator subsystems absent from the expression.
Fock and phase spaces still need a cutoff, and spin spaces need a spin length.

A basis of the backend can be passed in place of the space and its truncation, and then also fixes the backend. Symbolic level names map to basis indices in the order of the levels of the [`NLevelSpace`](@ref), so ``g`` is the first basis state and ``e`` the second.

```@example avg
bf = FockBasis(10)
bn = NLevelBasis(2)
b = bf ⊗ bn
to_numeric(σ(:g, :e), b) == embed(b, 2, transition(bn, 1, 2))
```

QuantumToolbox also accepts `to_numeric(expr, QuantumToolboxBackend())` without explicit dimensions for unambiguous systems of N-level or Pauli operators. Level counts come from the operators, and bound index ranges determine the number of copies of each subsystem. The ranges must be positive integers known at conversion time; symbolic counts can be fixed through `parameter`.

```julia
import QuantumToolbox

atom = NLevelSpace(:atom, 2)
P = Transition(atom, :P, 2, 2)
i = Index(atom, :i, 3, atom)
S = Σ(IndexedOperator(P, i), i)
Sn = to_numeric(cos(S), QuantumToolboxBackend())
size(Sn) # (8, 8): three two-level sites
```

The QuantumToolbox backend-only shorthand infers only the subsystems represented in the expression. Supply a Hilbert space to retain finite spectator subsystems, or an explicit basis or dimensions to resolve an ambiguous layout, and an explicit `sites` mapping when needed. Bosonic and phase-space operators still require a cutoff, and spin operators require a spin length, through the existing space-and-dimensions form.

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

### QuantumToolbox solver parameters

With QuantumToolbox, scalar parameters left unspecified by `parameter` or `time_parameter` automatically become symbolic keys in the solver's parameter dictionary `p`, and the result is a native `QobjEvo`. The default `default_time_dependence = true` generates callbacks of the form `p[η](t)`; set it to `false` to use numeric values through `p[η]` instead. Explicit parameter mappings take precedence.

```julia
import QuantumToolbox

Ht = to_numeric(H, h, (10, 2); backend = QuantumToolboxBackend(), parameter = Dict(Δ => 0.5, g => 0.1))
Ht(Dict(η => t -> 0.2 * sin(t)), 0.5)

Hp = to_numeric(H, h, (10, 2); backend = QuantumToolboxBackend(), default_time_dependence = false)
Hp(Dict(Δ => 0.5, g => 0.1, η => 0.2), 0.5)
```

The same dictionary can be passed to a QuantumToolbox solver through its `params` keyword. Scalar prefactors of static matrix functions, such as `x * cos(S)`, also support these callbacks; parameters inside matrix-function arguments still require static values via `parameter`.

Conversion uses QuantumToolbox's standard `QobjEvo` constructors and SciMLOperators representation. Coefficients are ordinary `(p, t)` callbacks, with the same evaluation and solver behavior as a hand-built `QobjEvo`. Evaluating `Ht(p, t)` materializes an operator and allocates its result.

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
expect(::QField, ::Any)
```
