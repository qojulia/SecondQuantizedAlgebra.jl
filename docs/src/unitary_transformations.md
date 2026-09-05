# Unitary transformations

A unitary transformation changes the operator basis while preserving the
operator algebra. In SecondQuantizedAlgebra, a [`UnitaryTransform`](@ref) stores
the transformed fundamental operators, their inverse transformation, and—when
the basis moves in time—the corresponding Hamiltonian gauge term.

The basic workflow is:

```julia
using SecondQuantizedAlgebra

h = FockSpace(:resonator)
@qnumbers a::Destroy(h)
@variables θ ω t

U = Rotation(a, θ)
conjugate(a, U)
```

Here [`conjugate(A, U)`](@ref conjugate) evaluates
``U^\dagger A U`` by substituting the transformed fundamental operators into
`A`. The same substitution works for products and sums:

```julia
n = a' * a
X = a + a'

conjugate(n, U)
conjugate(X^2, U)
```

## Available transformations

The named constructors select the transformation from the operator type and
their arguments.

| System | Constructor | Action |
|---|---|---|
| Fock mode | `Displace(a, α)` | ``a \mapsto a+\alpha`` |
| Fock mode | `Rotation(a, θ)` | ``a \mapsto e^{-i\theta}a`` |
| Fock mode | `Squeeze(a, r, ϕ=0)` | ``a \mapsto \cosh(r)a+e^{i\phi}\sinh(r)a^\dagger`` |
| Two Fock modes | `Rotation(a, b, θ)` | beamsplitter rotation |
| Two Fock modes | `Squeeze(a, b, r)` | two-mode squeezing |
| Canonical quadratures | `Displace(x, p, dx, dp)` | ``x\mapsto x+dx``, ``p\mapsto p+dp`` |
| Canonical quadratures | `Rotation(x, p, θ)` | phase-space rotation |
| Canonical quadratures | `Squeeze(x, p, r)` | ``x\mapsto e^r x``, ``p\mapsto e^{-r}p`` |
| Spin or Pauli operators | `Rotation(S, axis, θ)` | rotation around axis 1, 2, or 3 |
| N-level transitions | `Rotation(σ, W)` | basis change defined by the unitary matrix `W` |

Each constructor also defines the inverse transformation, so `inv(U)` can be
used without deriving another set of rules.

## Static and time-dependent transformations

Use [`conjugate`](@ref) for observables and static changes of basis. For a
Hamiltonian in a moving basis, use [`transform`](@ref SecondQuantizedAlgebra.transform):

```math
\operatorname{transform}(H,U)
= U^\dagger H U + i(\partial_t U^\dagger)U.
```

The second term accounts for the motion of the basis. It is available through
[`gauge_term(U)`](@ref gauge_term).

```julia
H = ω * a' * a

Ustatic = Rotation(a, θ)
transform(H, Ustatic) == conjugate(H, Ustatic)

Utimed = Rotation(a, ω * t, t)
Hmoving = transform(H, Utimed)
gauge_term(Utimed)
```

The final argument identifies the differentiation variable. It must be a real
symbolic variable. When a parameter depends on time, pass the time variable to
the constructor if the gauge term is needed; for example, use
`Rotation(a, ω*t, t)` for a rotating Hamiltonian frame.


## Inversion and composition

`inv(U)` reverses a transformation. Transform products follow application
order: the left transformation is applied first.

```julia
U1 = Displace(a, 1 // 3)
U2 = Rotation(a, θ)
U = U1 * U2

conjugate(a, U) == conjugate(conjugate(a, U1), U2)
iszero(simplify(conjugate(conjugate(a, U), inv(U)) - a))
```

Static and timed transformations can be composed. Timed transformations in a
single product must use the same time variable. The gauge terms are composed
in the same order as the operator maps, so `transform(H, U1 * U2)` agrees with
applying the two frame changes successively.

[`generators(U)`](@ref generators) returns the fundamental operators acted on
by a transformation.

## N-level basis rotations

For an ordinary [`NLevelSpace`](@ref), `Rotation(σ, W)` transforms all
transitions on the same site according to the basis matrix `W`.

```julia
h_atom = NLevelSpace(:atom, (:g, :e))
σ = Transition(h_atom, :σ, 1, 2)

W = [cos(θ) -sin(θ); sin(θ) cos(θ)]
Ulevels = Rotation(σ, W)
conjugate(σ, Ulevels)
```

`W` must be square, match the number of levels, and satisfy
``W^\dagger W=I`` symbolically. For a time-dependent matrix, use
`Rotation(σ, W, t)`. Its gauge is computed entrywise from
``i\dot W^\dagger W``.

## API reference

```@docs
conjugate
gauge_term
generators
BeamSplitter
TwoModeSqueeze
BasisRotation
Bogoliubov
RotatingFrame
DisplacementFrame
```
