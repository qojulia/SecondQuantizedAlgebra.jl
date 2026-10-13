# Unitary transformations

A unitary transformation changes the operator basis while preserving the
operator algebra. In SecondQuantizedAlgebra, a [`UnitaryTransform`](@ref) stores
an exact affine action together with compiled forward rules and—when the basis
moves in time—the corresponding Hamiltonian gauge term. Inverse rules are
derived exactly from the affine action when `inv(U)` is requested.

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
| Two Fock modes | `Rotation(a, b, θ)` | passive beam-splitter mixing |
| Two Fock modes | `Squeeze(a, b, r)` | two-mode squeezing |
| Bosonic modes | `Bogoliubov(modes, S)` | general Nambu-linear canonical map |
| Canonical quadratures | `Displace(x, p, dx, dp)` | ``x\mapsto x+dx``, ``p\mapsto p+dp`` |
| Canonical quadratures | `Rotation(x, p, θ)` | phase-space rotation |
| Canonical quadratures | `Squeeze(x, p, r)` | ``x\mapsto e^r x``, ``p\mapsto e^{-r}p`` |
| Spin or Pauli operators | `Rotation(S, axis, θ)` | rotation around axis 1, 2, or 3 |
| N-level transitions | `Rotation(σ, W)` | basis change defined by the unitary matrix `W` |

Each constructor defines an affine action with an exact structural inverse, so
`inv(U)` does not require the caller to supply inverse rules.

## Raw bosonic Bogoliubov maps

[`Bogoliubov`](@ref) uses Nambu ordering

```math
\Xi=(a_1,\ldots,a_N,a_1^\dagger,\ldots,a_N^\dagger)^T,
\qquad \Xi' = S\Xi.
```

The public spellings are

```julia
Bogoliubov(modes, S)
Bogoliubov(modes, U, V)
```

where the block form means ``a' = Ua + Va^\dagger``.

`Bogoliubov` is an exact algebraic API, not a symbolic theorem prover. Its
function contract requires the supplied matrix to preserve adjoints and the
bosonic commutator form

```math
SJS^\dagger=J.
```

Satisfying these mathematical canonicality conditions is the caller's
responsibility for both numeric and symbolic inputs. The constructor validates
structural requirements such as the selected modes and matrix dimensions, but
it does not introduce hidden assumptions, maintain canonicality states, use
numerical tolerances, or project a matrix onto the canonical group.

The structured `Squeeze` and two-mode `Rotation`/`Squeeze` overloads satisfy
their canonicality conditions by construction and remain the preferred spelling
when they apply.

Scalar parameter substitution acts on the affine transformation itself and then
recompiles its execution rules:

```julia
@variables u::Number v::Number
B = Bogoliubov(a, [u v; conj(v) conj(u)]) # caller asserts canonicality
Bnum = substitute(B, Dict(u => 5 // 3, v => 4 // 3))
```

The same contract applies after substitution: the caller is responsible for
preserving any symbolic preconditions they supplied.

## Ladder bases at different frequencies

One oscillator with quadratures ``x`` and ``p`` has a ladder operator for
every reference frequency, ``a_\omega = (\omega x + i p)/\sqrt{2\omega}``
(unit mass, ``\hbar = 1``). A driven oscillator is often best written in the
basis of the drive frequency ``\omega_r`` rather than its bare frequency
``\omega_0``. The two bases are related by a single-mode squeeze,

```math
a_{\omega_0} = u\, a_{\omega_r} + v\, a_{\omega_r}^\dagger, \qquad
u = \frac{\omega_0 + \omega_r}{2\sqrt{\omega_0\omega_r}}, \quad
v = \frac{\omega_0 - \omega_r}{2\sqrt{\omega_0\omega_r}},
```

which is `Squeeze(a, log(ω0/ωr)/2)` with algebraic coefficients.
[`quadratures(a, ω)`](@ref quadratures) returns ``(x, p)`` in terms of a Fock
mode taken at frequency `ω`, and [`ladder_operator(x, p, ω)`](@ref ladder_operator)
is its inverse. Writing the Hamiltonian in ``x`` and ``p`` selects the basis
through the frequency passed to `quadratures`:

```julia
import SecondQuantizedAlgebra: quadratures, ladder_operator, FrequencyChange

@qnumbers b::Destroy(FockSpace(:b))
@variables ω0::Real ωr::Real

x, p = quadratures(b, ωr)
H = simplify(p * p / 2 + ω0^2 * x * x / 2)
# (ω0²+ωr²)/(2ωr) (b'b + 1/2) + (ω0²-ωr²)/(4ωr) (b² + b'²)
```

An expression already written in the ``\omega_0`` basis is moved to the
``\omega_r`` basis by the unitary transformation
[`FrequencyChange(a, ω0, ωr)`](@ref FrequencyChange), and
`ladder_operator(quadratures(a, ω0)..., ωr)` gives ``a_{\omega_r}`` expressed
through ``a_{\omega_0}``. Integer and rational frequencies give exact radicals.
Symbolic frequencies are taken to be positive; since ``\sqrt{\omega^2}`` is not
folded to ``\omega`` for an unsigned symbol, the coefficients are built from
``\sqrt{2\omega_0}`` and ``\sqrt{2\omega_r}`` so that their squares cancel.

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

Internally, exact affine transformations remain partitioned into homogeneous
algebra blocks. Overlapping blocks compose by the affine group law; disjoint
Fock, phase-space, spin, or N-level blocks remain separate. This avoids generic
symbolic matrix inversion: each block uses the inverse formula of its canonical
algebra.

Every `UnitaryTransform` carries this affine representation. Compiled forward
rules are execution data derived from it rather than an alternate semantic
representation; inverse rules are derived from the structural inverse on demand.

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

`W` must be square, match the number of levels, and be unitary. Unitarity is a
mathematical precondition supplied by the caller; the constructor checks the
transition family and matrix dimensions rather than attempting a symbolic proof.
For a time-dependent matrix, use `Rotation(σ, W, t)`. Its gauge is computed
entrywise from ``i\dot W^\dagger W`` under the same unitary-matrix contract.

The complete public API is listed under [Unitary Transformations](@ref "API: Unitary").
