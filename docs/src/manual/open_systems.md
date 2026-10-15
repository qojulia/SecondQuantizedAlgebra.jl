```@meta
CollapsedDocStrings = true
CurrentModule = SecondQuantizedAlgebra
```

# Open quantum systems

A cavity that leaks photons or an atom that emits into free space does not evolve unitarily. For a memoryless environment the density matrix obeys a master equation ``\dot\rho = \mathcal{L}\rho``, whose generator has the Gorini–Kossakowski–Sudarshan–Lindblad form

```math
\mathcal{L}\rho = -i[H, \rho] + \sum \gamma\, \mathcal{D}[J, K]\rho,
\qquad
\mathcal{D}[J, K]\rho = J \rho K^\dagger - \tfrac12 \left\{K^\dagger J, \rho\right\},
```

with ``\mathcal{D}[J] = \mathcal{D}[J, J]`` for an independent decay channel. Consider a two-level atom in a driven cavity,

```math
H = \omega_c a^\dagger a + \omega_a \sigma^{ee} + g\left(a^\dagger \sigma^{ge} + \sigma^{eg} a\right) + \eta\left(e^{-i\omega_d t} a^\dagger + e^{i\omega_d t} a\right),
```

with cavity loss ``\kappa`` and spontaneous emission ``\gamma``.

```@contents
Pages = ["open_systems.md"]
Depth = 2:3
```

## Lindblad generators

[`hamiltonian_action`](@ref) builds the coherent part ``-i[H, \cdot]`` and [`dissipator`](@ref) a decay channel, with the rate as a scalar prefactor. The generator keeps its structure and prints as the Hamiltonian and its channels.

```@example open
using SecondQuantizedAlgebra
h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
a = Destroy(h, :a, 1)
σ(i, j) = Transition(h, :σ, i, j, 2)
@variables t::Real ωc::Real ωa::Real ωd::Real g::Real η::Real κ::Real γ::Real

H = ωc * a' * a + ωa * σ(:e, :e) + g * (a' * σ(:g, :e) + σ(:e, :g) * a) +
    η * (expim(-ωd * t) * a' + expim(ωd * t) * a)
L = hamiltonian_action(H) + κ * dissipator(a) + γ * dissipator(σ(:g, :e))
```

```@docs
Lindbladian
hamiltonian_action
dissipator
```

Different spellings of one generator compare equal. A numeric factor of a jump moves into the rate, and a constant in a jump is a coherent drive, ``\mathcal{D}[J + \alpha] = \mathcal{D}[J] - i\left[\tfrac{i}{2}(\bar\alpha J - \alpha J^\dagger), \cdot\right]``, so it moves into the Hamiltonian.

```@example open
dissipator(2a) == 4 * dissipator(a)
```

```@example open
dissipator(a + 3)
```

## Rotating frame

The drive oscillates at ``\omega_d``. [`transform`](@ref) moves the generator into the frame rotating with ``\omega_d (a^\dagger a + \sigma^{ee})``: the Hamiltonian is transformed and every jump is conjugated. The drive becomes static, the frequencies become the detunings ``\Delta_c = \omega_c - \omega_d`` and ``\Delta_a = \omega_a - \omega_d``, and the jumps are unchanged, since a loss channel only picks up a phase.

```@example open
Lrot = transform(L, RotatingFrame(ωd * (a' * a + σ(:e, :e)), t))
```

## Equations of motion

The transpose of a generator is its Heisenberg-picture form, ``\operatorname{Tr}(A\,\mathcal{L}\rho) = \operatorname{Tr}(\mathcal{L}^T\!A\,\rho)``, so ``\frac{d}{dt}\langle A \rangle = \langle \mathcal{L}^T A \rangle``. The field is damped at ``\kappa/2``, driven by ``\eta`` and coupled to the atomic coherence.

```@example open
average(transpose(Lrot)(a))
```

The coherence decays at ``\gamma/2`` and is driven by the field times the inversion.

```@example open
average(transpose(Lrot)(σ(:g, :e)))
```

For a weak drive the atom stays near its ground state, ``\langle a\sigma^{gg}\rangle \approx \langle a \rangle`` and ``\langle a\sigma^{ee}\rangle \approx 0``, and the two equations close. Their steady state is the vacuum-Rabi response of the cavity field,

```math
\langle a \rangle = \frac{-i\eta}{\kappa/2 + i\Delta_c + g^2/(\gamma/2 + i\Delta_a)}.
```

Trace preservation reads ``\mathcal{L}^T 1 = 0``.

```@example open
transpose(Lrot)(1)
```

```@docs
transpose(::SuperOperator)
```

## Numerics

[`to_numeric`](@ref) converts a generator to the backend's superoperator, with the same `parameter`, `time_parameter` and indexed `sites` arguments as for operators. The steady state of the rotating-frame Liouvillian at weak drive reproduces the vacuum-Rabi amplitude.

```@example open
import QuantumToolbox as QT
p = Dict(ωc => 1.0, ωa => 1.1, ωd => 0.95, g => 0.2, η => 0.002, κ => 0.1, γ => 0.05)
backend = QuantumToolboxBackend()
Lnum = to_numeric(Lrot, h, (10, 2); backend, parameter = p)
anum = to_numeric(a, h, (10, 2); backend)
Δc, Δa = 0.05, 0.15
QT.expect(anum, QT.steadystate(Lnum)), -0.002im / (0.05 + Δc * im + 0.2^2 / (0.025 + Δa * im))
```

Solvers that take a Hamiltonian and collapse operators get them from [`to_numeric_master`](@ref), here ``\sqrt\kappa\, a`` and ``\sqrt\gamma\, \sigma^{ge}``.

```@example open
using SecondQuantizedAlgebra: to_numeric_master
Hnum, c_ops = to_numeric_master(Lrot, h, (10, 2); backend, parameter = p)
QT.expect(anum, QT.steadystate(Hnum, c_ops))
```

```@docs
to_numeric_master
```

## Displacement frames

A damped linear mode is displaced to its steady state. [`DisplacementFrame`](@ref) takes the amplitude from the fixed point of the full drift, so loss and thermal heating both enter it, and the drive is gone in the displaced frame.

```@example displaced
using SecondQuantizedAlgebra
b = Destroy(FockSpace(:c), :b)
Lth = hamiltonian_action(2 * b' * b + 3 * (b + b')) + (1 // 2) * dissipator(b) + (1 // 10) * dissipator(b')
U = DisplacementFrame(b, Lth)
```

```@example displaced
transpose(transform(Lth, U))(b)
```

## Superoperators

A Lindbladian is one kind of linear map on operators. The general polynomial map ``\rho \mapsto \sum c\, A \rho B`` is a [`SuperOperator`](@ref), built from the sandwich [`sprepost`](@ref) and its one-sided forms [`spre`](@ref) and [`spost`](@ref). `SuperOperator(L)` lowers a Lindbladian to its sandwiches.

```@example open
using SecondQuantizedAlgebra: SuperOperator, spre, spost, sprepost
D = SuperOperator(dissipator(a))
```

```@example open
D == sprepost(a, a') - (1 // 2) * spre(a' * a) - (1 // 2) * spost(a' * a)
```

These names are public but not exported, because QuantumOpticsBase and QuantumToolbox export the same names for their numerical superoperators. With a backend loaded, its `spre` on a symbolic operator forwards to the symbolic one.

Maps compose with `∘` or `*`, and the right factor acts first. Composing, transposing or applying a Lindbladian gives a `SuperOperator`.

```@example open
(spre(a) ∘ spost(a'))(σ(:e, :e))
```

Every Lindbladian preserves the trace. The identity map does not, so it is not a Lindbladian: `L^0`, like a Lindbladian plus a scalar, is a general map. It prints as the bare placeholder of ``\rho \mapsto \rho``.

```@example open
L^0
```

`transpose` and `adjoint` are different maps. The adjoint belongs to the Hilbert–Schmidt product ``\operatorname{Tr}(X^\dagger S\rho)`` and takes the adjoint of both sides, so the two agree only for maps that preserve Hermiticity.

```@example open
transpose(spre(a))(σ(:e, :e)), adjoint(spre(a))(σ(:e, :e))
```

```@docs
SuperOperator
sprepost
spre
spost
compose
adjoint(::SuperOperator)
⊗(::SuperOperator, ::SuperOperator)
```

## GKSL form of a map

A map written in sandwiches, or one that comes out of a composition or a frame change, is read back as a `Lindbladian`. The extraction is exact in the operator algebra and needs no Hilbert-space cutoff.

```@example open
Lindbladian(κ * (sprepost(a, a') - (1 // 2) * spre(a' * a) - (1 // 2) * spost(a' * a)))
```

```@example open
Lindbladian(SuperOperator(Lrot)) == Lrot
```

The Kossakowski matrix depends on the basis of jumps. In the basis ``(a, a + a^\dagger a)``, independent loss and dephasing become a correlated matrix.

```@example open
kossakowski(κ * dissipator(a) + γ * dissipator(a' * a), JumpBasis(a, a + a' * a))
```

Basis operators are taken modulo the identity and with the completeness relation expanded.

```@example open
JumpBasis(σ(:g, :g), a + 1)
```

```@docs
Lindbladian(::SuperOperator, ::JumpBasis)
JumpBasis
kossakowski
support_basis
independent_basis
GKSLFormError
```

## Correlated and collective decay

A cross channel ``\mathcal{D}[J, K]`` with ``J \neq K`` describes two decay paths into a shared bath, such as the cavity and the atom radiating into the same outgoing mode. [`liouvillian`](@ref) builds a generator from a list of jumps and a rate matrix in the QuantumOptics convention, where ``\Gamma_{ij}`` multiplies ``J_i \rho J_j^\dagger``. The off-diagonal entries become cross channels.

```@example open
using SecondQuantizedAlgebra: liouvillian
@variables Γc
liouvillian(0, [a, σ(:g, :e)]; rates = [κ Γc; Γc γ])
```

```@docs
liouvillian
```

For an ensemble of emitters, a generator is summed over sites with [`Σ`](@ref). Independent emission at rate ``\gamma`` on every atom drains the total excitation at rate ``\gamma``.

```@example many
using SecondQuantizedAlgebra
ha = NLevelSpace(:atom, (:g, :e))
@variables N γ
i = Index(ha, :i, N, ha)
j = Index(ha, :j, N, ha)
k = Index(ha, :k, N, ha)
s(x, y, n) = IndexedOperator(Transition(ha, :s, x, y), n)

Lind = Σ(γ * dissipator(s(:g, :e, i)), i)
transpose(Lind)(Σ(s(:e, :e, j), j))
```

Emitters closer than a wavelength radiate collectively, with a rate matrix ``\Gamma_{ij}`` over pairs of sites. `dissipator(J, Γ, i => j)` builds ``\sum_{i,j} \Gamma_{ij} \mathcal{D}[J_i, J_j]`` and splits off the diagonal ``i = j``, so the generator shows its local and correlated parts.

```@example many
Lcol = dissipator(s(:g, :e, i), DoubleIndexedVariable(:Γ, i, j), i => j)
```

The population of one emitter ``k`` then couples to the coherences between ``k`` and every other emitter.

```@example many
average(transpose(Lcol)(s(:e, :e, k)))
```

A map shares one summation scope between its two sides, as in ``\sum_i J_i \rho J_i^\dagger``. Constructions that would make the scope ambiguous throw instead of guessing, and the algebra never invents an index to resolve them:

- composing two maps that bind the same index; rename one with [`change_index`](@ref);
- applying a map to an operator in which its bound index is free, since the sum would capture it, which is why the drift above is taken on `s(:e, :e, k)`;
- a jump that carries its own `Σ`, since ``J \rho J^\dagger`` with ``J = \sum_i J_i`` needs two independent copies of the index;
- powers of a map with bound indices, for the same reason.

```@repl many
Lind ∘ Lind
dissipator(Σ(s(:g, :e, i), i))
```

## Formal operator functions

A Hamiltonian or a jump with a formal function such as `cos`, `sin` or `expim` gives a `Lindbladian{QExpr}`, as for a transmon with its cosine potential and dephasing of the phase. [`taylor`](@ref) lowers it to a polynomial generator.

```@example transmon
using SecondQuantizedAlgebra
φ = Position(PhaseSpace(:transmon), :φ)
@variables γφ
Lt = hamiltonian_action(-cos(φ)) + γφ * dissipator(φ)
```

```@example transmon
taylor(Lt, 0:2)
```

A general map with a formal side is a [`SuperExpr`](@ref), and `taylor` lowers it in the same way.

```@example transmon
using SecondQuantizedAlgebra: sprepost
taylor(sprepost(cos(φ), φ), 0:2)
```

```@docs
SuperExpr
```

## Inspecting and building maps

The parts of a generator are read back with accessors. [`hamiltonian`](@ref) and `channels` split a Lindbladian into its Hamiltonian and its `Channel => rate` pairs, and [`jump`](@ref), `partner` and `scope` read a channel.

```@example open
using SecondQuantizedAlgebra: hamiltonian, channels, jump, partner, scope
ch, rate = first(channels(Lrot))
jump(ch), partner(ch), rate
```

Iterating a `SuperOperator` gives its `SuperTerm => coefficient` pairs, whose sides are read with [`left_operator`](@ref) and `right_operator`. A `SuperExpr` splits into its polynomial part and its formal sandwiches with [`poly`](@ref) and `formal_terms`.

```@example open
using SecondQuantizedAlgebra: left_operator, right_operator
term, c = first(D)
left_operator(term), right_operator(term), c
```

`sum` over an array of Lindbladians or `SuperOperator`s accumulates in place, like `sum` over operators, so a generator with many channels costs linear time. A loop uses a builder directly.

```@example open
using SecondQuantizedAlgebra: LindbladianBuilder, QAdd, accumulate!, build
hm = reduce(⊗, [FockSpace(Symbol(:mode, n)) for n in 1:3])
modes = [Destroy(hm, Symbol(:a, n), n) for n in 1:3]
builder = LindbladianBuilder{QAdd}()
for m in modes
    accumulate!(builder, κ * dissipator(m))
end
build(builder) == sum([κ * dissipator(m) for m in modes])
```

```@docs
hamiltonian
Channel
jump
SuperTerm
superterm
left_operator
poly
SuperBuilder
LindbladianBuilder
accumulate!
build
```
