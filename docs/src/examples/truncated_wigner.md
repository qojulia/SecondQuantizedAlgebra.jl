```@meta
EditURL = "../../../examples/truncated_wigner.jl"
```

# Truncated Wigner for a Driven-Dissipative Kerr Resonator

The **truncated Wigner approximation** (TWA) maps a driven-dissipative
quantum system onto an ensemble of stochastic *classical* trajectories in
phase space. It keeps the leading quantum fluctuations that mean-field theory
throws away, at the cost of a classical Langevin equation rather than a full
master equation. TWA earns its keep when the drive populates the modes with
many photons, so that a master equation needs a large Fock cutoff, and it
scales to lattices where the joint Fock space grows as the *product* of the
per-mode cutoffs. Here we take the minimal such system, a single mode, so the
recipe stays readable; nothing below is specific to one mode.

The physical system is a **driven Kerr resonator**: one lossy Kerr mode ``a``
with a coherent drive ``F``. In the frame rotating at the drive frequency
(detuning ``\Delta``),

```math
H = -\Delta\, a^\dagger a
  + \frac{U}{2}\, a^\dagger a^\dagger a a
  + F\,(a^\dagger + a),
```

with single-photon loss ``\kappa`` (Lindblad jump operator
``L = \sqrt{\kappa}\,a``). This model is **bistable**: over a range of drives
the mean-field equation has two coexisting stable amplitudes. We will see that
mean-field, launched from vacuum, gets trapped on the low-amplitude branch,
while TWA reveals the quantum fluctuations that drive switching to the
high-amplitude branch. That switching is invisible to any single classical
trajectory.

The star of this example is [`normal_to_symmetric`](@ref). Wigner-function
moments are exactly the **symmetrically-ordered** operator averages, so the
TWA c-number equations follow from one uniform recipe:

> rewrite the Heisenberg drift in symmetric order, take the [`average`](@ref),
> factorize into a product of first moments.

## Setup

````@example truncated_wigner
using SecondQuantizedAlgebra
using Symbolics, SymbolicUtils

h = FockSpace(:cavity)

@qnumbers a::Destroy(h)
@variables Δ U F κ

H = -Δ * (a' * a) + (U / 2) * (a' * a' * a * a) + F * (a' + a)
````

## The Wigner shift comes out of `normal_to_symmetric`

The whole approximation hinges on one algebraic fact. The normal-ordered Kerr
product, rewritten in symmetric (Weyl) order, picks up a shift:

````@example truncated_wigner
normal_to_symmetric(a' * a * a)
````

The trailing ``-a`` is the Wigner correction. Under TWA it becomes the famous
``(|\alpha|^2 - 1)\,\alpha`` term. The same mechanism corrects observables:
the occupation operator ``a^\dagger a`` maps to ``|\alpha|^2 - \tfrac12``,

````@example truncated_wigner
normal_to_symmetric(a' * a)
````

so a physical photon number is recovered from the phase-space cloud as
``\langle n\rangle = \langle |\alpha|^2\rangle - \tfrac12``.

## From operators to a c-number drift

The TWA drift is the Heisenberg-Langevin equation, symmetrized, averaged, then
factorized to first order (the mean-field truncation). The factorization
``\langle a^\dagger a a\rangle \to \langle a^\dagger\rangle
\langle a\rangle \langle a\rangle`` is the one step SecondQuantizedAlgebra
does not ship as a built-in, so we spell it out: walk the averaged expression
and split every multi-operator moment into a product of single-operator
averages.

````@example truncated_wigner
function factorize_meanfield(expr)
    rule = SymbolicUtils.PassThrough(
        function (x)
            is_average(x) || return nothing
            q = undo_average(x)
            terms = collect(q.arguments)
            length(terms) == 1 || return nothing
            ops = first(terms)[1].ops
            length(ops) <= 1 && return nothing
            return SymbolicUtils.unwrap(prod(average(op) for op in ops))
        end,
    )
    return SymbolicUtils.Postwalk(rule)(SymbolicUtils.unwrap(expr))
end
````

The loss enters as the standard drift ``-\tfrac{\kappa}{2} a`` of a lowering
operator. The mean-field and TWA drifts differ *only* by the
`normal_to_symmetric` call:

````@example truncated_wigner
heisenberg(o) = -1im * commutator(o, H) - (κ / 2) * o

meanfield_drift(o) = factorize_meanfield(average(heisenberg(o)))
twa_drift(o) = factorize_meanfield(average(normal_to_symmetric(heisenberg(o))))

twa_drift(a)
````

Subtracting the two confirms the surplus is exactly the Wigner shift
``+iU\langle a\rangle``:

````@example truncated_wigner
Symbolics.expand(twa_drift(a) - meanfield_drift(a))
````

## Porting to a ModelingToolkit system

The stochastic variable is the average itself. [`make_time_dependent`](@ref)
lifts ``\langle a\rangle`` into a ModelingToolkit unknown
``\langle a\rangle(t)``, and the conjugate average ``\langle a^\dagger\rangle``
becomes ``\mathrm{conj}(\langle a\rangle(t))``. We keep the unknown
**complex** and let StochasticDiffEq carry complex state, so no real/imaginary
split is needed.

````@example truncated_wigner
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

α = make_time_dependent(average(a), t)
````

SecondQuantizedAlgebra represents the imaginary unit inside an `average` as a
symbolic constant; fold it to a numeric `im`, and map the averages to the
time-dependent unknown, so ModelingToolkit sees a plain complex ODE.

````@example truncated_wigner
function to_rhs(drift)
    numeric = Symbolics.substitute(drift, Dict(SymbolicUtils.unwrap(Symbolics.IM) => 1im))
    return Symbolics.substitute(numeric, Dict(average(a) => α, average(a') => conj(α)))
end

eqs = [D(α) ~ to_rhs(twa_drift(a))]
````

Assemble the (deterministic) drift system. Its generated right-hand side is
reused below as the drift of the stochastic problem.

````@example truncated_wigner
@named twa = System(eqs, t, [α], [Δ, U, F, κ])
twa = mtkcompile(twa)
````

## The noise

Vacuum single-photon loss contributes a diagonal Wigner diffusion. In the
complex representation the mode gets a complex noise of amplitude
``\sqrt{\kappa/2}``. This is the one number fixed by convention rather than by
the algebra: it is set so that a bare lossy mode relaxes to the vacuum width
``\langle|\alpha|^2\rangle = \tfrac12``, matching the ``-\tfrac12`` in the
occupation map above.

````@example truncated_wigner
using OrdinaryDiffEq, StochasticDiffEq
using SymbolicIndexingInterface: getp

getκ = getp(twa, κ)

function noise!(du, u, p, t)
    du[1] = sqrt(getκ(p) / 2)
    return nothing
end
````

## Mean-field: two coexisting branches

The mean-field limit drops both quantum corrections (the Wigner shift and the
noise): it is the deterministic drift built from `meanfield_drift`. We
integrate it to steady state from two different starts, from vacuum and from a
large-amplitude seed.

````@example truncated_wigner
using Statistics

params = Dict(Δ => 2.0, U => 0.1, F => 3.0, κ => 1.0)
tspan = (0.0, 400.0)

mf_eqs = [D(α) ~ to_rhs(meanfield_drift(a))]
@named mf = System(mf_eqs, t, [α], [Δ, U, F, κ])
mf = mtkcompile(mf)

mfsolve(α0) = solve(ODEProblem(mf, merge(Dict(α => α0), params), tspan), Tsit5())
αlo = mfsolve(0.0im).u[end][1]        # from vacuum
αhi = mfsolve(6.0 + 0im).u[end][1]    # from a bright seed
(abs2(αlo), abs2(αhi))                # occupation on each branch
````

Two stable amplitudes coexist for the same parameters. A single classical
trajectory, released from vacuum, can only ever find the lower one.

## Solving the TWA ensemble

TWA replaces that one trajectory by an ensemble. Each realization starts from
a **vacuum-sampled** initial condition (the Wigner distribution of the vacuum
is a complex Gaussian with ``\langle|\alpha|^2\rangle = \tfrac12``) and is
driven by the loss noise. The drift is reused verbatim from the deterministic
problem via `SDEFunction`.

````@example truncated_wigner
prob = ODEProblem(twa, merge(Dict(α => 0.0im), params), tspan)
sdeprob = SDEProblem(SDEFunction(prob.f.f, noise!), prob.u0, tspan, prob.p)

vacuum_start(prob, ctx) = remake(prob; u0 = randn(ComplexF64, 1) ./ sqrt(2))

ntraj = 5_000
saveat = 0.0:2.0:400.0
ensemble = EnsembleProblem(sdeprob; prob_func = vacuum_start)
sim = solve(ensemble, LambaEulerHeun(), EnsembleThreads(); trajectories = ntraj, saveat)
nothing #hide
````

Each trajectory is a phase-space sample ``\alpha(t)``. The physical occupation
is the ensemble mean, minus the vacuum ``\tfrac12``:

````@example truncated_wigner
traj = [reduce(vcat, s.u) for s in sim.u]
n_twa = [mean(abs2(traj[k][i]) for k in 1:ntraj) - 0.5 for i in eachindex(saveat)]
nothing #hide
````

## Occupation: TWA versus the mean-field branches

````@example truncated_wigner
using CairoMakie

fig = Figure(; size = (470, 330))
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "⟨n⟩")
hlines!(ax, [abs2(αlo)]; color = :gray, linestyle = :dash, label = "mean-field (lower)")
hlines!(ax, [abs2(αhi)]; color = :gray, linestyle = :dot, label = "mean-field (upper)")
lines!(ax, collect(saveat), n_twa; color = :firebrick, label = "TWA ⟨n⟩")
axislegend(ax; position = :lt)
fig
````

Mean-field launched from vacuum sits forever on the lower branch (dashed). The
TWA ensemble mean (red) climbs steadily above it: quantum fluctuations kick a
growing fraction of trajectories over the barrier onto the upper branch
(dotted), so the average drifts up between the two, and it is still climbing at
``t = 400``. Nothing in the deterministic dynamics from vacuum could produce
this.

## The phase-space (Wigner) distribution

Because every trajectory is a phase-space point, the ensemble at the final
time is a sampled Wigner quasi-distribution. We bin it with
[FHist](https://github.com/Moelf/FHist.jl) and overlay both mean-field fixed
points.

````@example truncated_wigner
using FHist

reα = [real(traj[k][end]) for k in 1:ntraj]
imα = [imag(traj[k][end]) for k in 1:ntraj]
hist = Hist2D((reα, imα); binedges = (-7:0.2:5, -7:0.2:5))

fig2 = Figure(; size = (470, 410))
ax = Axis(fig2[1, 1]; xlabel = "Re α", ylabel = "Im α", aspect = 1, title = "Wigner samples (t = 400)")
heatmap!(ax, bincenters(hist)..., bincounts(hist); colormap = :magma)
scatter!(ax, [real(αlo)], [imag(αlo)]; color = :cyan, marker = :xcross, markersize = 16, label = "mean-field (lower)")
scatter!(ax, [real(αhi)], [imag(αhi)]; color = :lime, marker = :cross, markersize = 16, label = "mean-field (upper)")
axislegend(ax; position = :rb)
fig2
````

The distribution is unmistakably **bimodal**: one lobe sits on each mean-field
fixed point, the upper one smeared into an arc by phase diffusion, with a faint
bridge of trajectories caught mid-switch between them. This is the
picture mean-field cannot draw. A direct master-equation solve would need a
Fock cutoff comfortably above the upper-branch occupation, and the same recipe
on a lattice of such modes would need that cutoff on *every* site, whereas TWA
needed only classical trajectories.

Every ingredient downstream of the Hamiltonian, the drift, the observable map
``\langle n\rangle = \langle|\alpha|^2\rangle - \tfrac12``, and even the sign
of the Kerr shift, came out of symmetric ordering plus `average`, with no hand
conversion from operators to c-numbers.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

