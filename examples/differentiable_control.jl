# # Differentiable Quantum Optimal Control (QuantumToolbox + Mooncake)
#
# QuantumToolbox.jl can differentiate the Schrödinger and master-equation solvers
# with respect to their parameters, using the SciMLSensitivity adjoint machinery
# together with a reverse-mode backend such as Mooncake or Enzyme. That makes
# gradient-based **quantum optimal control** possible: shape a control pulse, run the
# dynamics, measure how far the final state is from a target, and descend the gradient
# of that error with respect to the pulse parameters. This capability is specific to the
# QuantumToolbox backend; the QuantumOpticsBase backend has no solver-parameter vector to
# differentiate against.
#
# `SecondQuantizedAlgebra` contributes the piece that is otherwise missing: it lets you
# write the entire time-dependent control Hamiltonian **symbolically** and hand back a
# differentiable numeric operator. A `time_parameter` value function may be written as
# `(p, t) -> value`, a function of the solver parameter vector `p` and time `t`. When any
# value function reads `p`, the resulting `QobjEvo` threads `p` into its coefficients and
# stays differentiable. No manual operator assembly is required.
#
# The task below is single-photon state preparation in a Kerr resonator. The Kerr
# nonlinearity shifts the ``1 \to 2`` transition out of resonance (photon blockade), so a
# resonant drive can pump ``0 \to 1`` but not beyond. We optimize the drive pulse to reach
# the Fock state ``|1\rangle`` from vacuum with high fidelity.

# ## Setup
#
# The Hamiltonian in the frame rotating at the drive frequency is
#
# ```math
# H(t) = \Delta\, a^\dagger a + \frac{K}{2}\, a^\dagger a^\dagger a a
#      + \varepsilon(t)\,(a^\dagger + a),
# ```
#
# with detuning ``\Delta``, Kerr strength ``K``, and a real control ``\varepsilon(t)``.

using SecondQuantizedAlgebra
using QuantumToolbox
using Symbolics

using SciMLSensitivity
import Mooncake                       # loads the Mooncake reverse-mode backend
import DifferentiationInterface as DI
using DifferentiationInterface: AutoMooncake
using Optimization, OptimizationOptimJL
using Optimization, OptimizationOptimisers

h = FockSpace(:c)
@qnumbers a::Destroy(h)
@variables Δ K ε

# Physical constants and control discretization. The pulse is a Fourier sine series, which
# is smooth and vanishes at both ends of the interval, and its amplitudes are the parameters
# we optimize.

N = 10                       # Fock cutoff (Hilbert dimension N + 1)
Δ0, K0 = 0.0, -2.0           # on resonance for 0 -> 1; Kerr blockade detunes 1 -> 2
T = 6.0                      # control duration
tlist = collect(range(0, T, 61))
M = 6                        # number of control parameters

pulse(p, t) = sum(p[k] * sinpi(k * t / T) for k in eachindex(p))

# ## Building the differentiable operator
#
# The static part of the Hamiltonian is fixed with `parameter` (substituted eagerly). The
# drive is registered through `time_parameter` as a two-argument function of `(p, t)`. That
# `(p, t)` arity is what marks the coefficient as solver-parameter aware, so the returned
# `QobjEvo` is differentiable with respect to `p`. The operator is built **once**; each
# `sesolve` call supplies the current `p`.

Hc = to_numeric(
    Δ * a' * a + (K / 2) * a' * a' * a * a + ε * (a' + a), h, N;
    backend = QuantumToolboxBackend(),
    parameter = Dict(Δ => Δ0, K => K0),
    time_parameter = Dict(ε => (p, t) -> pulse(p, t)),
)

# ## Cost function
#
# The cost is the infidelity ``1 - |\langle 1 | \psi(T) \rangle|^2``, i.e. one minus the
# population in ``|1\rangle`` at the final time. `sesolve` receives `params = p` and the
# adjoint choice verified to give correct gradients here, `BacksolveAdjoint` with the
# Mooncake vector-Jacobian product.

ψ0 = fock(N + 1, 0)
proj1 = fock(N + 1, 1) * fock(N + 1, 1)'
sensealg = BacksolveAdjoint(autojacvec = SciMLSensitivity.MooncakeVJP())

function infidelity(p)
    sol = sesolve(Hc, ψ0, tlist; params = p, sensealg = sensealg, progress_bar = Val(false))
    return 1 - real(expect(proj1, sol.states[end]))
end

# ## Gradient check
#
# Before optimizing, confirm the Mooncake gradient agrees with a finite-difference estimate.

p0 = fill(0.3, M)
g_ad = DI.gradient(infidelity, AutoMooncake(), p0)

δ = 1.0e-5
g_fd = [(infidelity(p0 .+ δ .* (1:M .== k)) - infidelity(p0 .- δ .* (1:M .== k))) / (2δ) for k in 1:M]

maximum(abs, g_ad .- g_fd) < 1.0e-6

# ## Optimization
#
# `Optimization.jl` is the idiomatic SciML interface: wrap the cost in an `OptimizationFunction`
# that carries the AD backend, build an `OptimizationProblem`, and `solve`. For a smooth,
# deterministic, low-dimensional control problem with exact gradients, the quasi-Newton `LBFGS`
# converges in far fewer solver evaluations than a first-order method such as Adam. A callback
# records the loss for the convergence plot.

losses = Float64[]
record = (state, l) -> (push!(losses, l); false)

optf = OptimizationFunction((u, _p) -> infidelity(u), AutoMooncake())
prob = OptimizationProblem(optf, fill(0.3, M))
result = solve(prob, OptimizationOptimJL.LBFGS(); maxiters = 100, callback = record)
p_opt = result.u

1 - result.objective   # final fidelity

# ## Results
#
# The optimized pulse drives the mean photon number from 0 to 1, and the final population
# lands almost entirely in ``|1\rangle``; the Kerr blockade suppresses ``|2\rangle``.

using CairoMakie

sol = sesolve(Hc, ψ0, tlist; params = p_opt, e_ops = [create(N + 1) * destroy(N + 1)], progress_bar = Val(false))
nt = real.(sol.expect[1, :])
εt = pulse.(Ref(p_opt), tlist)
finalpops = [real(expect(fock(N + 1, n) * fock(N + 1, n)', sol.states[end])) for n in 0:4]

fig = Figure(size = (1050, 320))

ax1 = Axis(fig[1, 1]; xlabel = "iteration", ylabel = "infidelity 1 − F", yscale = log10, title = "LBFGS convergence")
lines!(ax1, 1:length(losses), losses; color = :crimson)

ax2 = Axis(fig[1, 2]; xlabel = "t", ylabel = "⟨n⟩ / ε(t)", title = "optimized pulse & dynamics")
lines!(ax2, tlist, nt; label = "⟨n⟩(t)", color = :navy)
lines!(ax2, tlist, εt; label = "ε(t)", color = :darkorange, linestyle = :dash)
axislegend(ax2; position = :lt)

ax3 = Axis(fig[1, 3]; xlabel = "Fock state n", ylabel = "population", title = "final state", xticks = 0:4)
barplot!(ax3, 0:4, finalpops; color = :seagreen)

fig

# ## Notes
#
# * The same mechanism extends to open systems. A collapse operator built with a p-aware
#   rate, `to_numeric(sqrt(κ) * a, h, N; time_parameter = Dict(κ => (p, t) -> p[1]))`, feeds
#   into `liouvillian` and `mesolve` and stays differentiable, so dissipation rates can be
#   optimized alongside the drive.
# * `BacksolveAdjoint` with the Mooncake vector-Jacobian product is used here because it was
#   verified to match finite differences for both `sesolve` and `mesolve`. Enzyme is also
#   supported. Zygote is not: it returns silently wrong gradients for time evolution.
# * The `(p, t)` control form is specific to the QuantumToolbox backend. On the
#   QuantumOpticsBase backend a two-argument `time_parameter` value raises an informative
#   error, since that backend has no solver parameter vector to differentiate against.
