# # Dicke Superradiance with Collective Transitions
#
# `N` identical two-level atoms coupled to a single cavity mode form the
# **Tavis-Cummings / Dicke** system.  When all atoms see the same field, the
# dynamics preserves permutation symmetry — states stay in the *symmetric
# (bosonic) subspace*, whose dimension grows polynomially in `N` rather than
# as `2^N`.  This is exactly the regime where
# [`CollectiveTransition`](@ref SecondQuantizedAlgebra.CollectiveTransition)
# operators ``S^{ij} = \sum_k \sigma_k^{ij}`` are the natural variables.
#
# In this example we:
# 1. Build the Tavis-Cummings Hamiltonian with collective atom operators,
# 2. Verify the ``\mathfrak{su}(N)`` Lie-algebra commutators that drive the
#    dynamics,
# 3. Propagate the fully-inverted state ``|\text{vac}\rangle \otimes |e\rangle^{\otimes N}``
#    and observe the **superradiant burst** — the cavity occupation builds up
#    faster than for `N` independent atoms thanks to the ``\sqrt{N}`` enhanced
#    coupling.
#
# See [Collective Atoms — Two Regimes](@ref) for when to use
# `CollectiveTransition` versus indexed [`Σ`](@ref SecondQuantizedAlgebra.Σ)
# on a `ProductSpace`.

# ## Setup
#
# A single cavity mode plus an `NLevelSpace` for the (collective) two-level atoms.
# The atomic Hilbert space is a *single* `NLevelSpace` — the collective
# operator already represents the sum over `N` atoms; we don't build a product
# of `N` atom spaces.

using SecondQuantizedAlgebra

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g, :e))
h = hc ⊗ ha

@qnumbers a::Destroy(h, 1)

S(α, β) = CollectiveTransition(h, :S, α, β, 2)

@variables ωc ωa g

# ## Tavis-Cummings Hamiltonian
#
# In the rotating-wave approximation,
#
# ```math
# H = \omega_c\, a^\dagger a
#   + \omega_a\, S^{ee}
#   + g\bigl(a^\dagger S^{ge} + a\, S^{eg}\bigr).
# ```

H = ωc * a' * a + ωa * S(:e, :e) + g * (a' * S(:g, :e) + a * S(:e, :g))

# ## Lie-algebra structure
#
# The collective operators satisfy
# ``[S^{ij}, S^{kl}] = \delta_{jk}\, S^{il} - \delta_{li}\, S^{kj}``.  Under the
# default [`NormalOrder`](@ref SecondQuantizedAlgebra.NormalOrder) convention
# this commutator fires eagerly during multiplication, so derivations stay in
# closed collective form rather than expanding into explicit per-atom sums.
# A quick sanity check:

commutator(S(:e, :g), S(:g, :e))   # = S^{ee} − S^{gg}

# Equations of motion follow directly.  Cavity field:

i_dadt = 1im * commutator(H, a)

# Atomic coherence ``S^{ge}``:

i_dSgedt = 1im * commutator(H, S(:g, :e))

# Notice that the right-hand sides stay compact — every term is either a
# `CollectiveTransition` or a Fock operator times one.  The same derivation in
# the indexed representation (a single `IndexedOperator` per atom on a
# `ProductSpace`) would carry around `Σ` symbols with explicit non-equality
# constraints.

# ## Hilbert-space scaling
#
# The symmetric subspace dimension is ``\binom{N + n - 1}{n - 1}``, polynomial
# in `N`, whereas the full product space scales as ``n^N``.  For two-level
# atoms (`n = 2`) the contrast is striking: `N = 100` atoms is `2^100 ≈ 10^30`
# in the full product space — entirely intractable — but only `101` in the
# symmetric subspace.

# ## Numerical simulation: superradiant burst
#
# We propagate the fully inverted state under the resonant Hamiltonian
# (``\omega_c = \omega_a``) and watch the cavity occupation build up.
# The atomic basis is the symmetric subspace built from `bosonstates`.

using QuantumOpticsBase
using OrdinaryDiffEqTsit5
using CairoMakie

function build_problem(N; n_cav_max = N + 4, ω = 1.0, g_val = 0.05)
    b_cav = FockBasis(n_cav_max)
    b_atom = NLevelBasis(2)
    b_collective = ManyBodyBasis(b_atom, bosonstates(b_atom, N))
    b = b_cav ⊗ b_collective

    subs = Dict(ωc => ω, ωa => ω, g => g_val)
    H_num = to_numeric(substitute(H, subs), b)

    # Initial state: cavity vacuum, all atoms excited.
    # In `bosonstates` ordering, [N - k, k] has k atoms in level 2 (|e⟩).
    ψ0 = fockstate(b_cav, 0) ⊗ basisstate(b_collective, [0, N])

    n_cav_op = to_numeric(a' * a, b)

    return H_num, ψ0, n_cav_op
end

# We sweep three system sizes.  Time is given in units of the inverse coupling
# ``g^{-1}`` so that the single-atom Rabi period is ``2\pi``.

times = collect(range(0, 80.0; length = 401))
g_val = 0.05

results = Dict{Int, Vector{Float64}}()

for N in (1, 5, 20)
    H_num, ψ0, n_cav_op = build_problem(N; g_val = g_val)
    rhs(dψ, ψ, _, _) = (dψ .= -im .* (H_num * Ket(ψ0.basis, ψ)).data; nothing)
    prob = ODEProblem(rhs, ψ0.data, (times[1], times[end]))
    sol = solve(prob, Tsit5(); saveat = times, abstol = 1.0e-10, reltol = 1.0e-10)
    n_cav_t = [real(QuantumOpticsBase.expect(n_cav_op, Ket(ψ0.basis, ψ))) for ψ in sol.u]
    results[N] = n_cav_t
end

# Plot the cavity occupation:

fig = Figure(; size = (640, 360))
ax = Axis(
    fig[1, 1];
    xlabel = "time  (g·t)",
    ylabel = "⟨a†a⟩(t)",
    title = "Tavis-Cummings build-up: collective coupling enhances Rabi rate",
)
for N in (1, 5, 20)
    n_cav_t = results[N]
    lines!(ax, g_val .* times, n_cav_t; label = "N = $N", linewidth = 2)
end
axislegend(ax; position = :rb)
fig

# The first peak in the cavity occupation arrives faster as `N` grows, with the
# Rabi frequency scaling as ``g\sqrt{N}`` — the hallmark of collective
# light-matter coupling.  In the `CollectiveTransition` representation this is
# direct: one `NLevelSpace`, one cavity, polynomial-size symmetric subspace,
# regardless of `N`.
