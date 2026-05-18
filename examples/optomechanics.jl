# # Cavity Optomechanics and the Polaron Transformation
#
# A single optical cavity coupled to a mechanical oscillator by radiation
# pressure has the Hamiltonian
#
# ```math
# H = \omega_c\, a^\dagger a + \omega_m\, b^\dagger b
#   + g\, a^\dagger a\,(b + b^\dagger),
# ```
#
# where ``a, a^\dagger`` are cavity (optical) ladder operators, ``b, b^\dagger``
# the mechanical mode, and ``g`` the single-photon optomechanical coupling.
# Although ``H`` is non-linear, it has a hidden structure: the cavity photon
# number ``a^\dagger a`` commutes with everything in sight, so each
# photon-number sector sees the mechanical mode as a *displaced* harmonic
# oscillator.  The polaron (Lang-Firsov) transformation makes this exact and
# yields an effective **Kerr** nonlinearity for the cavity.
#
# The unitary that does the job is
#
# ```math
# U = \exp\!\left[\tfrac{g}{\omega_m}\,a^\dagger a\,(b - b^\dagger)\right].
# ```
#
# Because ``a^\dagger a`` commutes with the displacement generator,
# ``U^\dagger b U = b - (g/\omega_m)\,a^\dagger a`` is exact — the Hadamard
# series truncates at first order.  We exploit that by simply *substituting*
# the displaced operators into ``H`` and letting the package's eager
# normal-ordering do the bookkeeping.

# ## Setup

using SecondQuantizedAlgebra

hc = FockSpace(:cavity)
hm = FockSpace(:mech)
h = hc ⊗ hm

@qnumbers a::Destroy(h, 1)
@qnumbers b::Destroy(h, 2)

@variables ωc ωm g

H = ωc * a' * a + ωm * b' * b + g * a' * a * (b + b')

# ## Heisenberg dynamics: the radiation-pressure force
#
# The mechanical mode feels a force ``-g\,a^\dagger a`` from the optical
# field, while the optical mode is purely phase-modulated by mechanical
# position:

-1im * commutator(b, H)

#-

-1im * commutator(a, H)

# Reading off,
#
# ```math
# \dot b = -i\,\omega_m\, b - i\, g\, a^\dagger a, \qquad
# \dot a = -i\,\omega_c\, a - i\, g\, a\,(b + b^\dagger).
# ```
#
# The first equation makes the **displacement** structure explicit: at fixed
# ``a^\dagger a = n``, the mechanical mode oscillates around a shifted origin
# ``\langle b\rangle_\text{eq} = -(g/\omega_m)\,n``.

# ## Polaron transformation by operator substitution
#
# A unitary transformation acts on operators by conjugation: ``\tilde O =
# U^\dagger O U``.  Here that means
#
# ```math
# \tilde a = a, \qquad
# \tilde b = b - \tfrac{g}{\omega_m}\,a^\dagger a,
# \qquad
# \tilde b^\dagger = b^\dagger - \tfrac{g}{\omega_m}\,a^\dagger a,
# ```
#
# where the first identity holds because ``a^\dagger a`` is the generator
# of the transformation and commutes with ``a`` up to scalars that cancel
# (cf. the Hadamard series).  Substitute these into ``H`` and the package
# canonicalises the result:

α = g / ωm
b_tilde = b - α * a' * a
bd_tilde = b' - α * a' * a

H_pol = ωc * a' * a + ωm * bd_tilde * b_tilde + g * a' * a * (b_tilde + bd_tilde)

# ## Effective Kerr Hamiltonian
#
# Combining the canonical form of ``H_\mathrm{pol}`` with the identity
# ``(a^\dagger a)^2 = a^\dagger{}^2 a^2 + a^\dagger a`` gives
#
# ```math
# \boxed{\;
#   H_\mathrm{pol}
#   = \omega_c\, a^\dagger a + \omega_m\, b^\dagger b
#   - \frac{g^2}{\omega_m}\,(a^\dagger a)^2.
# \;}
# ```
#
# The optomechanical coupling has been **eliminated entirely** in favour of
# a Kerr nonlinearity of strength ``K = g^2/\omega_m`` on the cavity mode.
# Photons attract each other with energy ``-K\,n(n-1)`` — a self-Kerr
# anharmonicity that is now routinely measured in superconducting and
# membrane-in-the-middle optomechanical experiments.  The mechanical mode,
# in the polaron frame, has decoupled completely from the optics.

# ## Numerical verification
#
# We diagonalise both the original ``H`` and the Kerr-only effective
# Hamiltonian
# ``H_\mathrm{Kerr} = \omega_c\,a^\dagger a + \omega_m\,b^\dagger b - (g^2/\omega_m)\,(a^\dagger a)^2``
# and check that they share the same spectrum.

using QuantumOpticsBase, LinearAlgebra, CairoMakie

ωc_val, ωm_val = 1.0, 0.4
n_max_c, n_max_m = 4, 60
b_cav = FockBasis(n_max_c)
b_mech = FockBasis(n_max_m)
b_total = b_cav ⊗ b_mech

H_Kerr = ωc * a' * a + ωm * b' * b - (g^2 / ωm) * a' * a * a' * a

gs = range(0.0, 0.6 * ωm_val, length = 16)
E_full = Vector{Float64}[]
E_eff = Vector{Float64}[]
for g_val in gs
    subs = Dict(ωc => ωc_val, ωm => ωm_val, g => g_val)
    Hf = dense(to_numeric(substitute(H, subs), b_total))
    He = dense(to_numeric(substitute(H_Kerr, subs), b_total))
    push!(E_full, sort(real.(eigvals(Hermitian(Hf.data))))[1:8])
    push!(E_eff, sort(real.(eigvals(Hermitian(He.data))))[1:8])
end

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel = L"g / \omega_m",
    ylabel = L"(E_n - E_0) / \omega_c",
    title = "Optomechanics: lowest 8 levels",
)
colors = Makie.wong_colors()
for n in 2:8
    full = [E_full[i][n] - E_full[i][1] for i in eachindex(gs)]
    eff = [E_eff[i][n] - E_eff[i][1] for i in eachindex(gs)]
    scatter!(ax, collect(gs) ./ ωm_val, full; color = colors[n - 1], marker = :circle)
    lines!(ax, collect(gs) ./ ωm_val, eff; color = colors[n - 1], linestyle = :dash)
end
fig

# Dashed lines: Kerr-effective spectrum
# ``E_{n,m}^\mathrm{eff} = \omega_c n - (g^2/\omega_m)\,n^2 + \omega_m m``,
# visible as bundles of ``\omega_m``-spaced mechanical phonons attached to
# each Kerr-shifted cavity level.
# Markers: exact diagonalisation of the full radiation-pressure
# Hamiltonian.  They agree to truncation error across the whole coupling
# range — the polaron substitution is exact, and the package's eager
# canonicalisation reproduced the textbook derivation in three operator
# substitutions.

# ## Membrane-in-the-middle: quadratic coupling
#
# A different geometry — a dielectric membrane suspended *between* two
# fixed mirrors — replaces the linear ``a^\dagger a\,(b + b^\dagger)``
# coupling by a quadratic one,
#
# ```math
# H_\mathrm{mim} = \omega_c\, a^\dagger a + \omega_m\, b^\dagger b
#   + g_2\, a^\dagger a\,(b + b^\dagger)^2.
# ```
#
# Because the membrane sits at an optical node, the cavity frequency depends
# on the *square* of its displacement, enabling QND measurement of the
# phonon number.  Expanding ``(b + b^\dagger)^2`` is the package's job:

@variables g2
H_mim = ωc * a' * a + ωm * b' * b + g2 * a' * a * (b + b')^2

# Three distinct effects fall out: a constant cavity-frequency renormalisation
# ``+g_2``, a **photon-number-dependent mechanical frequency** ``2 g_2\,
# a^\dagger a`` accompanying ``b^\dagger b``, and a **photon-number-dependent
# two-phonon drive** ``g_2\, a^\dagger a\,(b^2 + b^{\dagger 2})``.  Restricting
# to a cavity Fock sector ``|n\rangle_c`` (which is preserved because
# ``a^\dagger a`` commutes with everything in ``H_\mathrm{mim}``), the
# mechanical Hamiltonian is
#
# ```math
# H_m(n) = (\omega_m + 2 g_2 n)\,b^\dagger b
#   + g_2 n\,(b^2 + b^{\dagger 2}) + (\omega_c + g_2)\,n
# ```
#
# — precisely the single-mode squeezing Hamiltonian solved in the
# [Bogoliubov example](bogoliubov.md), with effective oscillator frequency
# ``\omega(n) = \omega_m + 2 g_2 n`` and parametric drive
# ``\kappa(n) = 2 g_2 n``.  The Bogoliubov dispersion gives a
# **photon-number-conditional mechanical frequency**
#
# ```math
# \varepsilon(n) = \sqrt{\omega(n)^2 - \kappa(n)^2}
#   = \omega_m\,\sqrt{1 + 4\,g_2\,n / \omega_m}\,.
# ```
#
# The mechanical ground state in each cavity sector is a **squeezed vacuum**
# whose squeezing parameter grows monotonically with the cavity photon
# number.  Including the squeezing-induced zero-point shift, the full
# spectrum reads
#
# ```math
# E(n, m) = (\omega_c + g_2)\,n + \varepsilon(n)\,m
#   + \tfrac{1}{2}\bigl[\varepsilon(n) - \omega_m - 2 g_2 n\bigr].
# ```
#
# We verify by exact diagonalisation:

ωc_val2, ωm_val2, g2_val = 2.0, 1.0, 0.1
b_total2 = FockBasis(4) ⊗ FockBasis(40)
H_mim_op = dense(
    to_numeric(
        substitute(H_mim, Dict(ωc => ωc_val2, ωm => ωm_val2, g2 => g2_val)),
        b_total2,
    ),
)
F_mim = eigen(Hermitian(H_mim_op.data))

ε_mim(n) = ωm_val2 * sqrt(1 + 4 * g2_val * n / ωm_val2)
E_th(n, m) = (ωc_val2 + g2_val) * n + ε_mim(n) * m +
    (ε_mim(n) - ωm_val2 - 2 * g2_val * n) / 2

preds = sort(
    [(n, m, E_th(n, m)) for n in 0:3 for m in 0:10];
    by = t -> t[3],
)[1:12]
for k in 1:12
    n, m, eth = preds[k]
    enum = F_mim.values[k]
    println(
        "  (n=$n, m=$m)  E_th = $(round(eth; digits = 4))",
        "   E_num = $(round(enum; digits = 4))",
        "   |Δ| = $(round(abs(eth - enum); digits = 8))"
    )
end

# Every level matches the closed-form prediction to truncation precision.
# The first few photon-conditional mechanical frequencies are
# ``\varepsilon(0) = \omega_m``, ``\varepsilon(1) \approx 1.0954\,\omega_m``,
# ``\varepsilon(2) \approx 1.1832\,\omega_m`` — the **mechanical mode is
# stiffer in every higher-photon sector**, with the stiffening predictable
# in closed form.  This number-state-conditional spring constant is what
# makes the membrane-in-the-middle a QND phonon-counter: by reading the
# cavity transmission's frequency-tracked sidebands one literally counts
# mechanical quanta.
