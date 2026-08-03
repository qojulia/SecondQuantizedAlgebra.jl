# # Cavity Optomechanics and the Polaron Transformation
#
# A single optical cavity coupled to a mechanical oscillator by radiation
# pressure has the Hamiltonian
#
# ```math
# H = \omega_0\, a^\dagger a + \omega_m\, b^\dagger b
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
# ``U^\dagger b U = b - (g/\omega_m)\,a^\dagger a`` is exact, since the Hadamard
# series truncates at first order.  The displacement amplitude is an operator,
# but a conserved one, so [`Displace`](@ref) reproduces the transformation
# photon-number sector by photon-number sector.

# ## Setup

using SecondQuantizedAlgebra

hc = FockSpace(:cavity)
hm = FockSpace(:mech)
h = hc ⊗ hm

@qnumbers a::Destroy(h, 1)
@qnumbers b::Destroy(h, 2)

@variables ω₀ ωₘ g

H = ω₀ * a' * a + ωₘ * b' * b + g * a' * a * (b + b')

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
# \dot a = -i\,\omega_0\, a - i\, g\, a\,(b + b^\dagger).
# ```
#
# The first equation makes the **displacement** structure explicit: at fixed
# ``a^\dagger a = n``, the mechanical mode oscillates around a shifted origin
# ``\langle b\rangle_\text{eq} = -(g/\omega_m)\,n``.

# ## Polaron transformation
#
# A unitary transformation acts on operators by conjugation, ``\tilde O =
# U^\dagger O U``.  The polaron unitary displaces the mechanical mode by an
# amount that is itself an operator, the conserved cavity occupation:
#
# ```math
# \tilde b = b - \tfrac{g}{\omega_m}\,a^\dagger a .
# ```
#
# [`Displace`](@ref) takes that amplitude directly.  It must commute with the
# mode it displaces and with its own adjoint, which ``a^\dagger a`` does, and
# then the Hadamard series truncates at first order exactly as it would for a
# number.

U_pol = Displace(b, (-g / ωₘ) * (a' * a))

conjugate(b, U_pol)

# [`is_canonical`](@ref) certifies it against the canonical commutator of the
# site it acts on, with no matrix representation of ``U`` involved:

is_canonical(U_pol)

# ## Effective Kerr Hamiltonian
#
# Conjugating ``H`` eliminates the coupling outright:

H_pol = conjugate(H, U_pol)

# ```math
# \boxed{\;
#   H_\mathrm{pol}
#   = \omega_0\, a^\dagger a + \omega_m\, b^\dagger b
#   - \frac{g^2}{\omega_m}\,(a^\dagger a)^2,
# \;}
# ```
#
# written above with ``(a^\dagger a)^2 = a^\dagger{}^2 a^2 + a^\dagger a``, which
# is the normal-ordered form the algebra returns.  The optomechanical coupling
# has been **eliminated entirely** in favour of a Kerr nonlinearity of strength
# ``K = g^2/\omega_m`` on the cavity mode.  Photons attract each other with
# energy ``-K\,n(n-1)``, a self-Kerr anharmonicity that is now routinely
# measured in superconducting and membrane-in-the-middle optomechanical
# experiments.  The mechanical mode, in the polaron frame, has decoupled
# completely from the optics.
#
# The cavity is a different matter.  ``a`` does not commute with the amplitude,
# so ``U^\dagger a U`` is ``a`` times a displacement operator, which is not a
# polynomial in the generators at all.  Rather than return it untransformed,
# `conjugate` refuses:

try
    conjugate(a, U_pol)
catch err
    err
end

# Everything commuting with ``a^\dagger a`` does pass through, which is why the
# conjugation of ``H`` above went through untouched on the cavity side.

conjugate(a' * a, U_pol)

# ## Numerical verification
#
# We diagonalise both the original ``H`` and the polaron-frame ``H_\mathrm{pol}``
# computed above, and check that they share the same spectrum.  A frame change is
# a similarity transformation, so agreement here is a check on the truncation as
# much as on the algebra.

using QuantumOpticsBase, LinearAlgebra, CairoMakie

ω₀_val, ωm_val = 1.0, 0.4
n_max_c, n_max_m = 4, 60
b_cav = FockBasis(n_max_c)
b_mech = FockBasis(n_max_m)
b_total = b_cav ⊗ b_mech

gs = range(0.0, 0.6 * ωm_val, length = 16)
E_full = Vector{Float64}[]
E_eff = Vector{Float64}[]
for g_val in gs
    subs = Dict(ω₀ => ω₀_val, ωₘ => ωm_val, g => g_val)
    Hf = dense(to_numeric(substitute(H, subs), b_total))
    He = dense(to_numeric(substitute(H_pol, subs), b_total))
    push!(E_full, sort(real.(eigvals(Hermitian(Hf.data))))[1:8])
    push!(E_eff, sort(real.(eigvals(Hermitian(He.data))))[1:8])
end

fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel = L"g / \omega_m",
    ylabel = L"(E_n - E_0) / \omega_0",
    title = "Optomechanics: lowest 8 levels",
)
colors = Makie.wong_colors()
for k in 2:8
    full = [E_full[i][k] - E_full[i][1] for i in eachindex(gs)]
    eff = [E_eff[i][k] - E_eff[i][1] for i in eachindex(gs)]
    scatter!(ax, collect(gs) ./ ωm_val, full; color = colors[k - 1], marker = :circle)
    lines!(ax, collect(gs) ./ ωm_val, eff; color = colors[k - 1], linestyle = :dash)
end
fig

# Dashed lines: Kerr-effective spectrum
# ``E_{n,m}^\mathrm{eff} = \omega_0 n - (g^2/\omega_m)\,n^2 + \omega_m m``,
# visible as bundles of ``\omega_m``-spaced mechanical phonons attached to
# each Kerr-shifted cavity level.
# Markers: exact diagonalisation of the full radiation-pressure
# Hamiltonian.  They agree to truncation error across the whole coupling
# range: the polaron transformation is exact, and one `conjugate` call
# reproduced the textbook derivation.

# ## Membrane-in-the-middle: quadratic coupling
#
# A different geometry (a dielectric membrane suspended *between* two
# fixed mirrors) replaces the linear ``a^\dagger a\,(b + b^\dagger)``
# coupling by a quadratic one,
#
# ```math
# H_\mathrm{mim} = \omega_0\, a^\dagger a + \omega_m\, b^\dagger b
#   + g_2\, a^\dagger a\,(b + b^\dagger)^2.
# ```
#
# Because the membrane sits at an optical node, the cavity frequency depends
# on the *square* of its displacement, enabling QND measurement of the
# phonon number.  Expanding ``(b + b^\dagger)^2`` is the package's job:

@variables g₂
H_mim = ω₀ * a' * a + ωₘ * b' * b + g₂ * a' * a * (b + b')^2

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
#   + g_2 n\,(b^2 + b^{\dagger 2}) + (\omega_0 + g_2)\,n
# ```
#
# which is the single-mode squeezing Hamiltonian, of oscillator frequency
# ``\omega(n) = \omega_m + 2 g_2 n`` and parametric drive
# ``\kappa(n) = 2 g_2 n``.  Its two-mode cousin is the subject of the
# [Bogoliubov example](bogoliubov.md); here a single-mode [`Squeeze`](@ref)
# diagonalises it.  Written in a generic frequency ``\Omega`` and drive ``\kappa``:

@variables Ω κ r
conjugate(Ω * b' * b + κ / 2 * (b * b + b' * b'), Squeeze(b, r))

# The ``b^2`` coefficient is ``\tfrac{1}{2}(\Omega\sinh 2r + \kappa\cosh 2r)``
# and vanishes at ``\tanh 2r = -\kappa/\Omega``.  There the ``b^\dagger b``
# coefficient ``\Omega\cosh 2r + \kappa\sinh 2r`` collapses to
# ``\sqrt{\Omega^2 - \kappa^2}``, giving a **photon-number-conditional
# mechanical frequency**
#
# ```math
# \varepsilon(n) = \sqrt{\omega(n)^2 - \kappa(n)^2}
#   = \omega_m\,\sqrt{1 + 4\,g_2\,n / \omega_m}\,,
# ```
#
# and the leftover constant ``\tfrac{1}{2}[\varepsilon(n) - \omega(n)]`` is the
# zero-point shift of the squeezed vacuum.  The mechanical ground state in each
# cavity sector is therefore a **squeezed vacuum** whose squeezing parameter
# grows monotonically with the cavity photon number, and the full spectrum reads
#
# ```math
# E(n, m) = (\omega_0 + g_2)\,n + \varepsilon(n)\,m + \tfrac{1}{2}\bigl[\varepsilon(n) - \omega_m - 2 g_2 n\bigr].
# ```
#
# We verify by exact diagonalisation:

ω₀_val2, ωm_val2, g2_val = 2.0, 1.0, 0.1
b_total2 = FockBasis(4) ⊗ FockBasis(40)
H_mim_op = dense(
    to_numeric(
        substitute(H_mim, Dict(ω₀ => ω₀_val2, ωₘ => ωm_val2, g₂ => g2_val)),
        b_total2,
    ),
)
F_mim = eigen(Hermitian(H_mim_op.data))

ε_mim(n) = ωm_val2 * sqrt(1 + 4 * g2_val * n / ωm_val2)
E_th(n, m) = (ω₀_val2 + g2_val) * n + ε_mim(n) * m +
    (ε_mim(n) - ωm_val2 - 2 * g2_val * n) / 2

preds = sort(
    [(n, m, E_th(n, m)) for n in 0:3 for m in 0:10];
    by = t -> t[3],
)[1:12]
for k in 1:12
    nk, mk, eth = preds[k]
    enum = F_mim.values[k]
    println(
        "  (n=$nk, m=$mk)  E_th = $(round(eth; digits = 4))",
        "   E_num = $(round(enum; digits = 4))",
        "   |Δ| = $(round(abs(eth - enum); digits = 8))"
    )
end

# Every level matches the closed-form prediction to truncation precision.
# The first few photon-conditional mechanical frequencies are
# ``\varepsilon(0) = \omega_m``, ``\varepsilon(1) \approx 1.0954\,\omega_m``,
# ``\varepsilon(2) \approx 1.1832\,\omega_m``.  The **mechanical mode is
# stiffer in every higher-photon sector**, with the stiffening predictable
# in closed form.  This number-state-conditional spring constant is what
# makes the membrane-in-the-middle a QND phonon-counter: by reading the
# cavity transmission's frequency-tracked sidebands one literally counts
# mechanical quanta.
