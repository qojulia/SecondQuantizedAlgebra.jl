# # Two-Atom Dipole-Dipole: Bright and Dark States
#
# Two identical two-level atoms in close proximity exchange virtual photons,
# giving rise to a dipole-dipole interaction that splits the single-excitation
# manifold into a symmetric (bright, superradiant) and antisymmetric (dark,
# subradiant) eigenstate.  This is the elementary building block of
# *collective* atomic physics and the simplest non-trivial setting for the
# package's index machinery: two indexed atomic operators on the same
# Hilbert subspace, kept distinct via [`assume_distinct_index`](@ref).

# ## Setup
#
# A single two-level Hilbert space, populated by two atoms via abstract
# indices ``i`` and ``j``.  The convention `(:g, :e)` lets us name
# transitions by symbol rather than by integer level.

using SecondQuantizedAlgebra

ha = NLevelSpace(:atom, (:g, :e))

i = Index(ha, :i, 2, ha)
j = Index(ha, :j, 2, ha)

σge(idx) = IndexedOperator(Transition(ha, :σ, :g, :e), idx)
σeg(idx) = IndexedOperator(Transition(ha, :σ, :e, :g), idx)
σee(idx) = IndexedOperator(Transition(ha, :σ, :e, :e), idx)

@variables ω J

# ## Hamiltonian
#
# Two identical atoms at frequency ``\omega`` exchange excitations through a
# coherent dipole-dipole coupling ``J``:
#
# ```math
# H = \omega\,(\sigma_i^{ee} + \sigma_j^{ee})
#   + J\,(\sigma_i^{eg} \sigma_j^{ge} + \sigma_j^{eg} \sigma_i^{ge}).
# ```
#
# Building it with free indices leaves the site relationship `i = j` vs
# `i ≠ j` ambiguous.  Calling [`assume_distinct_index`](@ref) records the
# pair `(i, j)` as a permanent inequality constraint and triggers
# re-canonicalisation, including completeness expansion of any ground-state
# projectors that emerge from the new constraint:

H = ω * (σee(i) + σee(j)) + J * (σeg(i) * σge(j) + σeg(j) * σge(i))
H = assume_distinct_index(H, [(i, j)])

# ## Symmetric and antisymmetric raising operators
#
# Define
#
# ```math
# S_\pm^\dagger
#   = \frac{1}{\sqrt{2}}\,(\sigma_i^{eg} \pm \sigma_j^{eg}),
# ```
#
# so ``S_+^\dagger\, |gg\rangle`` is the bright single-excitation state and
# ``S_-^\dagger\, |gg\rangle`` is the dark one.  In the package we drop the
# ``1/\sqrt{2}`` to keep symbolic expressions clean — only the eigenvalues
# matter.

S_bright = assume_distinct_index(σeg(i) + σeg(j), [(i, j)])
S_dark = assume_distinct_index(σeg(i) - σeg(j), [(i, j)])

# ## Diagonalising the single-excitation manifold
#
# Acting on the ground state ``|gg\rangle`` we have ``H\,|gg\rangle = 0``,
# so ``H\,S_\pm^\dagger\,|gg\rangle = [H, S_\pm^\dagger]\,|gg\rangle``.  Each
# commutator is a single line in the package, but the raw form mixes
# single-excitation and two-excitation operators because ``\sigma^{ee}`` and
# ``\sigma^{eg}`` don't quite commute even at distinct sites.
# [`expand_completeness`](@ref) projects to the ``\{\sigma^{ij} : (i,j) \neq
# (g,g)\} \cup \{1\}`` basis and the structure becomes clear:

expand_completeness(commutator(H, S_bright))

#-

expand_completeness(commutator(H, S_dark))

# Reading off the single-excitation parts,
#
# ```math
# [H, S_+^\dagger] = (\omega + J)\,S_+^\dagger
#   - 2J\,(\sigma_i^{ee}\sigma_j^{eg} + \sigma_i^{eg}\sigma_j^{ee}), \qquad
# [H, S_-^\dagger] = (\omega - J)\,S_-^\dagger
#   - 2J\,(\sigma_i^{ee}\sigma_j^{eg} - \sigma_i^{eg}\sigma_j^{ee}).
# ```
#
# The residual ``\sigma^{ee}\sigma^{eg}`` terms vanish identically on the
# ground state (they require one site already excited).  Restricted to the
# single-excitation manifold:
#
# ```math
# H\,S_\pm^\dagger\,|gg\rangle = (\omega \pm J)\,S_\pm^\dagger\,|gg\rangle.
# ```
#
# The **bright** state ``S_+^\dagger\,|gg\rangle`` sits at ``\omega + J`` and
# carries the full collective dipole moment — when the atoms are coupled to a
# common radiation continuum, this state decays at rate ``2\Gamma``
# (superradiance).  The **dark** state ``S_-^\dagger\,|gg\rangle`` sits at
# ``\omega - J`` and has zero collective dipole moment — it is decoupled
# from the symmetric radiation mode and is effectively immortal
# (subradiance).

# ## Number conservation
#
# Total atomic excitation ``N_e = \sigma_i^{ee} + \sigma_j^{ee}`` commutes
# with ``H``:

commutator(assume_distinct_index(σee(i) + σee(j), [(i, j)]), H)

# So the four-dimensional Hilbert space (``2\times 2`` atoms = ``\{|gg\rangle,
# |eg\rangle, |ge\rangle, |ee\rangle\}``) decomposes into the
# zero-excitation singlet, the dipole-split single-excitation doublet
# ``\{|+\rangle, |-\rangle\}`` we just diagonalised, and the
# two-excitation singlet ``|ee\rangle`` at energy ``2\omega``.

# ## Numerical verification
#
# We diagonalise the 4-dimensional Hilbert space and check the spectrum
# directly:

using QuantumOpticsBase, LinearAlgebra

ω_val, J_val = 1.0, 0.25
b_atom = NLevelBasis(2)
b = b_atom ⊗ b_atom

## Build numeric operators by hand (two atoms on the same NLevelSpace are
## the simplest non-indexed scenario for to_numeric).
σ_e_a = (transition(b_atom, 2, 2)) ⊗ one(b_atom)
σ_e_b = one(b_atom) ⊗ transition(b_atom, 2, 2)
σ_eg_a = transition(b_atom, 2, 1) ⊗ one(b_atom)
σ_ge_a = transition(b_atom, 1, 2) ⊗ one(b_atom)
σ_eg_b = one(b_atom) ⊗ transition(b_atom, 2, 1)
σ_ge_b = one(b_atom) ⊗ transition(b_atom, 1, 2)

H_num = ω_val * (σ_e_a + σ_e_b) + J_val * (σ_eg_a * σ_ge_b + σ_eg_b * σ_ge_a)
E = sort(real.(eigvals(Matrix(H_num.data))))

println("spectrum             = ", round.(E; digits = 4))
println("expected {0, ω-J, ω+J, 2ω} = ",
        round.([0.0, ω_val - J_val, ω_val + J_val, 2*ω_val]; digits = 4))

# The four eigenvalues come out as ``\{0, \omega - J, \omega + J,
# 2\omega\}`` — the singlet, dark, bright, and doubly-excited levels — in
# exact agreement with the symbolic derivation.  All four sat in a
# four-dimensional Hilbert space derived from two indexed atoms with no
# physical Hilbert-space duplication: the same `NLevelSpace` carries both
# sites, with the index machinery handling site identity.

# ## Coupling to a cavity: vacuum Rabi splitting
#
# The bright-dark dichotomy gets its name from a *physical* observable:
# coupling to a radiation mode.  Adding a common cavity field at frequency
# ``\omega_c`` to the existing two-atom system,
#
# ```math
# H \;\to\; H + \omega_c\, a^\dagger a
#   + g\,a^\dagger\,(\sigma_i^{ge} + \sigma_j^{ge})
#   + g\,a\,(\sigma_i^{eg} + \sigma_j^{eg}),
# ```
#
# the atom-field coupling appears only through the **collective lowering**
# ``S_- = \sigma_i^{ge} + \sigma_j^{ge}`` and its hermitian conjugate.  The
# dark raising operator ``S_-^\dagger`` is orthogonal to ``S_+`` and the
# package confirms it directly — the cavity field commutes with the dark
# operator:

hcav = FockSpace(:cavity)
hatom = NLevelSpace(:atom, (:g, :e))
hh = hcav ⊗ hatom
@variables ωc g
@qnumbers a::Destroy(hh, 1)
ic = Index(hh, :i, 2, hatom)
jc = Index(hh, :j, 2, hatom)
σge_c(idx) = IndexedOperator(Transition(hatom, :σ, :g, :e), idx)
σeg_c(idx) = IndexedOperator(Transition(hatom, :σ, :e, :g), idx)

H_int = g * a' * (σge_c(ic) + σge_c(jc)) + g * a * (σeg_c(ic) + σeg_c(jc))
H_int = assume_distinct_index(H_int, [(ic, jc)])

S_dark_c = assume_distinct_index(σeg_c(ic) - σeg_c(jc), [(ic, jc)])

commutator(H_int, S_dark_c)

# The output is supported only on two-excitation operators of the form
# ``\sigma^{ee}\sigma^{eg}`` that **annihilate the ground state**.  So
# acting on ``|0_c, gg\rangle``, ``[H_\mathrm{int}, S_-^\dagger]`` returns
# zero — the dark single-excitation state is invisible to the cavity.
#
# The bright state, by contrast, hybridises with the cavity photon to form
# polaritons.  In the single-excitation basis
# ``\{|1_c, gg\rangle, |0_c, +\rangle, |0_c, -\rangle\}``, the package's
# canonical form factorises into a ``2\times 2`` polariton block and a
# decoupled dark eigenvalue.  At resonance
# ``\omega_c = \omega + J`` the polariton block is
# ``\begin{pmatrix}\omega_c & g\sqrt{2}\\ g\sqrt{2} & \omega + J\end{pmatrix}``,
# giving the **vacuum Rabi splitting** ``2g\sqrt{2}``:

ωc_val, g_val = ω_val + J_val, 0.2
b_cav = FockBasis(2)
b_full = b_cav ⊗ b_atom ⊗ b_atom
a_num = destroy(b_cav) ⊗ one(b_atom) ⊗ one(b_atom)
σ_ee_a2 = one(b_cav) ⊗ transition(b_atom, 2, 2) ⊗ one(b_atom)
σ_ee_b2 = one(b_cav) ⊗ one(b_atom) ⊗ transition(b_atom, 2, 2)
σ_eg_a2 = one(b_cav) ⊗ transition(b_atom, 2, 1) ⊗ one(b_atom)
σ_ge_a2 = one(b_cav) ⊗ transition(b_atom, 1, 2) ⊗ one(b_atom)
σ_eg_b2 = one(b_cav) ⊗ one(b_atom) ⊗ transition(b_atom, 2, 1)
σ_ge_b2 = one(b_cav) ⊗ one(b_atom) ⊗ transition(b_atom, 1, 2)

H_tc = ωc_val * a_num' * a_num + ω_val * (σ_ee_a2 + σ_ee_b2) +
       J_val * (σ_eg_a2 * σ_ge_b2 + σ_eg_b2 * σ_ge_a2) +
       g_val * (a_num' * σ_ge_a2 + a_num' * σ_ge_b2 +
                a_num * σ_eg_a2 + a_num * σ_eg_b2)

E_tc = sort(real.(eigvals(Hermitian(Matrix(H_tc.data)))))
pol_lower = ωc_val - g_val * sqrt(2)
pol_upper = ωc_val + g_val * sqrt(2)
dark_lvl = ω_val - J_val
println("single-excitation spectrum (numeric): ",
        round.(E_tc[2:4]; digits = 4))
println("polariton + dark prediction:           ",
        round.(sort([pol_lower, dark_lvl, pol_upper]); digits = 4))

# The three single-excitation eigenvalues fall right where the analytic
# story predicts: a lower polariton at ``\omega_c - g\sqrt{2}``, the
# untouched dark state at ``\omega - J``, and an upper polariton at
# ``\omega_c + g\sqrt{2}``.  Dropping a single atom from this experiment
# would leave only the ordinary Jaynes-Cummings ``g`` splitting — the
# ``\sqrt{2}`` is the collective enhancement.

# ## Generalising to N atoms with `Σ`
#
# Stepping up to ``N`` identical atoms uses the same operator
# constructions, but now with a *symbolic* range ``N`` and explicit
# [`Σ`](@ref) sums.  The dipole-dipole term becomes
# ``J\sum_{i \neq j} \sigma_i^{eg}\sigma_j^{ge}`` and the cavity coupling
# becomes ``g\sum_i (a^\dagger \sigma_i^{ge} + a\,\sigma_i^{eg})``:

@variables Nsym
iN = Index(hh, :i, Nsym, hatom)
jN = Index(hh, :j, Nsym, hatom)

H_dd_N = J * Σ(Σ(σeg_c(iN) * σge_c(jN), jN, [iN]), iN)
H_int_N = g * Σ(a' * σge_c(iN) + a * σeg_c(iN), iN)

H_dd_N

#-

H_int_N

# The single-excitation block decomposes exactly as for ``N = 2``: a
# ``2\times 2`` polariton sector for ``\{|1_c, gg\cdots g\rangle,
# |0_c, B\rangle\}`` with bright energy ``\omega + (N-1)J`` and enhanced
# coupling ``g\sqrt{N}``, plus an ``(N-1)``-fold degenerate dark manifold
# at ``\omega - J``.  At resonance ``\omega_c = \omega + (N-1)J`` the
# polariton splitting is the **collective vacuum Rabi splitting**
# ``2g\sqrt{N}`` — a hallmark of cavity QED with cold-atom ensembles.

function tavis_single_excitation(N::Int; ω = 1.0, J = 0.05, g = 0.2)
    b_c = FockBasis(2)
    b_a = NLevelBasis(2)
    bs = vcat([b_c], fill(b_a, N))
    b_tot = reduce(⊗, bs)
    op_at(idx, op_local) = reduce(
        ⊗,
        vcat(
            [idx == 0 ? op_local : one(b_c)],
            [k == idx ? op_local : one(b_a) for k in 1:N],
        ),
    )
    a_n = op_at(0, destroy(b_c))
    σee_n = [op_at(k, transition(b_a, 2, 2)) for k in 1:N]
    σeg_n = [op_at(k, transition(b_a, 2, 1)) for k in 1:N]
    σge_n = [op_at(k, transition(b_a, 1, 2)) for k in 1:N]
    ωc = ω + (N - 1) * J
    H_N = ωc * a_n' * a_n + ω * sum(σee_n) +
          J * sum(σeg_n[p] * σge_n[q] for p in 1:N, q in 1:N if p != q) +
          g * sum(a_n' * σge_n[k] + a_n * σeg_n[k] for k in 1:N)
    N_tot = a_n' * a_n + sum(σee_n)
    Hd = Hermitian(Matrix(H_N.data))
    F = eigen(Hd)
    Ns = real.(diag(F.vectors' * Matrix(N_tot.data) * F.vectors))
    keep = findall(n -> abs(n - 1) < 1.0e-6, Ns)
    return ωc, sort(real.(F.values[keep]))
end

for N in 2:4
    ωc_val_N, E_single = tavis_single_excitation(N)
    pol_lo_th = ωc_val_N - 0.2 * sqrt(N)
    pol_hi_th = ωc_val_N + 0.2 * sqrt(N)
    dark_th = 1.0 - 0.05
    th_sorted = sort(vcat([pol_lo_th, pol_hi_th], fill(dark_th, N - 1)))
    println("N = $N:  numeric single-excitation = ",
            round.(E_single; digits = 4))
    println("        theory                     = ",
            round.(th_sorted; digits = 4),
            "   (ε_pol± = ω_c ± g√N, ",
            "dark ×$(N - 1) at ω - J)")
end

# The lower polariton tracks ``\omega_c - g\sqrt{N}`` exactly for every
# ``N``, demonstrating the collective ``\sqrt{N}`` enhancement in a closed
# form.  Doubling the atom count gives roughly a ``\sqrt{2}`` increase in
# the bright-state cavity coupling — invisible to the dark manifold, which
# remains at ``\omega - J`` regardless of how many atoms join.  All of
# this came from the same indexed `NLevelSpace`, a single `\Sigma` for the
# atom-cavity coupling, and a doubly-nested `\Sigma` (with `i \neq j`) for
# the dipole-dipole term — no per-atom Hilbert-space book-keeping.
