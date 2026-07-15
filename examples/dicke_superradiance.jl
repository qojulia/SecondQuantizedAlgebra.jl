# # Exact finite-N Tavis–Cummings model in the symmetric subspace
#
# `CollectiveTransition` keeps the atomic Hilbert space polynomial in atom
# number while retaining exact finite-N dynamics. For two-level atoms the
# symmetric sector has dimension `N + 1`, instead of `2^N` for the full product
# space.

using SecondQuantizedAlgebra
using QuantumOpticsBase

N = 20
nphotons = 12

hc = FockSpace(:cavity)
ha = CollectiveNLevelSpace(:atoms, (:g, :e))
h = hc ⊗ ha

@qnumbers a::Destroy(h, 1)
Sge = CollectiveTransition(h, :S, :g, :e, 2)
Seg = Sge'
See = CollectiveTransition(h, :S, :e, :e, 2)

@variables ωc ωa g
H = ωc * a' * a + ωa * See + g * (a' * Sge + a * Seg)

# The symbolic algebra knows the collective su(2) commutator exactly.
@assert iszero(
    simplify(
        commutator(Sge, Seg) -
            CollectiveTransition(h, :S, :g, :g, 2) + See
    )
)

# Fixing total atomic occupation to N constructs the symmetric sector.
bc = FockBasis(nphotons)
b1 = NLevelBasis(2)
ba = ManyBodyBasis(b1, bosonstates(b1, N))
b = bc ⊗ ba

# Substitute physical parameters and materialize the exact sparse Hamiltonian.
Hnum = to_numeric(H, b; parameter = Dict(ωc => 1.0, ωa => 1.0, g => 0.1))

@assert length(ba) == N + 1
@assert size(Hnum.data) == (length(b), length(b))

# Diagonal collective populations sum to N, not one.
Sgg = CollectiveTransition(h, :S, :g, :g, 2)
number_identity = to_numeric(Sgg + See, b)
@assert number_identity ≈ N * one(b)
