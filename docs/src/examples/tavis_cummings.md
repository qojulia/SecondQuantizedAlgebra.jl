```@meta
EditURL = "../../../examples/tavis_cummings.jl"
```

# Exact finite-N Tavis–Cummings model in the symmetric subspace

`CollectiveTransition` keeps the atomic Hilbert space polynomial in atom
number while retaining exact finite-N dynamics. For two-level atoms the
symmetric sector has dimension `N + 1`, instead of `2^N` for the full product
space.

````@example tavis_cummings
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
Sgg = CollectiveTransition(h, :S, :g, :g, 2)

@variables ωc ωa g
H = ωc * a' * a + ωa * See + g * (a' * Sge + a * Seg)
````

The symbolic algebra applies the collective SU(2) commutator.

````@example tavis_cummings
iszero(simplify(commutator(Sge, Seg) - Sgg + See))
````

Fixing total atomic occupation to N constructs the symmetric sector.

````@example tavis_cummings
bc = FockBasis(nphotons)
b1 = NLevelBasis(2)
ba = ManyBodyBasis(b1, bosonstates(b1, N))
b = bc ⊗ ba
````

Substitute physical parameters and materialize the exact sparse Hamiltonian.

````@example tavis_cummings
Hnum = to_numeric(H, b; parameter = Dict(ωc => 1.0, ωa => 1.0, g => 0.1))

@assert length(ba) == N + 1
@assert size(Hnum.data) == (length(b), length(b))
````

Diagonal collective populations sum to N, not one.

````@example tavis_cummings
number_identity = to_numeric(Sgg + See, b)
@assert number_identity ≈ N * one(b)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

