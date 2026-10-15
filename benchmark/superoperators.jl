import SecondQuantizedAlgebra: SuperOperator, liouvillian, to_numeric_master
import QuantumOpticsBase as QOB

function benchmark_superoperators!(SUITE)
    M = 4
    h = ProductSpace(ntuple(k -> FockSpace(Symbol(:c, k)), M))
    a = [Destroy(h, Symbol(:a, k), k) for k in 1:M]
    @variables Δ κ J
    H = sum([Δ * a[k]' * a[k] for k in 1:M]) + sum([J * (a[k]' * a[k + 1] + a[k + 1]' * a[k]) for k in 1:(M - 1)])
    L = liouvillian(H, a; rates = fill(κ, M))
    S = SuperOperator(L)
    ρ = a[1]' * a[2] + a[3] * a[4]'

    ha = NLevelSpace(:atom, 2)
    @variables N
    i = Index(ha, :i, N, ha)
    j = Index(ha, :j, N, ha)
    σ12 = IndexedOperator(Transition(ha, :σ, 1, 2), i)
    Γ = DoubleIndexedVariable(:Γ, i, j)

    group = SUITE["Superoperators"]
    group["Construction"]["Lindbladian build, 4 modes"] = @benchmarkable liouvillian($H, $a; rates = $(fill(κ, M))) seconds = 3 evals = 1
    group["Construction"]["Lowering to SuperOperator"] = @benchmarkable SuperOperator($L) seconds = 3 evals = 1
    group["Symbolic action"]["Composition L ∘ L"] = @benchmarkable $L ∘ $L seconds = 3 evals = 1
    group["Symbolic action"]["Apply L(ρ)"] = @benchmarkable $L($ρ) seconds = 3 evals = 1
    group["Symbolic action"]["Heisenberg drift"] = @benchmarkable transpose($S)($ρ) seconds = 3 evals = 1
    group["Construction"]["Collective dissipator"] = @benchmarkable dissipator($σ12, $Γ, $i => $j) seconds = 3 evals = 1

    # Numeric conversion: 4 modes with cutoff 3 give a 256-dimensional Hilbert space.
    par = Dict(Δ => 0.5, κ => 0.2, J => 0.1)
    b = QOB.FockBasis(3)^M
    hc = FockSpace(:c)
    c = Destroy(hc, :c)
    Lc = hamiltonian_action(c' * c) + dissipator(c)
    group["Numeric conversion"]["to_numeric(L), 4 modes"] = @benchmarkable to_numeric($L, $b; parameter = $par) seconds = 3 evals = 1
    group["Numeric conversion"]["to_numeric(L), single mode N=20"] = @benchmarkable to_numeric($Lc, $(QOB.FockBasis(20))) seconds = 3 evals = 1
    group["Numeric conversion"]["to_numeric_master, 4 modes"] = @benchmarkable to_numeric_master(
        $L, $h, $(ntuple(_ -> 3, M)); backend = QuantumOpticsBackend(), parameter = $par,
    ) seconds = 3 evals = 1
    return SUITE
end
