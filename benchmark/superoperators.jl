import SecondQuantizedAlgebra: SuperOperator, liouvillian

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
    group["Lindbladian build, 4 modes"] = @benchmarkable liouvillian($H, $a; rates = $(fill(κ, M))) seconds = 3 evals = 1
    group["Lowering to SuperOperator"] = @benchmarkable SuperOperator($L) seconds = 3 evals = 1
    group["Composition L ∘ L"] = @benchmarkable $L ∘ $L seconds = 3 evals = 1
    group["Apply L(ρ)"] = @benchmarkable $L($ρ) seconds = 3 evals = 1
    group["Heisenberg drift"] = @benchmarkable transpose($S)($ρ) seconds = 3 evals = 1
    group["Collective dissipator"] = @benchmarkable dissipator($σ12, $Γ, $i => $j) seconds = 3 evals = 1
    return SUITE
end
