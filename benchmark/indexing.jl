using Symbolics: @variables

function benchmark_indexing!(SUITE)
    hf = FockSpace(:cavity)
    hn = NLevelSpace(:atom, 2, 1)
    h = hf ⊗ hn
    a = Destroy(h, :a, 1)

    @variables N ω_c ω_a
    i = Index(h, :i, N, hn)
    j = Index(h, :j, N, hn)
    gi = IndexedVariable(:g, i)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

    H_jc = ω_c * a' * a + Σ(ω_a * σ(2, 2, i) + gi * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    # ---- Diagonal collapse: [Σ_i H_i, O_j] ----
    SUITE["Indexing"]["Diagonal collapse"]["[H_JC, σ_j]"] = @benchmarkable commutator(
        $H_jc, $(σ(2, 2, j))
    ) seconds = 3

    # Dicke model with spins
    hd_f = FockSpace(:c)
    hd_s = SpinSpace(:s)
    hd = hd_f ⊗ hd_s
    b = Destroy(hd, :b, 1)
    id = Index(hd, :i, N, hd_s)
    jd = Index(hd, :j, N, hd_s)
    Sx(idx) = IndexedOperator(Spin(hd, :S, 1, 2), idx)
    Sz(idx) = IndexedOperator(Spin(hd, :S, 3, 2), idx)

    @variables ω λ
    H_dicke = ω * b' * b + Σ(ω * Sz(id) + λ * (b + b') * Sx(id), id)

    SUITE["Indexing"]["Diagonal collapse"]["[H_Dicke, S_j]"] = @benchmarkable commutator(
        $H_dicke, $(Sz(jd))
    ) seconds = 3

    # ---- Expand sums ----
    single_sum = Σ(gi * σ(2, 1, i) * σ(1, 2, j), i)
    SUITE["Indexing"]["Expand sums"]["single Σ_i(σ_i*σ_j)"] = @benchmarkable expand_sums(
        $single_sum
    ) seconds = 3

    Jij = DoubleIndexedVariable(:J, id, jd; identical = false)
    double_sum = Σ(Σ(Jij * Sx(id) * Sx(jd), id), jd)
    SUITE["Indexing"]["Expand sums"]["double Σ_ij(J_ij*S_i*S_j)"] = @benchmarkable expand_sums(
        $double_sum
    ) seconds = 3

    # ---- Simplify indexed expressions ----
    SUITE["Indexing"]["Simplify"]["indexed JC H"] = @benchmarkable simplify(
        $H_jc
    ) seconds = 3

    SUITE["Indexing"]["Simplify"]["double-sum spin-spin"] = @benchmarkable simplify(
        $double_sum
    ) seconds = 3

    return nothing
end
