using Symbolics: @variables

function benchmark_accumulation!(SUITE)
    # Additive reductions: `sum`/`reduce(+, …)` over a `Vector{QAdd}` accumulate
    # in place via the MutableArithmetics interface, vs the O(n²) `foldl(+, …)`
    # baseline that copies the growing accumulator on every step.

    # ---- Many-mode Hamiltonian: H = Σ_k ω_k a_k† a_k (distinct sites) ----
    for M in (8, 16, 24)
        hs = ⊗([FockSpace(Symbol(:m, k)) for k in 1:M]...)
        terms = [Float64(k) * (Destroy(hs, Symbol(:a, k), k)' * Destroy(hs, Symbol(:a, k), k)) for k in 1:M]
        SUITE["Accumulation"]["Many-mode H"]["foldl M=$M"] = @benchmarkable foldl(+, $terms) seconds = 3
        SUITE["Accumulation"]["Many-mode H"]["sum M=$M"] = @benchmarkable sum($terms) seconds = 3
    end

    # ---- Repeated same-site term: Σ c†c (collects to one term) ----
    hf = FockSpace(:f)
    c = Destroy(hf, :c)
    same_site = [c' * c for _ in 1:24]
    SUITE["Accumulation"]["Same-site"]["foldl"] = @benchmarkable foldl(+, $same_site) seconds = 3
    SUITE["Accumulation"]["Same-site"]["sum"] = @benchmarkable sum($same_site) seconds = 3

    return nothing
end
