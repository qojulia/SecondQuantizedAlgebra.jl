using Symbolics: @variables

function benchmark_simplify_and_normal_order!(SUITE)

    # ==== SIMPLIFY ====

    # ---- Jaynes-Cummings scaling ----
    hf = FockSpace(:c)
    ha = NLevelSpace(:a, 2, 1)
    hjc = hf ⊗ ha
    a = Destroy(hjc, :a, 1)
    σge = Transition(hjc, :σ, 1, 2, 2)
    σee = Transition(hjc, :σ, 2, 2, 2)

    @variables ω_c ω_a g_jc
    H_jc = ω_c * a' * a + ω_a * σee + g_jc * (a' * σge + a * σge')

    SUITE["Simplify"]["Jaynes-Cummings"]["H"] = @benchmarkable simplify($H_jc) seconds = 3
    SUITE["Simplify"]["Jaynes-Cummings"]["H²"] = @benchmarkable simplify(
        $(H_jc * H_jc)
    ) seconds = 3

    # ---- Λ-system scaling ----
    h3 = FockSpace(:c) ⊗ NLevelSpace(:Λ, 3, 1)
    a3 = Destroy(h3, :a, 1)
    σ(i, j) = Transition(h3, :σ, i, j, 2)

    @variables Ω_p Ω_c Δ_p Δ_c
    H_lambda = Δ_p * σ(2, 2) + (Δ_p - Δ_c) * σ(3, 3) +
               Ω_p * (a3' * σ(1, 2) + a3 * σ(2, 1)) +
               Ω_c * (σ(2, 3) + σ(3, 2))

    SUITE["Simplify"]["Λ-system"]["H"] = @benchmarkable simplify($H_lambda) seconds = 3
    SUITE["Simplify"]["Λ-system"]["H²"] = @benchmarkable simplify(
        $(H_lambda * H_lambda)
    ) seconds = 3

    # ---- Two coupled cavities scaling ----
    h2c = FockSpace(:a) ⊗ FockSpace(:b)
    a1 = Destroy(h2c, :a, 1)
    a2 = Destroy(h2c, :b, 2)

    @variables J ω1 ω2
    H_2cav = ω1 * a1' * a1 + ω2 * a2' * a2 + J * (a1' * a2 + a2' * a1)

    SUITE["Simplify"]["Two cavities"]["H"] = @benchmarkable simplify($H_2cav) seconds = 3
    SUITE["Simplify"]["Two cavities"]["H²"] = @benchmarkable simplify(
        $(H_2cav * H_2cav)
    ) seconds = 3

    # ==== NORMAL ORDER ====

    # ---- Fock chain scaling ----
    hfock = FockSpace(:f)
    @qnumbers c::Destroy(hfock)

    for n in 2:5
        expr = (c * c')^n
        SUITE["Normal Order"]["Fock (c*c')^n"]["n=$n"] = @benchmarkable normal_order(
            $expr
        ) seconds = 3
    end

    # ---- Multi-mode ----
    two_mode = a1 * a2 * a1' * a2' * a1 * a2'
    SUITE["Normal Order"]["Multi-mode"]["2-mode 6-op chain"] = @benchmarkable normal_order(
        $two_mode
    ) seconds = 3

    # ---- Ground state rewriting ----
    hn3 = NLevelSpace(:atom, 3, 1)
    t(i, j) = Transition(hn3, :σ, i, j)
    gs_expr = t(1, 2) * t(2, 1) + t(1, 3) * t(3, 1) + t(1, 1) * t(2, 3)
    SUITE["Normal Order"]["Ground state"]["3-level rewrite"] = @benchmarkable normal_order(
        $gs_expr, $hn3
    ) seconds = 3

    return nothing
end
