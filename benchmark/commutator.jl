using Symbolics: @variables

function nested_commutator(H, op, n)
    if n == 0
        return op
    else
        return commutator(H, nested_commutator(H, op, n - 1))
    end
end

function benchmark_commutator!(SUITE)
    # ---- Jaynes-Cummings nested ----
    hf = FockSpace(:c)
    ha = NLevelSpace(:a, 2, 1)
    hjc = hf ⊗ ha
    a = Destroy(hjc, :a, 1)
    σge = Transition(hjc, :σ, 1, 2, 2)
    σee = Transition(hjc, :σ, 2, 2, 2)

    @variables ω_c ω_a g
    H_jc = ω_c * a' * a + ω_a * σee + g * (a' * σge + a * σge')

    for n in 1:5
        SUITE["Commutator"]["Nested JC"]["depth=$n"] = @benchmarkable nested_commutator(
            $H_jc, $σge, $n
        ) seconds = 3
    end

    # ---- Schrieffer-Wolff: [S, [S, H0]] ----
    @variables ε g_sw
    H0 = ω_c * a' * a - ε * σee
    V = g_sw * (a' * σge' + a * σge)
    S = g_sw / (ω_c + ε) * (a' * σge' - a * σge)

    SUITE["Commutator"]["Schrieffer-Wolff"]["[S, V]"] = @benchmarkable commutator(
        $S, $V
    ) seconds = 3
    SUITE["Commutator"]["Schrieffer-Wolff"]["[S, [S, H0]]"] = @benchmarkable begin
        c1 = commutator($S, $H0)
        commutator($S, c1)
    end seconds = 3

    return nothing
end
