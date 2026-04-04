using Symbolics: @variables

function SW(H0, V, S)
    return H0 + (1 // 2) * commutator(S, V)
end
function benchmark_Schrieffer_Wolf!(SUITE)
    # 1. Define Hilbert spaces
    h_boson = FockSpace(:cavity)
    h_spin = NLevelSpace(:atom, 2, 1)
    h = h_boson ⊗ h_spin

    # 2. Define operators
    @qnumbers a::Destroy(h)
    σge = Transition(h, :σ, 1, 2, 2)    # |g⟩⟨e|
    σee = Transition(h, :σ, 2, 2, 2)    # |e⟩⟨e|

    # 3. Hamiltonian
    @variables ω0 ε g
    H0 = ω0 * a' * a - ε * σee
    V = g * (a' * σge' + a * σge)

    # https://arxiv.org/pdf/2004.06534
    S = g / (ω0 + ε) * (a' * σge' - a * σge)

    # 5. Compute H_eff = H0 + 1/2 [S, V]
    SUITE["Jaynes–Cummings"]["First order Schrieffer-Wolf"] = @benchmarkable SW($H0, $V, $S) seconds =
        10
    return nothing
end
