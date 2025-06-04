function SW(H0, V, S)
    H0 + (1//2) * commutator(S, V)
end
function benchmark_Schrieffer_Wolf!(SUITE)
    # 1. Define Hilbert spaces
    h_boson = FockSpace(:cavity)
    h_spin = NLevelSpace(:atom, (:g, :e))
    h = h_boson ⊗ h_spin

    # 2. Define operators
    @qnumbers a::Destroy(h)
    σp = Transition(h, :σ, :g, :e)
    σz = Transition(h, :σ, :e, :e)

    # 3. Hamiltonian
    @cnumbers ω0 ε g
    H0 = ω0 * a' * a - ε * σz
    V = g * (a' * σp' + a * σp)
    H = H0 + V

    # https://arxiv.org/pdf/2004.06534
    S = g/(ω0 + ε) * (a'*σp' - a*σp)

    # 5. Compute H_eff = H0 + 1/2 [S, V]
    SUITE["Jaynes–Cummings"]["First order Schrieffer-Wolf"] = @benchmarkable SW($H0, $V, $S) seconds =
        10
    return nothing
end
