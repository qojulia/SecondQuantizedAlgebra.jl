using SecondQuantizedAlgebra, QuantumOpticsBase

function benchmark_numeric_conversion!(SUITE)
    hfock = FockSpace(:cavity)
    hnlevel = NLevelSpace(:ThreeLevelAtom, (:a, :b, :c))
    h = hfock ⊗ hnlevel
    a = Destroy(h, :a)
    s = Transition(h, :s, :a, :c)
    levelmap = Dict(:a => 3, :b => 2, :c => 1)

    bfock = FockBasis(10)
    bnlevel = NLevelBasis(3)
    psi =
        coherentstate(bfock, 0.3) ⊗ (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) /
        sqrt(2)

    avg = average(a' * s)

    SUITE["Three Level Atom"]["Numeric conversion"] = @benchmarkable numeric_average(
        $avg, $psi; level_map=($levelmap)
    ) seconds = 10
end
