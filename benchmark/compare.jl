using BenchmarkTools
using SecondQuantizedAlgebra

# --- Setup ---

# Fock
hf = FockSpace(:fock)
@qnumbers a::Destroy(hf)

# NLevel
hn = NLevelSpace(:atom, 3, 1)
σ12 = Transition(hn, :σ, 1, 2)
σ13 = Transition(hn, :σ, 1, 3)
σ23 = Transition(hn, :σ, 2, 3)

# Pauli
hp = PauliSpace(:pauli)
σx = Pauli(hp, :σ, 1)
σy = Pauli(hp, :σ, 2)
σz = Pauli(hp, :σ, 3)

# Spin
hs = SpinSpace(:spin, 1 // 2)
Sx = Spin(hs, :S, 1)
Sy = Spin(hs, :S, 2)
Sz = Spin(hs, :S, 3)

# PhaseSpace
hps = PhaseSpace(:ps)
X = Position(hps, :X)
P = Momentum(hps, :P)

# Multi-space
hm = FockSpace(:c1) ⊗ FockSpace(:c2)
@qnumbers b::Destroy(hm, 1) c::Destroy(hm, 2)

# --- Benchmarks ---

SUITE = BenchmarkGroup()

# 1. Multiplication (exercises sort + vector construction)
SUITE["multiply"] = BenchmarkGroup()
SUITE["multiply"]["QSym*QSym"] = @benchmarkable $a * $(a')
SUITE["multiply"]["QMul*QMul small"] = @benchmarkable $(a' * a) * $(a' * a)
SUITE["multiply"]["QMul*QMul multi-space"] = @benchmarkable $(b' * c') * $(b * c)

# 2. Power (exercises ^)
SUITE["power"] = BenchmarkGroup()
SUITE["power"]["QSym^5"] = @benchmarkable $(a')^5
SUITE["power"]["QMul^3"] = @benchmarkable $(a' * a)^3

# 3. Addition (exercises QAdd construction)
SUITE["addition"] = BenchmarkGroup()
term1 = a' * a
term2 = 2 * a' * a
term3 = a * a'
SUITE["addition"]["build sum of 2"] = @benchmarkable $term1 + $term2
SUITE["addition"]["build sum of 4"] = @benchmarkable ($term1 + $term2) + ($term3 + $term1)

# 4. Simplify — Fock normal ordering
SUITE["simplify"] = BenchmarkGroup()
SUITE["simplify"]["a*a'"] = @benchmarkable simplify($(a * a'))
SUITE["simplify"]["(a*a')^2"] = @benchmarkable simplify($((a * a')^2))
SUITE["simplify"]["(a*a')^3"] = @benchmarkable simplify($((a * a')^3))
SUITE["simplify"]["(a*a')^4"] = @benchmarkable simplify($((a * a')^4))

# 5. Simplify — Transition composition
SUITE["simplify"]["Transition chain"] = @benchmarkable simplify($(σ12 * σ23 * σ13'))
SUITE["simplify"]["Transition sum"] = @benchmarkable simplify($(σ12 * σ12' + σ13 * σ13' + σ23 * σ23'))

# 6. Simplify — Pauli algebra
SUITE["simplify"]["Pauli σx*σy*σz"] = @benchmarkable simplify($(σx * σy * σz))
SUITE["simplify"]["Pauli sum"] = @benchmarkable simplify($(σx^2 + σy^2 + σz^2))

# 7. Simplify — Spin ordering
SUITE["simplify"]["Spin Sz*Sy*Sx"] = @benchmarkable simplify($(Sz * Sy * Sx))

# 8. Simplify — PhaseSpace
SUITE["simplify"]["P*X"] = @benchmarkable simplify($(P * X))
SUITE["simplify"]["(P*X)^2"] = @benchmarkable simplify($((P * X)^2))

# 9. Simplify — multi-space
SUITE["simplify"]["multi-space b*c*b'*c'"] = @benchmarkable simplify($(b * c * b' * c'))

# 10. Collect-like-terms stress (many duplicate terms)
many_terms = sum(a' * a for _ in 1:50) + sum(a * a' for _ in 1:50)
SUITE["simplify"]["collect 100 terms"] = @benchmarkable simplify($many_terms)

# 11. Normal order with ground state
SUITE["normal_order"] = BenchmarkGroup()
SUITE["normal_order"]["ground state rewrite"] = @benchmarkable normal_order(
    $(σ12 * σ12' + σ13 * σ13'), $hn
)

# 12. Commutator
SUITE["commutator"] = BenchmarkGroup()
SUITE["commutator"]["[a, a'] Fock"] = @benchmarkable commutator($a, $(a'))
SUITE["commutator"]["[a, a] zero"] = @benchmarkable commutator($a, $a)
SUITE["commutator"]["[a'a, a] QMul,QSym"] = @benchmarkable commutator($(a' * a), $a)
SUITE["commutator"]["[Sx, Sy] Spin"] = @benchmarkable commutator($Sx, $Sy)
SUITE["commutator"]["[X, P] PhaseSpace"] = @benchmarkable commutator($X, $P)
SUITE["commutator"]["[a+a', a] QAdd"] = @benchmarkable commutator($(a + a'), $a)
SUITE["commutator"]["nested depth=3"] = let term = a' * a
    function nested(t, n)
        n == 0 && return t
        return commutator(t, nested(t, n - 1))
    end
    @benchmarkable nested($term, 3)
end

# --- Run ---
BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)

println("\n" * "="^70)
println("RESULTS (median)")
println("="^70)
for (group_name, group) in sort(collect(results), by=first)
    println("\n[$group_name]")
    for (bench_name, trial) in sort(collect(group), by=first)
        m = median(trial)
        println("  $bench_name: $(BenchmarkTools.prettytime(time(m)))  ($(BenchmarkTools.prettymemory(memory(m))), $(allocs(m)) allocs)")
    end
end
