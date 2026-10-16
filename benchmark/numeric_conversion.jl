using BenchmarkTools: BenchmarkGroup, @benchmarkable
using SecondQuantizedAlgebra: FockSpace, NLevelSpace, Destroy, Transition, Index,
    IndexedOperator, Σ, ⊗, expim, to_numeric
using Symbolics: @variables
import QuantumOpticsBase as QOB

function benchmark_numeric_conversion!(SUITE::BenchmarkGroup)
    atom = NLevelSpace(:atom, 2)
    h = FockSpace(:cavity) ⊗ atom
    a = Destroy(h, :a, 1)
    σge = Transition(h, :σ, 1, 2, 2)
    σee = Transition(h, :σ, 2, 2, 2)
    @variables Δ::Real ω::Real g::Real η::Real t::Real
    H = Δ * a' * a + ω * σee + g * (a' * σge + a * σge')
    basis = QOB.FockBasis(20) ⊗ QOB.NLevelBasis(2)
    parameter = Dict(Δ => 0.5, ω => 1.0, g => 0.1)
    group = SUITE["Numeric conversion"]["Hamiltonians"]
    group["Static JC"] = @benchmarkable to_numeric($H, $basis; parameter = $parameter) seconds = 3 evals = 1

    Hd = H + η * (expim(ω * t) * a' + expim(-ω * t) * a)
    driven_parameter = merge(parameter, Dict(η => 0.2))
    time_parameter = Dict(t => identity)
    group["Time-dependent JC"] = @benchmarkable to_numeric(
        $Hd, $basis; parameter = $driven_parameter, time_parameter = $time_parameter,
    ) seconds = 3 evals = 1

    i = Index(h, :i, 4, atom)
    σge_i = IndexedOperator(σge, i)
    σee_i = IndexedOperator(σee, i)
    Htc = Δ * a' * a + Σ(ω * σee_i + g * (a' * σge_i + a * σge_i'), i)
    indexed_basis = QOB.FockBasis(5) ⊗ QOB.NLevelBasis(2)^4
    sites = Dict(1 => [1], 2 => collect(2:5))
    group["Indexed TC, 4 atoms"] = @benchmarkable to_numeric(
        $Htc, $indexed_basis, $sites; parameter = $parameter,
    ) seconds = 3 evals = 1
    return SUITE
end
