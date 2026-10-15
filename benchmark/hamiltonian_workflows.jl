using BenchmarkTools: BenchmarkGroup, @benchmarkable
using SecondQuantizedAlgebra: FockSpace, PhaseSpace, ProductSpace, Destroy,
    Position, Momentum, Displace, Squeeze, transform, substitute, taylor, ⊗
using Symbolics: @variables

function benchmark_hamiltonian_workflows!(SUITE::BenchmarkGroup)
    @variables ω J U
    for M in (4, 8, 16)
        h = ProductSpace(ntuple(k -> FockSpace(Symbol(:bh, k)), M))
        a = [Destroy(h, Symbol(:a, k), k) for k in 1:M]
        H = sum(ω * ak' * ak + (U / 2) * ak'^2 * ak^2 for ak in a) +
            J * sum(a[k]' * a[k + 1] + a[k + 1]' * a[k] for k in 1:(M - 1))
        SUITE["Multiplication"]["Bose-Hubbard H²"]["M=$M"] =
            @benchmarkable $H * $H seconds = 3 evals = 1
    end

    a = Destroy(FockSpace(:nonlinear), :a)
    @variables Δ::Real K::Real η::Real r::Real α::Number
    displacement = Displace(a, α)
    squeezing = Squeeze(a, r)
    group = SUITE["Unitary and exact phase workflows"]["Nonlinear transformations"]
    for (name, degree) in (("Kerr", 2), ("Sextic", 3))
        H = Δ * a' * a + K * a'^degree * a^degree + η * (a + a')
        group["$name displacement"] = @benchmarkable transform($H, $displacement) seconds = 3 evals = 1
        group["$name squeezing"] = @benchmarkable transform($H, $squeezing) seconds = 3 evals = 1
    end

    h = PhaseSpace(:transmon) ⊗ FockSpace(:resonator)
    φ = Position(h, :φ, 1)
    n = Momentum(h, :n, 1)
    a = Destroy(h, :a, 2)
    hb = FockSpace(:transmon_boson) ⊗ FockSpace(:resonator)
    b = Destroy(hb, :b, 1)
    @variables E_C::Real E_J::Real ωr::Real g::Real φ_zpf::Real n_zpf::Real
    H = 4 * E_C * n^2 - E_J * cos(φ) + ωr * a' * a + g * n * (a' + a)
    rules = Dict(φ => φ_zpf * (b + b'), n => im * n_zpf * (b' - b))
    for order in (4, 6, 8)
        SUITE["Formal operator functions"]["Transmon lowering"]["order=$order"] =
            @benchmarkable taylor(substitute($H, $rules), 0:$order) seconds = 3 evals = 1
    end
    return SUITE
end
