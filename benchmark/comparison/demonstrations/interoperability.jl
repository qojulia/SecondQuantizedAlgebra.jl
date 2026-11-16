#=
    julia --project=benchmark/comparison/demonstrations benchmark/comparison/demonstrations/interoperability.jl

One SQA driven, damped mode feeds QuantumCumulants → ModelingToolkitBase →
OrdinaryDiffEq and a direct QuantumOptics master-equation simulation.
The retained moments close exactly for this quadratic Hamiltonian and linear
jump. Finite Fock truncation is checked at two dimensions. This is an
interoperability check, not a comparison of solver performance.
=#
import SecondQuantizedAlgebra as SQA
import QuantumCumulants as QC
import QuantumOptics as QO
import ModelingToolkitBase as MTK
using OrdinaryDiffEqLowOrderRK: RK4, solve
using ModelingToolkitBase: ODEProblem
using SecondQuantizedAlgebra: @variables
using SHA: sha256, bytes2hex
using Dates: today
import JSON

function model()
    h = SQA.FockSpace(:cavity)
    a = SQA.Destroy(h, :a)
    @variables Δ::Real η::Real κ::Real
    H = Δ * a' * a + η * (a + a')
    jumps = [a]
    rates = [κ]
    parameters = Dict(Δ => 0.7, η => 0.3, κ => 0.4)
    generator = SQA.hamiltonian_action(H) + κ * SQA.dissipator(a)
    return (; h, a, number = a' * a, H, jumps, rates, parameters, generator, Δ, η, κ)
end

function cumulant_solution(m, times)
    equations = QC.meanfield([m.a, m.number], m.H, m.jumps; rates = m.rates, order = 2)
    system = MTK.mtkcompile(MTK.System(equations; name = :driven_mode))
    initial = QC.initial_values(equations, ComplexF64[0, 0])
    problem = ODEProblem(system, merge(initial, m.parameters), (first(times), last(times)))
    solution = solve(problem, RK4(); abstol = 1.0e-12, reltol = 1.0e-10, saveat = times)
    amplitude = QC.get_solution(solution, m.a, equations).(times)
    occupation = real.(QC.get_solution(solution, m.number, equations).(times))
    return (; amplitude, occupation)
end

function quantum_solution(m, times, dimension::Int)
    backend = SQA.QuantumOpticsBackend()
    H, collapse = SQA.to_numeric_master(m.generator, m.h, dimension; backend, parameter = m.parameters)
    a = SQA.to_numeric(m.a, m.h, dimension; backend, parameter = m.parameters)
    number = SQA.to_numeric(m.number, m.h, dimension; backend, parameter = m.parameters)
    vacuum = QO.fockstate(QO.basis(H), 0)
    returned_times, states = QO.timeevolution.master(times, vacuum, H, collapse; abstol = 1.0e-12, reltol = 1.0e-10)
    @assert returned_times ≈ times
    return (;
        amplitude = [QO.expect(a, state) for state in states],
        occupation = [real(QO.expect(number, state)) for state in states],
    )
end

function main()
    m = model()
    times = collect(range(0, 10; length = 101))
    moments = cumulant_solution(m, times)
    direct = quantum_solution(m, times, 12)
    converged = quantum_solution(m, times, 20)
    decay = m.parameters[m.κ] / 2 + im * m.parameters[m.Δ]
    exact = [-im * m.parameters[m.η] / decay * (1 - exp(-decay * t)) for t in times]
    checks = Dict(
        "cumulant_amplitude_error" => maximum(abs.(moments.amplitude .- exact)),
        "cumulant_occupation_error" => maximum(abs.(moments.occupation .- abs2.(exact))),
        "quantum_amplitude_error" => maximum(abs.(converged.amplitude .- exact)),
        "quantum_occupation_error" => maximum(abs.(converged.occupation .- abs2.(exact))),
        "cross_method_amplitude_error" => maximum(abs.(moments.amplitude .- converged.amplitude)),
        "cross_method_occupation_error" => maximum(abs.(moments.occupation .- converged.occupation)),
        "cutoff_amplitude_error" => maximum(abs.(direct.amplitude .- converged.amplitude)),
        "cutoff_occupation_error" => maximum(abs.(direct.occupation .- converged.occupation))
    )
    @assert all(error < 1.0e-8 for error in values(checks)) checks
    meta = Dict(
        "date" => string(today()), "julia" => string(VERSION),
        "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "environment_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
        "versions" => Dict(
            "SQA" => string(pkgversion(SQA)), "QuantumCumulants" => string(pkgversion(QC)),
            "QuantumOptics" => string(pkgversion(QO)), "ModelingToolkitBase" => string(pkgversion(MTK))
        ),
        "quantumcumulants_revision" => "bfe9b198a375f81be061f9af6d852dfcbcaa18fe",
        "quantumcumulants_status" => "pinned development revision; registered releases do not resolve with SQA 0.14",
        "model" => "H = Δ a†a + η(a+a†); collapse sqrt(κ)a; vacuum initial state",
        "parameters" => Dict("Δ" => 0.7, "η" => 0.3, "κ" => 0.4),
        "fock_dimensions" => [12, 20], "moment_order" => 2, "closure" => "exact for retained moments"
    )
    results = Dict(
        "meta" => meta, "checks" => checks, "verified" => true,
        "time" => times, "exact_amplitude_real" => real.(exact), "exact_amplitude_imag" => imag.(exact),
        "cumulant_amplitude_real" => real.(moments.amplitude), "cumulant_amplitude_imag" => imag.(moments.amplitude),
        "quantum_amplitude_real" => real.(converged.amplitude), "quantum_amplitude_imag" => imag.(converged.amplitude),
        "cumulant_occupation" => moments.occupation, "quantum_occupation" => converged.occupation
    )
    output = get(ENV, "SQA_DEMO_OUTPUT", joinpath(@__DIR__, "..", "results"))
    mkpath(output)
    open(joinpath(output, "interoperability.json"), "w") do io
        JSON.print(io, results, 2)
    end
    return println("Interoperability checks passed: ", checks)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
