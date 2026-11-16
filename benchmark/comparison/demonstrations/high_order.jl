#=
    julia --project=benchmark benchmark/comparison/demonstrations/high_order.jl

Order-by-order Schrieffer–Wolff calculation for the dispersive Jaynes–Cummings
model, using the same graded BCH algorithm over SQA and QuantumAlgebra.
The generator is the known JC rotation expanded in ε = g/Δ, not a general
algorithm for finding Schrieffer–Wolff generators. All coefficients are exact.
=#
import SecondQuantizedAlgebra as SQA
import QuantumAlgebra as QA
import JSON
import CairoMakie as Makie
using BenchmarkTools: @benchmarkable, run, minimum, time, memory, allocs
using Test: @inferred
using SHA: sha256, bytes2hex
using Dates: today

function canonical_power(base, exponent::Int, normal::F) where {F}
    value = base^0
    for step in 1:exponent
        value = normal(value * base)
    end
    return value
end

function graded_commutator(left::Vector{T}, right::Vector{T}, order::Int, comm::F) where {T, F}
    result = fill(zero(first(right)), order + 1)
    for i in 1:order, j in 0:(order - i)
        (iszero(left[i + 1]) || iszero(right[j + 1])) && continue
        result[i + j + 1] += comm(left[i + 1], right[j + 1])
    end
    return result
end

function derive_sw(model::NamedTuple, order::Int)
    zero_expr = zero(model.excited)
    generator = fill(zero_expr, order + 1)
    number_power = model.excitation^0
    for k in 0:((order - 1) ÷ 2)
        coefficient = (-1)^(k + 1) * big(4)^k // (2k + 1)
        generator[2k + 2] = model.normal(coefficient * number_power * model.rotation)
        number_power = model.normal(number_power * model.excitation)
    end
    term = fill(zero_expr, order + 1)
    term[1] = model.excited
    term[2] = model.interaction
    result = copy(term)
    for depth in 1:order
        term = graded_commutator(generator, term, order, model.comm)
        term = [(big(1) // depth) * coefficient for coefficient in term]
        result = result .+ term
    end
    return model.normal.(result)
end

function models()
    h = SQA.tensor(SQA.FockSpace(:cavity), SQA.NLevelSpace(:atom, 2, 1))
    a = SQA.Destroy(h, :a, 1)
    lowering = SQA.Transition(h, :σ, 1, 2, 2)
    excited = lowering' * lowering
    ground = lowering * lowering'
    number = a' * a
    absorption = a * lowering'
    emission = a' * lowering
    sqa = (;
        name = "SQA", version = string(pkgversion(SQA)), excited, ground, number,
        excitation = number + excited, interaction = emission + absorption,
        rotation = emission - absorption, normal = identity, comm = SQA.commutator,
        is_zero = x -> iszero(SQA.expand_completeness(x)), terms = length,
    )
    QA.auto_normal_form(false)
    qa_number = QA.normal_form(QA.a'() * QA.a())
    qa_excited = QA.normal_form(QA.σp() * QA.σm())
    qa_ground = QA.normal_form(QA.σm() * QA.σp())
    qa_emission = QA.normal_form(QA.a'() * QA.σm())
    qa_absorption = QA.normal_form(QA.a() * QA.σp())
    qa = (;
        name = "QuantumAlgebra", version = string(pkgversion(QA)),
        excited = qa_excited, ground = qa_ground, number = qa_number,
        excitation = QA.normal_form(qa_number + qa_excited),
        interaction = qa_emission + qa_absorption, rotation = qa_emission - qa_absorption,
        normal = QA.normal_form, comm = (x, y) -> QA.normal_form(QA.comm(x, y)),
        is_zero = x -> iszero(QA.normal_form(x)), terms = x -> length(x.terms),
    )
    return (sqa, qa)
end

function verify_sw(model::NamedTuple, coefficients::Vector, order::Int)
    @assert model.is_zero(coefficients[1] - model.excited)
    for degree in 1:order
        expected = zero(model.excited)
        if iseven(degree)
            k = degree ÷ 2
            catalan = binomial(big(2k - 2), k - 1) ÷ k
            excited_power = canonical_power(model.number + 1, k, model.normal)
            ground_power = canonical_power(model.number, k, model.normal)
            expected = (-1)^(k + 1) * catalan * (excited_power * model.excited - ground_power * model.ground)
        end
        @assert model.is_zero(coefficients[degree + 1] - expected) "$(model.name), degree $degree"
    end
    # The independently known lower eigenvalue of each JC excitation block.
    epsilon = BigFloat("0.1")
    error = maximum(0:4) do n
        predicted = zero(epsilon)
        for k in 1:(order ÷ 2)
            catalan = binomial(big(2k - 2), k - 1) ÷ k
            predicted += (-1)^k * catalan * epsilon^(2k) * n^k
        end
        exact = (1 - sqrt(1 + 4epsilon^2 * n)) / 2
        abs(predicted - exact)
    end
    return Float64(error)
end

function measure_model(model::NamedTuple, orders::StepRange{Int, Int})
    @inferred derive_sw(model, 2)
    rows = []
    for order in orders
        result = derive_sw(model, order)
        error = verify_sw(model, result, order)
        trial = run(@benchmarkable(derive_sw($model, $order)); seconds = 2, samples = 10, evals = 1)
        estimate = minimum(trial)
        push!(
            rows, Dict(
                "package" => model.name, "order" => order,
                "time_ns" => time(estimate), "bytes" => memory(estimate), "allocs" => allocs(estimate),
                "samples_ns" => trial.times, "native_terms" => sum([model.terms(coefficient) for coefficient in result]),
                "max_energy_error" => error, "verified" => true
            )
        )
        println(model.name, " order ", order, ": ", time(estimate) / 1.0e6, " ms; error ", error)
        flush(stdout)
    end
    return rows
end

function verify_orders(model::NamedTuple, orders::StepRange{Int, Int})
    @inferred derive_sw(model, 2)
    for order in orders
        result = derive_sw(model, order)
        println(model.name, " verified order ", order, "; energy error ", verify_sw(model, result, order))
    end
    return
end

function main()
    maximum_order = parse(Int, get(ENV, "SQA_SW_MAX_ORDER", "20"))
    maximum_order >= 2 && iseven(maximum_order) || error("Choose an even order >= 2")
    pair = models()
    if get(ENV, "SQA_SW_VALIDATE_ONLY", "false") == "true"
        for model in pair
            verify_orders(model, 2:2:maximum_order)
        end
        return nothing
    end
    rows = reduce(vcat, (measure_model(model, 2:2:maximum_order) for model in pair))
    output = get(ENV, "SQA_DEMO_OUTPUT", joinpath(@__DIR__, "..", "results"))
    mkpath(output)
    meta = Dict(
        "date" => string(today()), "julia" => string(VERSION),
        "cpu" => Sys.cpu_info()[1].model, "threads" => Threads.nthreads(),
        "versions" => Dict(model.name => model.version for model in pair),
        "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "algorithm" => "graded BCH with analytic JC rotation generator",
        "epsilon" => 0.1, "verified_photon_numbers" => collect(0:4),
        "coefficient_domain" => "exact rational; formal epsilon degree",
        "energy_error_branch" => "ground-like",
        "statistic" => "minimum warm runtime; compilation and verification excluded"
    )
    open(joinpath(output, "high_order.json.tmp"), "w") do io
        JSON.print(io, Dict("meta" => meta, "results" => rows), 2)
    end
    mv(joinpath(output, "high_order.json.tmp"), joinpath(output, "high_order.json"); force = true)
    fig = Makie.Figure(; size = (1050, 430))
    runtime = Makie.Axis(fig[1, 1]; xlabel = "perturbative order", ylabel = "runtime (s)", yscale = log10)
    accuracy = Makie.Axis(
        fig[1, 2]; xlabel = "perturbative order", ylabel = "ground-like energy error / Δ", yscale = log10,
        title = "g/Δ = 0.1; photon numbers 0:4"
    )
    for model in pair
        selected = filter(row -> row["package"] == model.name, rows)
        Makie.scatterlines!(runtime, [row["order"] for row in selected], [row["time_ns"] / 1.0e9 for row in selected]; label = model.name, color = model.name == "SQA" ? "#2a78d6" : "#008300", marker = model.name == "SQA" ? :circle : :rect)
    end
    selected = filter(row -> row["package"] == "SQA", rows)
    Makie.scatterlines!(accuracy, [row["order"] for row in selected], [row["max_energy_error"] for row in selected])
    Makie.axislegend(runtime; position = :lt)
    for extension in ("svg", "png")
        Makie.save(joinpath(output, "high_order.$extension"), fig)
    end
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
