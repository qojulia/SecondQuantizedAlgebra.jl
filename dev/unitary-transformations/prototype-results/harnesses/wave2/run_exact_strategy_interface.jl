module ExactStrategyRunner

using Test
using LinearAlgebra
using Statistics
using JET

include(joinpath(@__DIR__, "exact_strategy_interface.jl"))
using .ExactStrategyPrototype

function fixtures()
    static = NoScalarLift()
    timed = ScalarLift(7 // 11)

    diagonal = exact_certificate(
        DiagonalPhaseBlock(
            Complex{Rational{Int}}[3 // 5 + 4 // 5 * im, 5 // 13 + 12 // 13 * im],
            Complex{Rational{Int}}[3 // 5 - 4 // 5 * im, 5 // 13 - 12 // 13 * im],
            timed,
        ),
    )
    rotation = exact_certificate(RotationBlock(3 // 5, 4 // 5, static))
    squeeze = exact_certificate(SqueezeBlock(5 // 4, 3 // 4, timed))

    complex_structure = Rational{Int}[0 -1; 1 0]
    involution_minus = exact_certificate(
        NormalizedInvolutionBlock(Val(-1), complex_structure, 3 // 5, 4 // 5),
    )
    reflection = Rational{Int}[1 0; 0 -1]
    involution_plus = exact_certificate(
        NormalizedInvolutionBlock(Val(1), reflection, 5 // 4, 3 // 4),
    )

    nilpotent_action = Rational{Int}[0 1 0; 0 0 1; 0 0 0]
    nilpotent = exact_certificate(NilpotentBlock(Val(3), nilpotent_action, 2 // 3))

    user_forward = Rational{Int}[3//5 -4//5; 4//5 3//5]
    user_inverse = transpose(user_forward) |> Matrix
    user = exact_certificate(UserMapBlock(user_forward, user_inverse))

    unsupported = exact_certificate(UnsupportedBlock(Rational{Int}[1 2; 3 4]))
    branch_unsafe = exact_certificate(
        ScaledInvolutionBlock(Rational{Int}[0 2; 2 0], 4 // 1),
    )
    return (; diagonal, rotation, squeeze, involution_minus, involution_plus,
        nilpotent, user, unsupported, branch_unsafe)
end

function run_correctness()
    f = fixtures()
    @testset "method-based exact strategy" begin
        for certificate in (
                f.diagonal, f.rotation, f.squeeze, f.involution_minus,
                f.involution_plus, f.nilpotent, f.user,
            )
            n = size(certificate.forward, 1)
            @test certificate.forward * certificate.inverse ==
                Matrix{eltype(certificate.forward)}(I, n, n)
            @test certificate.inverse * certificate.forward ==
                Matrix{eltype(certificate.forward)}(I, n, n)
            @test certificate.verification.checks == 0x03
            @test isconcretetype(typeof(certificate))
            @test all(isconcretetype, fieldtypes(typeof(certificate)))
        end
        @test f.diagonal.scalar_lift == ScalarLift(7 // 11)
        @test f.squeeze.scalar_lift == ScalarLift(7 // 11)
        @test f.rotation.scalar_lift isa NoScalarLift
        @test f.unsupported === nothing
        @test f.branch_unsafe === nothing

        derivative = Rational{Int}[0 -1; 1 0]
        @test body_velocity(f.rotation, derivative) == f.rotation.inverse * derivative

        # Wave 1 canonical order is (a1,a1',a2,a2'). A split view is temporary.
        site_map = Rational{Int}[
            1 0 2 0;
            0 1 0 2;
            3 0 4 0;
            0 3 0 4
        ]
        split_permutation = [1, 3, 2, 4]
        split = split_nambu_map(site_map, split_permutation)
        @test site_interleaved_map(split, split_permutation) == site_map

        @test_throws ArgumentError exact_certificate(RotationBlock(1 // 2, 1 // 2))
        @test_throws ArgumentError exact_certificate(SqueezeBlock(1 // 1, 1 // 1))
        @test_throws ArgumentError exact_certificate(
            NormalizedInvolutionBlock(Val(-1), Rational{Int}[1 0; 0 1], 1, 0),
        )
        @test_throws ArgumentError exact_certificate(
            NilpotentBlock(Val(2), Rational{Int}[0 1; 1 0], 1),
        )
        @test_throws ArgumentError exact_certificate(
            UserMapBlock(Rational{Int}[1 1; 0 1], Rational{Int}[1 0; 0 1]),
        )
    end
    return f
end

function inference_checks()
    phase = Complex{Rational{Int}}[3 // 5 + 4 // 5 * im]
    inverse_phase = conj.(phase)
    diagonal = DiagonalPhaseBlock(phase, inverse_phase)
    rotation = RotationBlock(3 // 5, 4 // 5)
    squeeze = SqueezeBlock(5 // 4, 3 // 4)
    action = Rational{Int}[0 -1; 1 0]
    involution = NormalizedInvolutionBlock(Val(-1), action, 3 // 5, 4 // 5)
    nilpotent = NilpotentBlock(Val(2), Rational{Int}[0 1; 0 0], 2 // 3)
    user_forward = Rational{Int}[3//5 -4//5; 4//5 3//5]
    user = UserMapBlock(user_forward, transpose(user_forward) |> Matrix)
    unsupported = UnsupportedBlock(action)

    @test (@inferred exact_certificate(diagonal)) isa ExactCertificate
    @test (@inferred exact_certificate(rotation)) isa ExactCertificate
    @test (@inferred exact_certificate(squeeze)) isa ExactCertificate
    @test (@inferred exact_certificate(involution)) isa ExactCertificate
    @test (@inferred exact_certificate(nilpotent)) isa ExactCertificate
    @test (@inferred exact_certificate(user)) isa ExactCertificate
    @test (@inferred exact_certificate(unsupported)) === nothing
    return (; diagonal, rotation, squeeze, involution, nilpotent, user, unsupported)
end

function warm_measure(f; samples = 31)
    f()
    times = Vector{Float64}(undef, samples)
    bytes = Vector{Int}(undef, samples)
    for i in eachindex(times)
        bytes[i] = @allocated begin
            start = time_ns()
            f()
            times[i] = time_ns() - start
        end
    end
    return (median_ns = median(times), median_bytes = median(bytes))
end

function measurements(blocks)
    return (
        diagonal_2 = warm_measure(() -> exact_certificate(blocks.diagonal)),
        rotation_2 = warm_measure(() -> exact_certificate(blocks.rotation)),
        squeeze_2 = warm_measure(() -> exact_certificate(blocks.squeeze)),
        involution_2 = warm_measure(() -> exact_certificate(blocks.involution)),
        nilpotent_2 = warm_measure(() -> exact_certificate(blocks.nilpotent)),
        user_2 = warm_measure(() -> exact_certificate(blocks.user)),
        unsupported = warm_measure(() -> exact_certificate(blocks.unsupported)),
    )
end

function jet_reports(blocks)
    reports = (
        rotation = JET.report_call(exact_certificate, (typeof(blocks.rotation),)),
        nilpotent = JET.report_call(exact_certificate, (typeof(blocks.nilpotent),)),
        user = JET.report_call(exact_certificate, (typeof(blocks.user),)),
        unsupported = JET.report_call(exact_certificate, (typeof(blocks.unsupported),)),
    )
    optimization = (
        rotation = JET.report_opt(exact_certificate, (typeof(blocks.rotation),)),
        nilpotent = JET.report_opt(exact_certificate, (typeof(blocks.nilpotent),)),
        unsupported = JET.report_opt(exact_certificate, (typeof(blocks.unsupported),)),
    )
    return reports, optimization
end

function cnum_checks()
    C = ExactStrategyPrototype.SQA._to_cnum
    block = RotationBlock(C(3 // 5), C(4 // 5))
    certificate = @inferred exact_certificate(block)
    @test certificate isa ExactCertificate
    @test isconcretetype(typeof(certificate))
    @test all(isconcretetype, fieldtypes(typeof(certificate)))
    correctness = JET.@report_call target_modules = (ExactStrategyPrototype,) exact_certificate(block)
    optimization = JET.@report_opt target_modules = (ExactStrategyPrototype,) exact_certificate(block)
    return (
        type = typeof(certificate),
        measurement = warm_measure(() -> exact_certificate(block)),
        correctness_reports = length(JET.get_reports(correctness)),
        optimization_reports = length(JET.get_reports(optimization)),
    )
end

f = run_correctness()
blocks = inference_checks()
perf = measurements(blocks)
ambiguities = Test.detect_ambiguities(ExactStrategyPrototype; recursive = true)
reports, optimization = jet_reports(blocks)
cnum = cnum_checks()

println("fixture_types=")
foreach(x -> println("  ", typeof(x)), values(f))
println("measurements=", perf)
println("ambiguities=", ambiguities)
println("method_count=", length(methods(exact_certificate)))
println("certificate_summary_bytes=", map(Base.summarysize, (
    f.diagonal, f.rotation, f.squeeze, f.involution_minus, f.nilpotent, f.user,
)))
println("cnum=", cnum)
println("jet_reports=")
foreach(pair -> println("  ", first(pair), " => ", last(pair)), pairs(reports))
println("jet_optimization_reports=")
foreach(pair -> println("  ", first(pair), " => ", last(pair)), pairs(optimization))

end

