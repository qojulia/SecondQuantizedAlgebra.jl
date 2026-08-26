using Test
using BenchmarkTools
using SecondQuantizedAlgebra
include(joinpath(@__DIR__, "product_api_spike.jl"))
using .ProductAPIPrototype

function onlyterm(expr)
    q = expr isa Op ? expr + 0 : expr
    terms = [term for (term, _) in q if !isempty(term.ops)]
    return only(terms)
end

function fixtures()
    # Duplicate factor types and renamed operators must remain position-addressable.
    hf = FockSpace(:left) ⊗ FockSpace(:right)
    a = Destroy(hf, :renamed_left, 1)
    b = Destroy(hf, :renamed_right, 2)
    Gg = a' * b + b' * a + a * b + a' * b'
    gaussian = FlatCoordinates(hf, a, b, a', b')

    # Finite cross-site closed adjoint: {σx₁, σy₁σz₂} under σz₁σz₂.
    hp = PauliSpace(:p_left) ⊗ PauliSpace(:p_right)
    sx1, sy1, sz1 = (Pauli(hp, :leftσ, k, 1) for k in 1:3)
    sz2 = Pauli(hp, :rightτ, 3, 2)
    Gp = sz1 * sz2
    cross = FlatCoordinates(hp, sx1, sy1 * sz2)

    # Ordinary transitions: non-default ground, concrete indexed site.
    hn = NLevelSpace(:atoms, 3, 2)
    i = Index(hn, :i, 4, hn)
    j = SecondQuantizedAlgebra.rename(i, :j)
    σ13 = Transition(hn, :renamedσ, 1, 3)
    σi = IndexedOperator(σ13, i(2))
    σj = IndexedOperator(σ13, j(2))

    # Collective transition remains unindexed and validates from explicit space.
    hc = CollectiveNLevelSpace(:ensemble, 3)
    S12 = CollectiveTransition(hc, :renamedS, 1, 2)

    return (; hf, a, b, Gg, gaussian, hp, sx1, sy1, sz1, sz2, Gp, cross,
        hn, i, j, σ13, σi, σj, hc, S12)
end

const FX = fixtures()

@testset "flat ProductSpace and API pressure" begin
    Ag = @inferred project_action(FX.Gg, FX.gaussian)
    @test size(Ag) == (4, 4)
    @test count(!iszero, Ag) == 8
    expected_g = similar(Ag)
    fill!(expected_g, 0)
    expected_g[2, 1] = -im; expected_g[4, 1] = -im
    expected_g[1, 2] = -im; expected_g[3, 2] = -im
    expected_g[2, 3] = im;  expected_g[4, 3] = im
    expected_g[1, 4] = im;  expected_g[3, 4] = im
    @test Ag == expected_g
    @test FX.gaussian.support == [[1], [2], [1], [2]]

    Ap = @inferred project_action(FX.Gp, FX.cross)
    @test size(Ap) == (2, 2)
    @test count(!iszero, Ap) == 2
    @test FX.cross.support == [[1], [1, 2]]
    @test iszero(Ap[1, 1]) && iszero(Ap[2, 2])
    @test Ap[1, 2] == 2 && Ap[2, 1] == -2
    @test Ap[1, 2] * Ap[2, 1] == -4

    @test @inferred(validate_hilbert(FX.hn, FX.σ13, FX.σi)) === FX.hn
    @test @inferred(validate_hilbert(FX.hc, FX.S12)) === FX.hc
    @test index_slot(operator_index(FX.σi)) == 2
    @test_throws ArgumentError validate_hilbert(FX.hn, IndexedOperator(FX.σ13, FX.i))

    wrong_transition_space = NLevelSpace(:wrong_dimension, 4, 2)
    @test_throws ArgumentError validate_hilbert(wrong_transition_space, FX.σ13)

    ProductAPIPrototype.reset_analysis_count!()
    result = @inferred analyze_frame(FX.Gg; hilbert = FX.hf)
    @test result isa AnalysisResult{typeof(FX.hf)}
    @test ProductAPIPrototype.analysis_count() == 1

    wrong = FockSpace(:only_one)
    ProductAPIPrototype.reset_analysis_count!()
    @test_throws ArgumentError analyze_frame(FX.Gg; hilbert = wrong)
    @test ProductAPIPrototype.analysis_count() == 0

    ProductAPIPrototype.reset_analysis_count!()
    @test_throws ArgumentError analyze_frame(FX.Gg)
    @test ProductAPIPrototype.analysis_count() == 0

    swapped = PauliSpace(:wrong1) ⊗ FockSpace(:wrong2)
    ProductAPIPrototype.reset_analysis_count!()
    @test_throws ArgumentError analyze_frame(FX.Gg; hilbert = swapped)
    @test ProductAPIPrototype.analysis_count() == 0

    # Raw QTerm identity intentionally sees labels; a future explicit binder descriptor can
    # canonicalize them without weakening Dict equality.
    ti, tj = onlyterm(FX.σi), onlyterm(FX.σj)
    @test ti != tj
    @test alpha_signature(ti) == alpha_signature(tj)

    # Inference cannot recover untouched factors, Hilbert names, or transition level names.
    incomplete = FlatCoordinates(FX.hf, FX.a)
    @test incomplete.support == [[1]]
    @test length(FX.hf) == 2
end

function metrics()
    gaussian_trial = @benchmark project_action($FX.Gg, $FX.gaussian) samples=20 evals=1
    cross_trial = @benchmark project_action($FX.Gp, $FX.cross) samples=20 evals=1
    adapter_trial = @benchmark analyze_frame($FX.Gg; hilbert=$FX.hf) samples=20 evals=1
    support_vector_bytes = Base.summarysize(FX.gaussian.support)
    support_mask_bytes = Base.summarysize(UInt64[sum(UInt64(1) << (Int(si) - 1) for si in s) for s in FX.gaussian.support])
    basis_bytes = Base.summarysize(FX.gaussian)
    action_bytes = Base.summarysize(project_action(FX.Gg, FX.gaussian))
    return (; gaussian_time = median(gaussian_trial).time,
        gaussian_memory = median(gaussian_trial).memory,
        gaussian_allocs = median(gaussian_trial).allocs,
        cross_time = median(cross_trial).time,
        cross_memory = median(cross_trial).memory,
        cross_allocs = median(cross_trial).allocs,
        adapter_time = median(adapter_trial).time,
        adapter_memory = median(adapter_trial).memory,
        adapter_allocs = median(adapter_trial).allocs,
        support_vector_bytes, support_mask_bytes, basis_bytes, action_bytes)
end

println("METRICS=", metrics())
println("GAUSSIAN_ACTION=", project_action(FX.Gg, FX.gaussian))
println("CROSS_ACTION=", project_action(FX.Gp, FX.cross))
println("METHODS=", methods(analyze_frame))
println("PRECOMPILE_POSITIONAL=", precompile(ProductAPIPrototype._analyze_frame,
    (typeof(FX.Gg), typeof(FX.hf))))

try
    using JET
    report = JET.report_call(project_action, (typeof(FX.Gg), typeof(FX.gaussian));
        target_modules = (ProductAPIPrototype,))
    println("JET_PROJECT=", report)
    report2 = JET.report_call(ProductAPIPrototype._adapter_call,
        (typeof(FX.Gg), typeof(FX.hf)); target_modules = (ProductAPIPrototype,))
    println("JET_ADAPTER=", report2)
catch err
    println("JET_SKIPPED=", sprint(showerror, err))
end

try
    using Aqua
    Aqua.test_ambiguities([ProductAPIPrototype]; recursive = false)
    Aqua.test_piracies(ProductAPIPrototype)
    println("AQUA=pass")
catch err
    println("AQUA_SKIPPED_OR_FAILED=", sprint(showerror, err))
end


println("AMBIGUITIES=", Test.detect_ambiguities(ProductAPIPrototype; recursive = false))
println("DOC_DISCOVERY=", Base.Docs.doc(analyze_frame) !== nothing)

