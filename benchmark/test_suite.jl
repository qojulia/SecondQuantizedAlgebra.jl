using Test: @test, @testset, @test_throws
using BenchmarkTools: BenchmarkGroup, @benchmarkable, leaves
using JSON: JSON

include("suite.jl")

@testset "Benchmark hierarchy" begin
    suite = benchmark_suite()
    paths = Set(join(path, '/') for (path, _) in leaves(suite))
    renames = JSON.parsefile(joinpath(@__DIR__, "hierarchy.json"))
    targets = [join(path, '/') for path in values(renames)]
    @test length(unique(targets)) == length(targets)
    @test all(path -> path in paths, targets)
    @test all(path -> length(split(path, '/')) == 3, paths)

    shallow = BenchmarkGroup()
    shallow["topic"]["case"] = @benchmarkable sin(1.0)
    @test_throws ErrorException validate_suite(shallow)
    deep = BenchmarkGroup()
    deep["topic"]["case"]["legend"]["extra"] = @benchmarkable sin(1.0)
    @test_throws ErrorException validate_suite(deep)
    slash = BenchmarkGroup()
    slash["topic"]["case"]["invalid/legend"] = @benchmarkable sin(1.0)
    @test_throws ErrorException validate_suite(slash)
    @test_throws ErrorException validate_suite(BenchmarkGroup())
end
