using SecondQuantizedAlgebra, Test

@testset "Concretely typed" begin
    using SecondQuantizedAlgebra
    using CheckConcreteStructs
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        JET.test_package(SecondQuantizedAlgebra; target_defined_modules=true)
        rep = report_package("SecondQuantizedAlgebra")
        @show rep
        @test length(JET.get_reports(rep)) <= 30
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) == nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) == nothing
end

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([SecondQuantizedAlgebra]; broken=false)
    Aqua.test_all(SecondQuantizedAlgebra; ambiguities=false)
end
