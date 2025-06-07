using SecondQuantizedAlgebra, Test

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([SecondQuantizedAlgebra]; broken=true)
    Aqua.test_piracies(SecondQuantizedAlgebra; broken=false)
    Aqua.test_all(SecondQuantizedAlgebra; ambiguities=false, piracies=false)
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(SecondQuantizedAlgebra) == nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) == nothing
    # @test check_all_explicit_imports_are_public(SecondQuantizedAlgebra) == nothing
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) == nothing
    @test check_all_qualified_accesses_via_owners(SecondQuantizedAlgebra) == nothing
    # @test check_all_qualified_accesses_are_public(SecondQuantizedAlgebra) == nothing
    @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) == nothing
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        # JET.test_package(SecondQuantizedAlgebra; target_defined_modules=true)
        rep = report_package("SecondQuantizedAlgebra")
        @show rep
        @test length(JET.get_reports(rep)) <= 259
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "Concretely typed" begin
    using SecondQuantizedAlgebra
    using CheckConcreteStructs
end
