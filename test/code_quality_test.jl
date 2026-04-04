using SecondQuantizedAlgebra
using Test

@testset "best practices" begin
    using Aqua
    Aqua.test_all(SecondQuantizedAlgebra; ambiguities=false, piracies=false)
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_qualified_accesses_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        rep = report_package(SecondQuantizedAlgebra)
        @show rep
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "Concretely typed" begin
    using CheckConcreteStructs
    all_concrete(SecondQuantizedAlgebra.FockSpace)
    all_concrete(SecondQuantizedAlgebra.ProductSpace)
    all_concrete(SecondQuantizedAlgebra.Destroy)
    all_concrete(SecondQuantizedAlgebra.Create)
    all_concrete(SecondQuantizedAlgebra.QMul)
    all_concrete(SecondQuantizedAlgebra.QAdd)
end
