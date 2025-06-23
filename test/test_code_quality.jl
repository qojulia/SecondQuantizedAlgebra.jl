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
        @test length(JET.get_reports(rep)) <= 262
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "Concretely typed" begin
    import SecondQuantizedAlgebra as SQA
    using CheckConcreteStructs

    all_concrete(SQA.QMul)
    all_concrete(SQA.QAdd)

    all_concrete(SQA.AvgSym)

    all_concrete(SQA.Parameter)
    all_concrete(SQA.RealParameter)

    all_concrete(SQA.FockSpace)
    all_concrete(SQA.Destroy)
    all_concrete(SQA.Create)

    all_concrete(SQA.PauliSpace)
    all_concrete(SQA.Pauli)
    all_concrete(SQA.SpinSpace)
    all_concrete(SQA.Spin)

    all_concrete(SQA.NLevelSpace)
    all_concrete(SQA.Transition)
    all_concrete(SQA.CallableTransition)

    all_concrete(SQA.ProductSpace)
    all_concrete(SQA.ClusterSpace)
    all_concrete(SQA.ClusterAon)

    all_concrete(SQA.IndexedAverageSum)
    all_concrete(SQA.IndexedAverageDoubleSum)
    all_concrete(SQA.SingleNumberedVariable)
    all_concrete(SQA.DoubleNumberedVariable)
    all_concrete(SQA.SpecialIndexedAverage)

    all_concrete(SQA.DoubleSum)

    all_concrete(SQA.NumberedOperator)

    all_concrete(SQA.Index)
    all_concrete(SQA.IndexedVariable)
    all_concrete(SQA.DoubleIndexedVariable)
    all_concrete(SQA.IndexedOperator)
    all_concrete(SQA.SingleSum)
    all_concrete(SQA.SpecialIndexedTerm)
end
