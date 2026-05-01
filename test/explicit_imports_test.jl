using SecondQuantizedAlgebra
using Test

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_qualified_accesses_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
end
