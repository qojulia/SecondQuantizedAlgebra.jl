using SecondQuantizedAlgebra
using ScopedValues: ScopedValues
using Test

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
    @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
    @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
    # ScopedValues.jl re-exports `Base.ScopedValues` on Julia ≥1.11 via a const
    # alias (`const get = Base.ScopedValues.get`), so ExplicitImports resolves the
    # owner of `ScopedValues.get` to `Base.ScopedValues` rather than the package
    # we depend on. The skip records that this re-export is the legitimate path.
    # On Julia 1.10, `Base.ScopedValues` does not exist (the standalone package
    # is the only source), so the skip is unnecessary and would itself fail.
    skip_tuple = if VERSION ≥ v"1.11"
        (Base => Core, ScopedValues => Base.ScopedValues)
    else
        (Base => Core,)
    end
    @test check_all_qualified_accesses_via_owners(
        SecondQuantizedAlgebra; skip = skip_tuple,
    ) === nothing
    @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
end
