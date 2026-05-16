using SecondQuantizedAlgebra
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports
using QuantumOpticsBase: expect

@testset "Quality gates" begin
    @testset "Aqua" begin
        # TODO: remove `treat_as_own = [expect]` once the piracies in src/numeric.jl
        # are resolved — `expect(::Num, ...)` and `expect(::BasicSymbolic, ...)` are
        # methods on a QuantumOpticsBase function with externally-owned argument types.
        Aqua.test_all(SecondQuantizedAlgebra; piracies = (treat_as_own = [expect],))
    end

    @testset "ExplicitImports" begin
        @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
        @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
        @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
        @test check_all_qualified_accesses_via_owners(
            SecondQuantizedAlgebra; skip = (Base => Core,),
        ) === nothing
        @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
    end

    @testset "CheckConcreteStructs" begin
        for name in names(SecondQuantizedAlgebra; all = true)
            isdefined(SecondQuantizedAlgebra, name) || continue
            T = getfield(SecondQuantizedAlgebra, name)
            T isa Type || continue
            isabstracttype(T) && continue
            T isa UnionAll && continue
            isstructtype(T) || continue
            parentmodule(T) === SecondQuantizedAlgebra || continue
            @testset "$name" begin
                @test all_concrete(T; verbose = false)
            end
        end
    end

    # JET.report_package is infeasible: methods taking Symbolics types (BasicSymbolic,
    # Complex{Num}) send JET deep into Symbolics.jl internals, causing unbounded analysis
    # time. Even max_methods=2 doesn't help. Use @test_call on concrete entry points
    # instead — JET then sees concrete types (e.g. Destroy, not abstract QSym) and
    # avoids the combinatorial subtype × Symbolics blowup.
    # @static if isempty(VERSION.prerelease)
    #     @testset "JET" begin
    #         using JET
    #         JET.test_package(SecondQuantizedAlgebra;
    #             target_modules = (SecondQuantizedAlgebra,))
    #     end
    # end
end
