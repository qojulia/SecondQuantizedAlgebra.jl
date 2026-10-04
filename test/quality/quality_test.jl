using SecondQuantizedAlgebra
using Test
using Aqua
using CheckConcreteStructs: all_concrete
using ExplicitImports
# ExplicitImports only inspects extensions that are loaded. Load their triggers here so the
# result does not depend on which tests ParallelTestRunner scheduled earlier on this worker.
import QuantumOpticsBase
import QuantumToolbox
import SciMLOperators

# QuantumToolbox ≥ 0.49 re-exports names owned by its `QuantumToolboxCore` submodule; the
# extension reaches them through the public `QuantumToolbox` binding to keep 0.47/0.48 support.
const QTB_OWNER = parentmodule(QuantumToolbox.position)

@testset "Quality gates" begin
    @testset "Aqua" begin
        Aqua.test_all(SecondQuantizedAlgebra)
    end

    @testset "ExplicitImports" begin
        @test check_no_implicit_imports(SecondQuantizedAlgebra) === nothing
        @test check_all_explicit_imports_via_owners(SecondQuantizedAlgebra) === nothing
        @test check_no_stale_explicit_imports(SecondQuantizedAlgebra) === nothing
        @test check_all_qualified_accesses_via_owners(
            SecondQuantizedAlgebra; skip = (Base => Core, QuantumToolbox => QTB_OWNER),
        ) === nothing
        @test check_no_self_qualified_accesses(SecondQuantizedAlgebra) === nothing
    end

    @testset "isbits Op/Index" begin
        # The whole point of the interned-id encoding (#137): Op/Index are isbits,
        # so Vector{Op} stores inline with no GC scanning. Guard the layout.
        @test isbitstype(SecondQuantizedAlgebra.Index)
        @test isbitstype(SecondQuantizedAlgebra.Op)
        @test Base.allocatedinline(SecondQuantizedAlgebra.Op)
        # Pinned to the documented 44 bytes so a layout regression can't slip through.
        @test sizeof(SecondQuantizedAlgebra.Op) == 44
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
end
