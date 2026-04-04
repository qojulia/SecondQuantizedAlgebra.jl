using Test
using CheckConcreteStructs: all_concrete
using SecondQuantizedAlgebra

const _CONCRETE_SKIP = Set{Symbol}()

@testset "CheckConcreteStructs" begin
    for name in names(SecondQuantizedAlgebra; all=true)
        name in _CONCRETE_SKIP && continue
        isdefined(SecondQuantizedAlgebra, name) || continue
        T = getfield(SecondQuantizedAlgebra, name)
        T isa Type || continue
        isabstracttype(T) && continue
        T isa UnionAll && continue
        isstructtype(T) || continue
        parentmodule(T) === SecondQuantizedAlgebra || continue
        @testset "$name" begin
            @test all_concrete(T; verbose=false)
        end
    end
end
