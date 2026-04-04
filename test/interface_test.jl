using SecondQuantizedAlgebra
using SymbolicUtils
using TermInterface
using Test

@testset "TermInterface" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "QSym is not callable" begin
        @test SymbolicUtils.iscall(a) == false
        @test SymbolicUtils.iscall(ad) == false
    end

    @testset "QMul TermInterface" begin
        m = 3 * a * ad
        @test SymbolicUtils.iscall(m) == true
        @test SymbolicUtils.operation(m) == (*)
        args = SymbolicUtils.arguments(m)
        @test args[1] == 3
        @test length(args) == 3
        @test TermInterface.metadata(m) === nothing
    end

    @testset "QAdd TermInterface" begin
        s = a + ad
        @test SymbolicUtils.iscall(s) == true
        @test SymbolicUtils.operation(s) == (+)
        args = SymbolicUtils.arguments(s)
        @test length(args) == 2
        @test all(x -> x isa QMul, args)
        @test TermInterface.metadata(s) === nothing
    end

    @testset "maketerm QMul" begin
        m = 2 * a * ad
        args = SymbolicUtils.arguments(m)
        m2 = TermInterface.maketerm(typeof(m), *, args, nothing)
        @test isequal(m, m2)
    end

    @testset "maketerm QAdd" begin
        s = a + ad
        args = SymbolicUtils.arguments(s)
        s2 = TermInterface.maketerm(typeof(s), +, args, nothing)
        @test isequal(s, s2)
    end

    @testset "symtype" begin
        @test SymbolicUtils.symtype(a) == Destroy
        @test SymbolicUtils.symtype(ad) == Create
    end

    @testset "Type stability" begin
        m = 2 * a * ad
        s = a + ad

        @inferred SymbolicUtils.iscall(a)
        @inferred SymbolicUtils.iscall(m)
        @inferred SymbolicUtils.iscall(s)
        @inferred SymbolicUtils.operation(m)
        @inferred SymbolicUtils.operation(s)
        @inferred SymbolicUtils.arguments(m)
        @inferred SymbolicUtils.arguments(s)
        @inferred TermInterface.metadata(a)
        @inferred TermInterface.metadata(m)
        @inferred TermInterface.metadata(s)
        @inferred TermInterface.head(a)
        @inferred SymbolicUtils.symtype(a)
        @inferred Base.one(a)
        @inferred Base.zero(a)
        @inferred Base.isone(a)
    end

    @testset "Allocations" begin
        m = 2 * a * ad
        s = a + ad

        # iscall: zero alloc for leaf checks
        @test @allocations(SymbolicUtils.iscall(a)) == 0
        @test @allocations(SymbolicUtils.iscall(m)) == 0
        @test @allocations(SymbolicUtils.iscall(s)) == 0

        # operation: zero alloc
        @test @allocations(SymbolicUtils.operation(m)) == 0
        @test @allocations(SymbolicUtils.operation(s)) == 0

        # metadata: zero alloc
        @test @allocations(TermInterface.metadata(a)) == 0
        @test @allocations(TermInterface.metadata(m)) == 0

        # maketerm roundtrip
        args_m = SymbolicUtils.arguments(m)
        TermInterface.maketerm(typeof(m), *, args_m, nothing)  # warmup
        @test @allocations(TermInterface.maketerm(typeof(m), *, args_m, nothing)) <= 5

        args_s = SymbolicUtils.arguments(s)
        TermInterface.maketerm(typeof(s), +, args_s, nothing)  # warmup
        @test @allocations(TermInterface.maketerm(typeof(s), +, args_s, nothing)) <= 15
    end
end
