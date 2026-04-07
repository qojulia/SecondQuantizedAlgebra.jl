using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QAdd, QSym
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

    @testset "QAdd from product — TermInterface" begin
        # With eager ordering, 3 * a * a† = 3(1 + a†a) = 3 + 3a†a
        m = 3 * a * ad
        @test m isa QAdd
        @test SymbolicUtils.iscall(m) == true
        @test SymbolicUtils.operation(m) == (+)
        args = SymbolicUtils.arguments(m)
        @test length(args) == 2  # scalar 3 + 3*a†a
        @test all(x -> x isa QAdd, args)
        @test TermInterface.metadata(m) === nothing
    end

    @testset "QAdd TermInterface" begin
        s = a + ad
        @test SymbolicUtils.iscall(s) == true
        @test SymbolicUtils.operation(s) == (+)
        args = SymbolicUtils.arguments(s)
        @test length(args) == 2
        @test all(x -> x isa QAdd, args)
        @test TermInterface.metadata(s) === nothing
    end

    @testset "maketerm QAdd from product" begin
        m = 3 * ad * a  # single-term QAdd (already normal ordered)
        args = SymbolicUtils.arguments(m)
        m2 = TermInterface.maketerm(typeof(m), +, args, nothing)
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
        m = ad * a  # QAdd
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

    @static if VERSION >= v"1.12"
        @testset "Allocations" begin
            m = ad * a  # QAdd
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
            args_s = SymbolicUtils.arguments(s)
            TermInterface.maketerm(typeof(s), +, args_s, nothing)  # warmup
            @test @allocations(TermInterface.maketerm(typeof(s), +, args_s, nothing)) <= 300
        end
    end
end
