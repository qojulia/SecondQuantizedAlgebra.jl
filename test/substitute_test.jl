using SecondQuantizedAlgebra
using Test
using Symbolics: Symbolics, @variables
import SecondQuantizedAlgebra: substitute, QMul, QAdd, QSym, CNum, _substitute_qmul, _split_sub_dict

@testset "Substitute" begin
    hf = FockSpace(:c)
    a = Destroy(hf, :a)
    @variables x::Real y::Real

    @testset "Number passthrough" begin
        @test substitute(42, Dict(x => 0)) == 42
        @test substitute(0, Dict(x => 1)) == 0
    end

    @testset "QSym substitution" begin
        @test isequal(substitute(a, Dict(a => a')), a')
        @test isequal(substitute(a, Dict(a' => a)), a)  # no match, unchanged
    end

    @testset "QMul — symbolic variable substitution" begin
        @test iszero(substitute(x * a, Dict(x => 0)))
        @test isequal(substitute(x * a, Dict(x => 1)), a)
        @test isequal(substitute(x * a, Dict(x => y)), y * a)
    end

    @testset "QMul — operator substitution" begin
        @test iszero(substitute(x * a, Dict(a => 0)))
        @test isequal(substitute(x * a, Dict(a => y)), x * y)
        @test isequal(substitute(x * a, Dict(a => a')), x * a')
    end

    @testset "QAdd substitution" begin
        @test isequal(substitute(x * (a + a'), Dict(x => y)), y * (a + a'))
    end

    @testset "Composite space" begin
        hf2 = FockSpace(:c2)
        h = hf ⊗ hf2
        b = Destroy(h, :b, 2)
        expr = x * a' * b
        @test isequal(substitute(expr, Dict(x => 2)), 2 * a' * b)
    end

    @testset "Type stability — internal _substitute_qmul" begin
        d = Dict(x => y)
        sym_dict, op_dict = _split_sub_dict(d)
        m = x * a
        @inferred _substitute_qmul(m, sym_dict, op_dict)
    end
end
