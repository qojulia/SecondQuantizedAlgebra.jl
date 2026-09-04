using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QSym
using QuantumOpticsBase
using Test

dat(x) = dense(x).data

@testset "Numeric materialization and performance" begin
    @testset "Type stability" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        ψ = fockstate(b, 2)

        @test @inferred(to_numeric(a, h, 7; backend = QuantumOpticsBackend())) isa Operator
        @test @inferred(to_numeric(a' * a + 2 * a, h, 7; backend = QuantumOpticsBackend())) isa Operator
        hp = SecondQuantizedAlgebra.var"⊗"(FockSpace(:a), FockSpace(:b))
        ap = Destroy(hp, :a, 1)
        @test @inferred(to_numeric(ap, hp, (3, 4); backend = QuantumOpticsBackend())) isa Operator

        @test to_numeric(a, b) isa AbstractOperator
        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(2 * a + 3 * a' + 5, b) isa Operator
        @test to_numeric(a, b, Dict{QSym, Union{}}()) isa AbstractOperator

        @test @inferred(numeric_average(a' * a, ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a) + average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(2 * average(a) * average(a'), ψ)) isa ComplexF64
        @test @inferred(numeric_average(average(a)^2, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3, ψ)) isa ComplexF64
        @test @inferred(numeric_average(3.5 + 1im, ψ)) isa ComplexF64
    end

    @testset "op_type materialization" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(2 * a + 3 * a', b) isa Operator
        @test to_numeric(a, b) isa Operator
        @test to_numeric(a' * a, b; op_type = identity) isa LazySum
        @test to_numeric(2 * a + 3 * a', b; op_type = identity) isa LazySum
        @test to_numeric(a, b) == sparse(to_numeric(a, b; op_type = identity))
        @test dat(to_numeric(a' * a, b; op_type = dense)) == dat(to_numeric(a' * a, b))
    end

    @testset "Allocations: to_numeric" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        to_numeric(a, b)
        to_numeric(a', b)
        to_numeric(a' * a, b)

        @test @allocations(to_numeric(a, b)) < 50
        @test @allocations(to_numeric(a', b)) < 50
        @test @allocations(to_numeric(a' * a, b)) < 1500
    end

    @testset "op_type is shape-independent" begin
        h = SecondQuantizedAlgebra.var"⊗"(FockSpace(:a), FockSpace(:b))
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = QuantumOpticsBase.var"⊗"(FockBasis(3), FockBasis(3))
        exprs = (a1, a1' * a1, a1' * a1 + a2' * a2)

        for expr in exprs
            @test to_numeric(expr, b) isa Operator
            @test to_numeric(expr, b) == to_numeric(expr, b; op_type = sparse)
            @test to_numeric(expr, b; op_type = dense) isa Operator
        end
        @test to_numeric(exprs[1], b; op_type = identity) isa LazySum
        @test to_numeric(exprs[3], b; op_type = identity) isa LazySum
        for expr in exprs
            @test to_numeric(expr, b) == sparse(to_numeric(expr, b; op_type = identity))
            @test dense(to_numeric(expr, b)) ≈ dense(to_numeric(expr, b; op_type = dense))
        end
    end
end
