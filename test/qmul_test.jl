using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QMul, QSym
using Test

@testset "QMul" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Construction" begin
        m = a * ad
        @test m isa QMul
        @test m.arg_c == 1
        @test length(m.args_nc) == 2
        # Same space: order preserved (a then a†)
        @test m.args_nc[1] isa Destroy
        @test m.args_nc[2] isa Create
    end

    @testset "Stable sort: same-site order preserved" begin
        # a·a† and a†·a must produce distinct QMul (not commuted)
        m1 = a * ad   # [Destroy, Create]
        m2 = ad * a   # [Create, Destroy]
        @test m1.args_nc[1] isa Destroy
        @test m1.args_nc[2] isa Create
        @test m2.args_nc[1] isa Create
        @test m2.args_nc[2] isa Destroy
        @test !isequal(m1, m2)

        # Multi-operator: a·a†·a must keep that order within the site
        m3 = a * ad * a
        @test m3.args_nc[1] isa Destroy
        @test m3.args_nc[2] isa Create
        @test m3.args_nc[3] isa Destroy

        # Cross-site operators reorder freely, but within-site order is preserved
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        # a2 * a1' * a1: cross-site sort groups space-1 ops together,
        # but a1' must stay left of a1 within space 1
        m4 = a2 * a1' * a1
        space1_ops = filter(op -> op.space_index == 1, m4.args_nc)
        @test space1_ops[1] isa Create
        @test space1_ops[2] isa Destroy
    end

    @testset "Scalar multiplication" begin
        m1 = 3 * a
        @test m1 isa QMul
        @test m1.arg_c == 3
        @test m1.args_nc == [a]

        m2 = a * 2.0
        @test m2 isa QMul
        @test m2.arg_c == 2.0

        m3 = 0 * a
        @test m3 isa QMul
        @test m3.arg_c == 0
    end

    @testset "QSym * QSym" begin
        h2 = FockSpace(:a) ⊗ FockSpace(:b)
        a1 = Destroy(h2, :a, 1)
        a2 = Destroy(h2, :b, 2)
        m = a1 * a2
        @test m isa QMul
        @test m.args_nc[1].space_index <= m.args_nc[2].space_index
    end

    @testset "QMul * QSym and QSym * QMul" begin
        m1 = (2 * a) * ad
        @test m1 isa QMul
        @test m1.arg_c == 2
        @test length(m1.args_nc) == 2

        m2 = ad * (3 * a)
        @test m2 isa QMul
        @test m2.arg_c == 3
        @test length(m2.args_nc) == 2
    end

    @testset "QMul * QMul" begin
        m1 = 2 * a
        m2 = 3 * ad
        m3 = m1 * m2
        @test m3 isa QMul
        @test m3.arg_c == 6
        @test length(m3.args_nc) == 2
    end

    @testset "QMul * Number" begin
        m = (2 * a) * 3
        @test m isa QMul
        @test m.arg_c == 6
    end

    @testset "Division" begin
        m = a / 2
        @test m isa QMul
        @test m.arg_c == 1 // 2
    end

    @testset "Power" begin
        m = a^3
        @test m isa QMul
        @test length(m.args_nc) == 3
    end

    @testset "Negation" begin
        m = -a
        @test m isa QMul
        @test m.arg_c == -1
    end

    @testset "Lazy — no commutation" begin
        m = a * ad
        @test m isa QMul
        @test length(m.args_nc) == 2
    end

    @testset "Equality and hashing" begin
        m1 = 2 * a * ad
        m2 = 2 * a * ad
        @test isequal(m1, m2)
        @test hash(m1) == hash(m2)
    end

    @testset "iszero" begin
        m = 0 * a
        @test iszero(m)
        @test !iszero(2 * a)
    end

    @testset "Adjoint" begin
        m = 2 * ad * a
        ma = m'
        @test ma isa QMul
        @test ma.arg_c == 2
        @test ma.args_nc[1] isa Create
        @test ma.args_nc[2] isa Destroy
    end

    @testset "Concrete types" begin
        using CheckConcreteStructs
        all_concrete(QMul)
    end

    @testset "Type stability" begin
        # QSym * QSym
        @inferred a * ad
        # Number * QSym / QSym * Number
        @inferred 3 * a
        @inferred a * 2.0
        # QSym * QMul / QMul * QSym
        @inferred (2 * a) * ad
        @inferred ad * (2 * a)
        # QMul * QMul
        @inferred (2 * a) * (3 * ad)
        # QMul * Number / Number * QMul
        @inferred (2 * a) * 3
        @inferred 3 * (2 * a)
        # Division
        @inferred a / 2
        @inferred a / 2.0
        @inferred (2 * a) / 3
        @inferred (2 * a) / 3.0
        # Power
        @inferred a^3
        @inferred (2 * a)^3
        # Negation
        @inferred -a
        @inferred -(2 * a)
        # Adjoint
        @inferred adjoint(2 * ad * a)
        # Equality / hashing
        @inferred isequal(2 * a, 2 * a)
        @inferred hash(2 * a, UInt(0))
        # iszero / zero
        @inferred iszero(0 * a)
        @inferred zero(2 * a)
    end

    @testset "Allocations" begin
        m_aa = a * ad
        m_2a = 2 * a
        m_3ad = 3 * ad

        # Warmup all code paths
        a * ad; 3 * a; a * 2; m_2a * 3; 3 * m_2a
        a * m_3ad; m_2a * ad; m_2a * m_3ad
        a^3; m_2a^3; -a; -m_2a; adjoint(m_aa)
        m_2a_copy = 2 * a
        isequal(m_2a, m_2a_copy); hash(m_2a, UInt(0))

        # QSym * QSym
        @test @allocations(a * ad) <= 1000
        # Number * QSym
        @test @allocations(3 * a) <= 1000
        @test @allocations(a * 2) <= 1000
        # QMul * Number
        @test @allocations(m_2a * 3) <= 250
        @test @allocations(3 * m_2a) <= 250
        # QSym * QMul
        @test @allocations(a * m_3ad) <= 1000
        @test @allocations(m_2a * ad) <= 1000
        # QMul * QMul
        @test @allocations(m_2a * m_3ad) <= 250
        # Power
        @test @allocations(a^3) <= 1000
        @test @allocations(m_2a^3) <= 1000
        # Negation
        @test @allocations(-a) <= 1000
        @test @allocations(-m_2a) <= 250
        # Adjoint
        @test @allocations(adjoint(m_aa)) <= 1000
        # Equality / hashing
        @test @allocations(isequal(m_2a, m_2a_copy)) <= 5
        @test @allocations(hash(m_2a, UInt(0))) <= 5
    end
end
