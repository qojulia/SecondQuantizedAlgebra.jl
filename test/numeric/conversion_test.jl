using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

dat(x) = dense(x).data

@testset "Numeric conversion" begin
    @testset "Single space basic" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
    end

    @testset "Products" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a' * a, b) isa Operator
        @test to_numeric(a' * a, b; op_type = identity) isa LazySum
        @test dat(to_numeric(a' * a, b)) == dat(create(b) * destroy(b))
        @test dat(to_numeric(2 * a, b)) == dat(2 * destroy(b))
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        @test result isa Operator
        @test to_numeric(a + a', b; op_type = identity) isa LazySum
        @test dat(result) == dat(destroy(b) + create(b))
    end

    @testset "Scalar" begin
        b = FockBasis(7)
        @test to_numeric(3, b) == 3 * one(b)
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2)
        b = FockBasis(3) ⊗ FockBasis(3)

        a1_num = to_numeric(a1, b)
        @test a1_num isa Operator
        @test to_numeric(a1, b; op_type = identity) isa LazySum
    end

    @testset "NLevel numeric" begin
        h = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(h, :σ, 1, 2)
        b = NLevelBasis(3)
        @test to_numeric(σ12, b) == transition(b, 1, 2)
        @test to_numeric(σ12', b) == transition(b, 2, 1)
    end

    @testset "Pauli numeric" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)
        b = SpinBasis(1 // 2)
        @test to_numeric(σx, b) == sigmax(b)
        @test to_numeric(σy, b) == sigmay(b)
        @test to_numeric(σz, b) == sigmaz(b)
    end

    @testset "Spin numeric" begin
        h = SpinSpace(:s)
        Sx = Spin(h, :S, 1)
        Sy = Spin(h, :S, 2)
        Sz = Spin(h, :S, 3)
        b = SpinBasis(5 // 2)
        @test to_numeric(Sx, b) == 0.5 * sigmax(b)
        @test to_numeric(Sy, b) == 0.5 * sigmay(b)
        @test to_numeric(Sz, b) == 0.5 * sigmaz(b)
    end

    @testset "Composite NLevel + Fock" begin
        h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1)
        @qnumbers a::Destroy(h, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        bf = FockBasis(3)
        bn = NLevelBasis(3)
        bc = bf ⊗ bn
        @test to_numeric(σ12, bc) isa Operator
        @test to_numeric(σ12, bc; op_type = identity) isa LazySum
    end

    @testset "Composite basis with gaps" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3, 1)
        hprod_gap = hfock ⊗ hnlevel ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod_gap = bfock ⊗ bnlevel ⊗ bnlevel

        a = Destroy(hprod_gap, :a, 1)
        σprod_gap(i, j) = Transition(hprod_gap, :σ, i, j, 3)

        for i in 1:3, j in 1:3
            i == j == 1 && continue
            op1 = a * σprod_gap(i, j)
            op2 = a' * σprod_gap(i, j)
            ref1 = LazyTensor(
                bprod_gap, [1, 3],
                (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            ref2 = LazyTensor(
                bprod_gap, [1, 3],
                (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            @test dat(to_numeric(op1, bprod_gap)) == dat(ref1)
            @test dat(to_numeric(op2, bprod_gap)) == dat(ref2)
        end
    end

    @testset "Large Hilbert space" begin
        hfock = FockSpace(:fock)
        @qnumbers a::Destroy(hfock)
        bfock = FockBasis(100)

        ref = 2 * create(bfock) + 2 * destroy(bfock)
        got = to_numeric(2 * a + 2 * a', bfock)
        @test isequal(dat(ref), dat(got))
        @test iszero(dat(ref) - dat(got))
        @test isequal(dat(to_numeric(2 * a, bfock)), dat(2 * destroy(bfock)))
    end

    @testset "Round-trip with coherent state" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
        @test dat(to_numeric(a' * a, b)) == dat(create(b) * destroy(b))

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)

        expr = a + a' * a
        @test numeric_average(expr, ψ) ≈ α + abs2(α)

        scalar_qadd = 3 * commutator(a, a')
        @test dense(to_numeric(scalar_qadd, b)) == 3 * one(b)
    end
end
