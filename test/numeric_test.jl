using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QMul, QAdd
using QuantumOpticsBase
using Test

@testset "numeric conversion" begin
    @testset "Single space — basic" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
    end

    @testset "QMul" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a' * a, b) == create(b) * destroy(b)
        @test to_numeric(2 * a, b) == 2 * destroy(b)
    end

    @testset "QAdd" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        result = to_numeric(a + a', b)
        expected = destroy(b) + create(b)
        @test result == expected
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
        @test a1_num isa LazyTensor
    end

    @testset "numeric_average" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
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
        h = SpinSpace(:s, 5 // 2)
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
        @test to_numeric(σ12, bc) isa LazyTensor
    end

    # TODO: to_numeric is not type-stable due to QuantumOpticsBase dispatch.
    # Fix upstream in QuantumOpticsBase.
    # @testset "Type stability" begin
    #     h = FockSpace(:fock)
    #     @qnumbers a::Destroy(h)
    #     b = FockBasis(7)
    #     @inferred to_numeric(a, b)
    #     @inferred to_numeric(a', b)
    # end

    @testset "numeric_average — Average expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        # numeric_average on average(a) should give ⟨a⟩ = α
        avg_a = average(a)
        @test numeric_average(avg_a, ψ) ≈ α

        # Sum of averages
        avg_sum = average(a) + average(a')
        @test numeric_average(avg_sum, ψ) ≈ α + conj(α)

        # Scalar times average
        avg_scaled = 2 * average(a)
        @test numeric_average(avg_scaled, ψ) ≈ 2α
    end

    @testset "to_numeric — Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        custom_op = 2 * destroy(b)
        d = Dict(a => custom_op)

        # Dict substitution replaces operator
        @test to_numeric(a, b, d) == custom_op
        # Adjoint not in dict → falls back to normal
        @test to_numeric(a', b, d) == create(b)

        # QMul with Dict
        result_mul = to_numeric(a' * a, b, d)
        expected_mul = create(b) * custom_op
        @test result_mul == expected_mul
    end

    @testset "numeric_average — Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        d = Dict{Any, Any}()  # empty dict — same as no dict
        @test numeric_average(a, ψ, d) ≈ α

        # Average expression with dict
        avg_a = average(a)
        @test numeric_average(avg_a, ψ, d) ≈ α

        # Number passthrough
        @test numeric_average(3, ψ, d) === 3
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
            @test to_numeric(op1, bprod_gap) == LazyTensor(
                bprod_gap,
                [1, 3],
                (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
            @test to_numeric(op2, bprod_gap) == LazyTensor(
                bprod_gap,
                [1, 3],
                (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
            )
        end
    end

    @testset "Large Hilbert space" begin
        hfock = FockSpace(:fock)
        @qnumbers a::Destroy(hfock)
        bfock = FockBasis(100)

        @test isequal(
            2 * create(bfock) + 2 * destroy(bfock),
            to_numeric(2 * a + 2 * a', bfock),
        )
        @test iszero(
            (2 * create(bfock) + 2 * destroy(bfock)) -
            to_numeric(2 * a + 2 * a', bfock),
        )
        @test isequal(to_numeric(2 * a, bfock), 2 * to_numeric(a, bfock))
    end

    @testset "numeric_average — product state" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3, 1)
        hprod = hfock ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel

        α = 0.1 + 0.2im
        ψ = coherentstate(bfock, α)
        ψprod = ψ ⊗ nlevelstate(bnlevel, 1)

        σprod(i, j) = Transition(hprod, :σ, i, j, 2)

        idfock = one(bfock)
        for i in 1:3, j in 1:3
            op = σprod(i, j)
            op_num = idfock ⊗ QuantumOpticsBase.transition(bnlevel, i, j)
            @test numeric_average(op, ψprod) ≈ expect(op_num, ψprod)
        end
    end

    @testset "numeric_average — comprehensive expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a + a'a, ψ) ≈ α + abs2(α)
        @test numeric_average(average(a) + average(a'a), ψ) ≈ α + abs2(α)
        @test numeric_average(average(a + a'a), ψ) ≈ α + abs2(α)
        @test numeric_average(average(a) * average(a'a), ψ) ≈ α * α' * α
        @test numeric_average(average(a)^2, ψ) ≈ α^2
        @test numeric_average(3, ψ) ≈ 3
    end

    @testset "numeric_average — Dict comprehensive" begin
        nQDs = 2
        h_qc1 = FockSpace(:ada)
        h_qc2 = FockSpace(:n)
        h_qc = h_qc1 ⊗ h_qc2
        a = Destroy(h_qc, :a, 1)
        n = Destroy(h_qc, :n, 2)
        ad = a'

        bs = NLevelBasis(2)
        b_all = tensor([bs for i in 1:nQDs]...)
        s(α, i, j) = embed(b_all, α, transition(bs, i, j))
        b_test = FockBasis(2) ⊗ FockBasis(3)

        dd = Dict([ad, a] .=> [s(2, 2, 1), s(2, 1, 2)])

        @test to_numeric(a, b_test, Dict()) == to_numeric(a, b_test)
        @test to_numeric(ad * n, b_test, Dict()) == to_numeric(ad * n, b_test)
        @test to_numeric(2 * ad * a, b_test, Dict()) == to_numeric(2 * ad * a, b_test)
        @test to_numeric(ad, b_all, dd) == s(2, 2, 1)
        @test to_numeric(2 * ad * a, b_all, dd) == 2 * s(2, 2, 1) * s(2, 1, 2)
        @test dense(to_numeric(3, b_all, dd)) == one(b_all) * 3

        ψ0 = tensor([nlevelstate(bs, 2) for i in 1:nQDs]...)
        @test numeric_average(average(ad * a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0)
        @test numeric_average(average(ad) * average(ad * a) + 3, ψ0, dd) ==
            expect(s(2, 2, 1), ψ0) * expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + 3
        @test numeric_average(3 * average(ad)^2, ψ0, dd) ==
            3 * expect(s(2, 2, 1), ψ0)^2
        @test numeric_average(average(ad * a) + average(a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + expect(s(2, 1, 2), ψ0)
    end

    @testset "Allocations — to_numeric" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        # Warmup
        to_numeric(a, b)
        to_numeric(a', b)
        to_numeric(a' * a, b)

        # Single operators should be bounded
        @test @allocations(to_numeric(a, b)) < 50
        @test @allocations(to_numeric(a', b)) < 50
        @test @allocations(to_numeric(a' * a, b)) < 1500
    end
end
