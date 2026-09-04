using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: QSym
using QuantumOpticsBase
using Test

dat(x) = dense(x).data

@testset "Numeric expectation values" begin
    @testset "numeric_average" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)
    end

    @testset "numeric_average: Average expressions" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        avg_a = average(a)
        @test numeric_average(avg_a, ψ) ≈ α

        avg_sum = average(a) + average(a')
        @test numeric_average(avg_sum, ψ) ≈ α + conj(α)

        avg_scaled = 2 * average(a)
        @test numeric_average(avg_scaled, ψ) ≈ 2α
    end

    @testset "numeric_average: Dict substitution" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)

        d = Dict{QSym, Any}()
        @test numeric_average(a, ψ, d) ≈ α

        avg_a = average(a)
        @test numeric_average(avg_a, ψ, d) ≈ α

        @test numeric_average(3, ψ, d) === ComplexF64(3)
    end

    @testset "numeric_average: product state" begin
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

    @testset "numeric_average: comprehensive expressions" begin
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

    @testset "numeric_average: vector of states" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im, 0.5 + 0.0im)
        ψs = [coherentstate(b, α) for α in αs]

        @test numeric_average(a, ψs) ≈ [α for α in αs]
        @test numeric_average(a' * a, ψs) ≈ [abs2(α) for α in αs]
        @test numeric_average(average(a), ψs) ≈ [α for α in αs]

        @test expect(a, ψs) ≈ numeric_average(a, ψs)
        @test expect(a' * a, ψs[1]) ≈ numeric_average(a' * a, ψs[1])
        @test numeric_average(average(a), ψs[1]) ≈ αs[1]

        empty_ψs = typeof(ψs[1])[]
        @test_throws ArgumentError numeric_average(a, empty_ψs)
    end

    @testset "numeric_average: vector of states with Dict" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        αs = (0.1 + 0.2im, -0.3 + 0.4im)
        ψs = [coherentstate(b, α) for α in αs]
        d = Dict{QSym, Any}()

        @test numeric_average(a, ψs, d) ≈ numeric_average(a, ψs)
        @test numeric_average(average(a' * a), ψs, d) ≈ [abs2(α) for α in αs]
        @test expect(a, ψs, d) ≈ numeric_average(a, ψs, d)
    end

    @testset "numeric_average: Dict comprehensive" begin
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

        @test to_numeric(a, b_test, Dict{QSym, Any}()) == to_numeric(a, b_test)
        @test dat(to_numeric(ad * n, b_test, Dict{QSym, Any}())) == dat(to_numeric(ad * n, b_test))
        @test dat(to_numeric(2 * ad * a, b_test, Dict{QSym, Any}())) == dat(to_numeric(2 * ad * a, b_test))
        @test to_numeric(ad, b_all, dd) == s(2, 2, 1)
        @test dat(to_numeric(2 * ad * a, b_all, dd)) == dat(2 * s(2, 2, 1) * s(2, 1, 2))
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

    @testset "numeric_average — LazyKet state" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, (:a, :b, :c))
        hprod = hfock ⊗ hnlevel
        bfock = FockBasis(10)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel

        @qnumbers a::Destroy(hprod, 1)
        σ(i, j) = Transition(hprod, :σ, i, j, 2)

        α = 0.3 + 0.0im
        ket_fock = coherentstate(bfock, α)
        ket_nlevel = (nlevelstate(bnlevel, 1) + nlevelstate(bnlevel, 3)) / sqrt(2)
        ψ_lazy = LazyKet(bprod, (ket_fock, ket_nlevel))
        ψ_dense = ket_fock ⊗ ket_nlevel

        for op in (a, a' * σ(:a, :c), a + a' * σ(:a, :c))
            @test numeric_average(op, ψ_lazy) ≈ expect(to_numeric(op, bprod), ψ_dense)
        end
        @test numeric_average(average(a' * a), ψ_lazy) ≈
            expect(to_numeric(a' * a, bprod), ψ_dense)
    end

    @testset "numeric_average — unsupported symbolic operation" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(5)
        ψ = fockstate(b, 2)
        @test_throws ArgumentError numeric_average(sqrt(average(a' * a)), ψ)
        @test_throws ArgumentError numeric_average(average(a' * a) / (1 + average(a)), ψ)
    end
end
