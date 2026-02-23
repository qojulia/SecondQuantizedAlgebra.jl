using SecondQuantizedAlgebra
using SymbolicUtils
using QuantumOpticsBase
using Test
using Random;
Random.seed!(0)

@testset "numeric-conversion" begin
    @testset "Basic Fock Space Conversion" begin
        hfock = FockSpace(:fock)
        @qnumbers a::Destroy(hfock)
        bfock = FockBasis(7)

        @test to_numeric(a, bfock) == destroy(bfock)
        @test to_numeric(a', bfock) == create(bfock)
    end

    @testset "NLevel Space Conversion" begin
        @testset "Integer Levels" begin
            hnlevel = NLevelSpace(:nlevel, 3)
            σ(i, j) = Transition(hnlevel, :σ, i, j)
            bnlevel = NLevelBasis(3)

            for i in 1:3, j in 1:3
                op = σ(i, j)
                @test to_numeric(op, bnlevel) == QuantumOpticsBase.transition(bnlevel, i, j)
                @test to_numeric(op, bnlevel, Dict()) == QuantumOpticsBase.transition(bnlevel, i, j)
            end
        end

        @testset "Symbolic Levels" begin
            levels = (:g, :e, :a)
            hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
            σ_sym(i, j) = Transition(hnlevel_sym, :σ, i, j)
            bnlevel = NLevelBasis(3)

            @test_throws ArgumentError to_numeric(σ_sym(:e, :g), bnlevel)

            level_map = Dict((levels .=> (1, 2, 3))...)
            for i in 1:3, j in 1:3
                lvl1 = levels[i]
                lvl2 = levels[j]
                op = σ_sym(lvl1, lvl2)
                @test to_numeric(op, bnlevel; level_map=level_map) ==
                    QuantumOpticsBase.transition(bnlevel, i, j)
            end
        end
        @testset "Pauli Operators" begin
            hp = PauliSpace(:p)
            σx = Pauli(hp, :σ, 1)
            σy = Pauli(hp, :σ, 2)
            σz = Pauli(hp, :σ, 3)

            bp = SpinBasis(1//2)
            @test to_numeric(σx, bp) == sigmax(bp)
            @test to_numeric(σy, bp) == sigmay(bp)
            @test to_numeric(σz, bp) == sigmaz(bp)
        end
        @testset "Spin Operators" begin
            hp = SpinSpace(:p)
            Sx = Spin(hp, :S, 1)
            Sy = Spin(hp, :S, 2)
            Sz = Spin(hp, :S, 3)

            bs = SpinBasis(5//2)
            @test to_numeric(Sx, bs) == 0.5*sigmax(bs)
            @test to_numeric(Sy, bs) == 0.5*sigmay(bs)
            @test to_numeric(Sz, bs) == 0.5*sigmaz(bs)
        end
        @testset "PhaseSpace" begin
            h = PhaseSpace(:motion)
            x = Position(h, :x, 1)
            p = Momentum(h, :p, 1)

            xmin = -2.0
            xmax = 4.0
            N = 100
            b_position = PositionBasis(xmin, xmax, N)
            x_qo = to_numeric(x, b_position)
            p_qo = to_numeric(p, b_position)

            x0 = 1.2;
            p0 = 0.4;
            sigma = 0.2
            ψ0 = gaussianstate(b_position, x0, p0, sigma)

            x_qo_ = position(b_position)
            p_qo_ = momentum(b_position)
            @test isequal(dense(x_qo), dense(x_qo_))
            @test isequal(dense(p_qo), dense(p_qo_))

            @test abs(expect(x_qo, ψ0) - x0) < 1e-6
            @test abs(expect(p_qo, ψ0) - p0) < 1e-6
            @test abs(√(variance(x_qo, ψ0)*2) - sigma) < 1e-6

            @test numeric_average(x, ψ0) ≈ x0
            @test numeric_average(p, ψ0) ≈ p0
            @test abs(numeric_average((x-x0)^2, ψ0) - variance(x_qo, ψ0)) < 1e-6
        end
    end

    @testset "Composite Basis Conversion" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3)
        hprod = hfock ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel

        a = Destroy(hprod, :a)
        σprod(i, j) = Transition(hprod, :σ, i, j)

        @testset "Regular Product Operations" begin
            for i in 1:3, j in 1:3
                i == j == 1 && continue  # rewritten as sum, see below
                op1 = a * σprod(i, j)
                op2 = a' * σprod(i, j)
                @test to_numeric(op1, bprod) == LazyTensor(
                    bprod,
                    [1, 2],
                    (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
                )
                @test to_numeric(op2, bprod) == LazyTensor(
                    bprod,
                    [1, 2],
                    (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
                )
            end

            @test to_numeric(a' * a, bprod) ==
                LazyTensor(bprod, [1], (create(bfock) * destroy(bfock),))
        end

        @testset "Special Cases (LazySum)" begin
            op1_num = to_numeric(a * σprod(1, 1), bprod)
            @test op1_num isa LazySum
            @test sparse(op1_num) ==
                destroy(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

            op2_num = to_numeric(a' * σprod(1, 1), bprod)
            @test op2_num isa LazySum
            @test sparse(op2_num) ==
                create(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)
        end
    end

    @testset "Composite Basis with Symbolic Levels" begin
        hfock = FockSpace(:fock)
        levels = (:g, :e, :a)
        hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod = bfock ⊗ bnlevel
        level_map = Dict((levels .=> (1, 2, 3))...)

        σsym_prod(i, j) = Transition(hfock ⊗ hnlevel_sym, :σ, i, j)
        a = Destroy(hfock ⊗ hnlevel_sym, :a)

        @test_throws ArgumentError to_numeric(a * σsym_prod(:e, :g), bprod)

        for i in 1:3, j in 1:3
            i == j == 1 && continue  # see below
            op1 = a * σsym_prod(levels[i], levels[j])
            op2 = a' * σsym_prod(levels[i], levels[j])
            @test to_numeric(op1, bprod; level_map=level_map) == LazyTensor(
                bprod, [1, 2], (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j))
            )
            @test to_numeric(op2, bprod; level_map=level_map) == LazyTensor(
                bprod, [1, 2], (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j))
            )
        end

        op1_num = to_numeric(a * σsym_prod(:g, :g), bprod; level_map=level_map)
        @test op1_num isa LazySum
        @test sparse(op1_num) ==
            destroy(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

        op2_num = to_numeric(a' * σsym_prod(:g, :g), bprod; level_map=level_map)
        @test op2_num isa LazySum
        @test sparse(op2_num) == create(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)
    end

    @testset "Composite Basis with Gaps" begin
        hfock = FockSpace(:fock)
        hnlevel = NLevelSpace(:nlevel, 3)
        hprod_gap = hfock ⊗ hnlevel ⊗ hnlevel
        bfock = FockBasis(7)
        bnlevel = NLevelBasis(3)
        bprod_gap = bfock ⊗ bnlevel ⊗ bnlevel

        a = Destroy(hprod_gap, :a)
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

        op1_num = to_numeric(a * σprod_gap(1, 1), bprod_gap)
        @test op1_num isa LazySum
        @test sparse(op1_num) ==
            destroy(bfock) ⊗ one(bnlevel) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

        op2_num = to_numeric(a' * σprod_gap(1, 1), bprod_gap)
        @test op2_num isa LazySum
        @test sparse(op2_num) ==
            create(bfock) ⊗ one(bnlevel) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)
    end

    @testset "Numeric Average Values" begin
        @testset "Basic Fock State Averages" begin
            hfock = FockSpace(:fock)
            bfock = FockBasis(7)
            α = 0.1 + 0.2im
            ψ = coherentstate(bfock, α)
            a = Destroy(hfock, :a)

            @test numeric_average(a, ψ) ≈ numeric_average(average(a), ψ) ≈ α
            @test numeric_average(a' * a, ψ) ≈ numeric_average(average(a' * a), ψ) ≈ abs2(α)

            @test numeric_average(a + a'a, ψ) ≈
                numeric_average(average(a) + average(a'a), ψ) ≈
                numeric_average(average(a + a'a), ψ) ≈
                α + abs2(α)
            @test numeric_average(average(a) * average(a'a), ψ) ≈ α * α' * α

            @test numeric_average(average(a)^2, ψ) ≈ α^2
            @test numeric_average(average(a) * average(a'a) + 3 * average(a)^2 + 0.6, ψ) ≈
                α * α' * α + 3α^2 + 0.6
            @test numeric_average(3, ψ) ≈ 3
        end

        @testset "Product State Averages" begin
            hfock = FockSpace(:fock)
            hnlevel = NLevelSpace(:nlevel, 3)
            levels = (:g, :e, :a)
            hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
            bfock = FockBasis(7)
            bnlevel = NLevelBasis(3)
            bprod = bfock ⊗ bnlevel
            level_map = Dict((levels .=> (1, 2, 3))...)

            α = 0.1 + 0.2im
            ψ = coherentstate(bfock, α)
            ψprod = ψ ⊗ nlevelstate(bnlevel, 1)

            σprod(i, j) = Transition(hfock ⊗ hnlevel, :σ, i, j)
            σsym_prod(i, j) = Transition(hfock ⊗ hnlevel_sym, :σ, i, j)

            @test_throws ArgumentError numeric_average(σsym_prod(:e, :g), ψprod)

            idfock = one(bfock)
            for i in 1:3, j in 1:3
                op = σprod(i, j)
                op_sym = σsym_prod(levels[i], levels[j])
                op_num = idfock ⊗ QuantumOpticsBase.transition(bnlevel, i, j)
                @test numeric_average(op, ψprod) ≈ expect(op_num, ψprod)
                @test numeric_average(op_sym, ψprod; level_map=level_map) ≈
                    expect(op_num, ψprod)
            end
        end

        @testset "LazyKet Support" begin
            if isdefined(QuantumOpticsBase, :LazyKet)
                hfock = FockSpace(:fock)
                hnlevel = NLevelSpace(:nlevel, 3)
                levels = (:g, :e, :a)
                hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
                bfock = FockBasis(7)
                bnlevel = NLevelBasis(3)
                bprod = bfock ⊗ bnlevel
                level_map = Dict((levels .=> (1, 2, 3))...)

                α = 0.1 + 0.2im
                ψ = coherentstate(bfock, α)
                ψlazy = LazyKet(bprod, (ψ, nlevelstate(bnlevel, 1)))

                σprod(i, j) = Transition(hfock ⊗ hnlevel, :σ, i, j)
                σsym_prod(i, j) = Transition(hfock ⊗ hnlevel_sym, :σ, i, j)

                @test_throws ArgumentError numeric_average(σsym_prod(:e, :g), ψlazy)
                for i in 1:3, j in 1:3
                    op = σprod(i, j)
                    op_sym = σsym_prod(levels[i], levels[j])
                    op_num = LazyTensor(
                        bprod, [2], (QuantumOpticsBase.transition(bnlevel, i, j),)
                    )
                    @test numeric_average(op, ψlazy) ≈ expect(op_num, ψlazy)
                    @test numeric_average(op_sym, ψlazy; level_map=level_map) ≈
                        expect(op_num, ψlazy)
                end
            end
        end
    end

    @testset "to_numeric and numeric_average with dictionary" begin
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

        q_add = 3ad + a # test for QAdd metadata
        @test isequal(substitute(q_add, dd), 3*s(2, 2, 1) + s(2, 1, 2))

        @test to_numeric(a, b_test, Dict()) == to_numeric(a, b_test)
        @test to_numeric(ad*n, b_test, Dict()) == to_numeric(ad*n, b_test)
        @test to_numeric(2*ad*a, b_test, Dict()) == to_numeric(2*ad*a, b_test)
        @test to_numeric(ad, b_all, dd) == s(2, 2, 1)
        @test to_numeric(2 * ad * a, b_all, dd) == 2 * s(2, 2, 1) * s(2, 1, 2)
        @test dense(to_numeric(3, b_all, dd)) == one(b_all) * 3

        ψ0 = tensor([nlevelstate(bs, 2) for i in 1:nQDs]...)
        @test numeric_average(average(ad * a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0)
        @test numeric_average(average(ad) * average(ad * a) + 3, ψ0, dd) ==
            expect(s(2, 2, 1), ψ0) * expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + 3
        @test numeric_average(3 * average(ad)^2, ψ0, dd) == 3 * expect(s(2, 2, 1), ψ0)^2
        @test numeric_average(average(ad * a) + average(a), ψ0, dd) ==
            expect(s(2, 2, 1) * s(2, 1, 2), ψ0) + expect(s(2, 1, 2), ψ0)
    end

    @testset "Edge Cases and Bug Fixes" begin
        @testset "Large Hilbert Space (Issue #109)" begin
            hfock = FockSpace(:fock)
            @qnumbers a::Destroy(hfock)
            bfock = FockBasis(100)

            diff =
                (2 * create(bfock) + 2 * destroy(bfock)) -
                to_numeric((2 * (a) + 2 * (a')), bfock)
            @test isequal(
                2 * create(bfock) + 2 * destroy(bfock),
                to_numeric((2 * (a) + 2 * (a')), bfock),
            )
            @test iszero(diff)

            @test isequal(to_numeric(2 * a, bfock), 2 * to_numeric(a, bfock))
            @test iszero(to_numeric(2 * a, bfock) - 2 * to_numeric(a, bfock))
        end

        @testset "Indexed Initial State (Superradiant Pulse)" begin
            order = 2 #order of the cumulant expansion
            @cnumbers κ g Γ Δ N
            hc = FockSpace(:cavity)
            ha = NLevelSpace(:atom, 2)
            h = hc ⊗ ha
            i = Index(h, :i, N, ha)
            j = Index(h, :j, N, ha)
            k = Index(h, :k, N, ha)
            @qnumbers a::Destroy(h, 1)
            σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

            bc = FockBasis(3)
            basis_a = NLevelBasis(2)
            b = tensor(bc, [basis_a for i in 1:order]...)
            ψc = fockstate(bc, 0)
            ψa = normalize(nlevelstate(basis_a, 1) + nlevelstate(basis_a, 2))
            ψ = tensor(ψc, [ψa for i in 1:order]...)

            a_ = LazyTensor(b, [1], (destroy(bc),))
            σ_(i, j, k) = LazyTensor(
                b, [1 + k], (QuantumOpticsBase.transition(basis_a, i, j),)
            )
            ranges = [1, 2]

            @test to_numeric(σ(1, 2, 1), b; ranges=ranges) == σ_(1, 2, 1)
            @test to_numeric(σ(2, 2, 2), b; ranges=ranges) == σ_(2, 2, 2)
            @test to_numeric(a, b; ranges=ranges) == a_
            @test to_numeric(a * σ(2, 2, 2), b; ranges=ranges) == σ_(2, 2, 2) * a_
            @test numeric_average(σ(2, 2, 2), ψ; ranges=ranges) ≈ 0.5
            @test numeric_average(average(σ(2, 2, 1)), ψ; ranges=ranges) ≈ 0.5
            @test numeric_average(average(a'a), ψ; ranges=ranges) ≈ 0.0
            @test numeric_average(average(a * σ(2, 2, 1)), ψ; ranges=ranges) ≈ 0.0
            @test_throws ArgumentError numeric_average(average(a'a), ψ)

            if isdefined(QuantumOpticsBase, :LazyKet)
                ψlazy = LazyKet(b, (ψc, (ψa for i in 1:order)...))
                @test numeric_average(σ(2, 2, 2), ψlazy; ranges=ranges) ≈ 0.5
                @test numeric_average(average(σ(2, 2, 1)), ψlazy; ranges=ranges) ≈ 0.5
                @test numeric_average(average(a'a), ψlazy; ranges=ranges) ≈ 0.0
                @test numeric_average(average(a * σ(2, 2, 1)), ψlazy; ranges=ranges) ≈ 0.0
                @test_throws ArgumentError numeric_average(average(a'a), ψlazy)
            end
        end
    end
end
