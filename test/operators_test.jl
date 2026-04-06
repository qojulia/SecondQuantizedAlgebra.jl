using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: QMul, QAdd, QSym, QField, _conj, _inconj, _adjoint,
    AvgFunc, _average

@testset "Operators" begin

    @testset "fundamental_operators — FockSpace" begin
        h = FockSpace(:c)
        ops = fundamental_operators(h)
        @test length(ops) == 1
        @test ops[1] isa Destroy
        @test ops[1].name == :a
        @test ops[1].space_index == 1
    end

    @testset "fundamental_operators — NLevelSpace" begin
        # 2-level with ground state 1: only σ₁₂ and σ₂₂
        h = NLevelSpace(:atom, 2, 1)
        ops = fundamental_operators(h)
        @test all(op -> op isa Transition, ops)
        @test length(ops) == 2  # σ₁₂ and σ₂₂ (skips σ₁₁ = ground state projector)

        # 3-level with ground state 1
        h3 = NLevelSpace(:atom, 3, 1)
        ops3 = fundamental_operators(h3)
        # Pairs: (1,2),(1,3),(2,2),(2,3),(3,3) — skip (1,1) = ground
        @test length(ops3) == 5
    end

    @testset "fundamental_operators — PauliSpace" begin
        h = PauliSpace(:p)
        ops = fundamental_operators(h)
        @test length(ops) == 3
        @test all(op -> op isa Pauli, ops)
        @test [op.axis for op in ops] == [1, 2, 3]
    end

    @testset "fundamental_operators — SpinSpace" begin
        h = SpinSpace(:s, 1 // 2)
        ops = fundamental_operators(h)
        @test length(ops) == 3
        @test all(op -> op isa Spin, ops)
        @test [op.axis for op in ops] == [1, 2, 3]
    end

    @testset "fundamental_operators — PhaseSpace" begin
        h = PhaseSpace(:q)
        ops = fundamental_operators(h)
        @test length(ops) == 2
        @test ops[1] isa Position
        @test ops[2] isa Momentum
    end

    @testset "fundamental_operators — ProductSpace" begin
        h = FockSpace(:f) ⊗ NLevelSpace(:n, 2, 1)
        ops = fundamental_operators(h)
        # Fock: 1 (Destroy), NLevel 2-level: 2 (σ₁₂, σ₂₂)
        @test length(ops) == 3
        @test ops[1].space_index == 1  # Fock operator
        @test ops[2].space_index == 2  # NLevel operators
        @test ops[3].space_index == 2
    end

    @testset "fundamental_operators — custom names" begin
        h = FockSpace(:c)
        ops = fundamental_operators(h; names = [:b])
        @test ops[1].name == :b
    end

    @testset "unique_ops" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # a and a† represent the same degree of freedom
        ops = unique_ops([a, ad])
        @test length(ops) == 1

        # Hermitian operator: σx and σx' are the same
        hp = PauliSpace(:p)
        σx = Pauli(hp, :σ, 1)
        @test length(unique_ops([σx, σx'])) == 1

        # Different operators stay
        σy = Pauli(hp, :σ, 2)
        @test length(unique_ops([σx, σy])) == 2
    end

    @testset "find_operators — Fock order 1" begin
        h = FockSpace(:c)
        ops = find_operators(h, 1)
        # Order 1: a and a† but unique_ops removes one
        @test length(ops) == 1
    end

    @testset "find_operators — Fock order 2" begin
        h = FockSpace(:c)
        ops = find_operators(h, 2)
        # Order 1: a (a† removed as adjoint duplicate)
        # Order 2: a*a, a†*a, a†*a† (a*a† = adjoint of a†*a → removed)
        @test length(ops) >= 3
    end

    @testset "find_operators — ProductSpace auto-naming" begin
        # Two FockSpaces → duplicate types → auto-named a, b
        h = FockSpace(:f1) ⊗ FockSpace(:f2)
        ops = find_operators(h, 1)
        names = [op isa QMul ? op.args_nc[1].name : op.name for op in ops]
        @test :a in names
        @test :b in names
    end

    @testset "_conj" begin
        @test _conj(3 + 2im) == 3 - 2im
        @test _conj(5) == 5

        h = FockSpace(:c)
        a = Destroy(h, :a)
        @test _conj(a) == a'  # adjoint of Destroy is Create
    end

    @testset "_inconj" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        ad = a'

        # _inconj on an average: takes adjoint of inner operator
        avg_a = average(a)
        result = _inconj(avg_a)
        @test is_average(result)
        # Inner should be a† (adjoint of a)
        using SymbolicUtils: SymbolicUtils
        inner_wrapped = SymbolicUtils.arguments(result)[1]
        inner = SymbolicUtils.isconst(inner_wrapped) ? inner_wrapped.val : inner_wrapped
        @test inner == ad
    end

    @testset "_adjoint" begin
        h = FockSpace(:c)
        a = Destroy(h, :a)
        @test _adjoint(a) == a'
        @test _adjoint(3 + 2im) == 3 - 2im
    end

end
