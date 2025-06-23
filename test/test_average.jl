using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "average" begin
    # Common setup
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha

    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @testset "Basic Average Properties" begin
        @test isequal(average(2a), 2 * average(a))
        @test isequal(average(2 * (a + a)), 2 * (average(a) + average(a)))
        @test isequal(average(a^2), average(a * a))
    end

    @testset "Average Arithmetic" begin
        @test isequal(simplify(average(σ) + average(σ)), average(2σ))
        @test iszero(simplify(average(σ) - average(σ)))
    end

    @testset "C-Number Handling" begin
        ωc, ωa = cnumbers("ω_c ω_a")
        @test isequal(average(ωc), ωc)
        @test isequal(average(ωc * a), ωc * average(a))
        @test isequal(average(ωc * (a + a')), ωc * average(a) + ωc * average(a'))
    end

    @testset "Commutativity Properties" begin
        @test iszero(average(a * σ) * average(a) - average(a) * average(a * σ))
    end

    @testset "Double Average" begin
        # Reference: https://github.com/qojulia/QuantumCumulants.jl/issues/242
        @cnumbers Δ_ g κ η

        # Multi-level system setup
        hf_complex = FockSpace(:cavity)
        ha1 = NLevelSpace(:atom1, 2)
        ha2 = NLevelSpace(:atom2, 2)
        h_complex = hf_complex ⊗ ha1 ⊗ ha2

        a_complex = Destroy(h_complex, :a)
        s1(i, j) = Transition(h_complex, :s1, i, j, 2)
        s2(i, j) = Transition(h_complex, :s2, i, j, 3)

        # Hamiltonian construction
        H =
            Δ_ * a_complex' * a_complex +
            g * (a_complex' + a_complex) * (s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2)) +
            η * (a_complex' + a_complex)

        J = [a_complex]
        rates = [κ]

        # Test double average operations
        imH = im * H
        op_ = a_complex'a_complex
        rhs_ = commutator(imH, op_)
        rhs_avg = average(rhs_)
        rhs_avg_simplified = SymbolicUtils.simplify(rhs_avg)
        terms = SecondQuantizedAlgebra.undo_average(rhs_avg_simplified)

        for arg in arguments(terms)
            @test arg isa SecondQuantizedAlgebra.QMul
        end
    end
end

@testset "QuantumCumulants PR 242" begin
    @testset "spin" begin
        ha1 = NLevelSpace(:atom1, 2)
        ha2 = NLevelSpace(:atom2, 2)
        h = ha1 ⊗ ha2
        s1(i, j) = Transition(h, :s1, i, j, 1)
        s2(i, j) = Transition(h, :s2, i, j, 2)

        expr = average(s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2))
        @test expr isa SymbolicUtils.BasicSymbolic{SecondQuantizedAlgebra.CNumber} broken =
            true
        expr = simplify(average(s1(2, 1) + s1(1, 2)))
        @test expr isa SymbolicUtils.BasicSymbolic{SecondQuantizedAlgebra.CNumber}
    end

    @testset "fock space" begin
        h₁ = FockSpace(:cavity1)
        h₂ = FockSpace(:cavity2)
        h₃ = FockSpace(:cavity3)
        h = h₁ ⊗ h₂ ⊗ h₃

        @qnumbers a1::Destroy(h, 1) a2::Destroy(h, 2) a3::Destroy(h, 3)

        expr = average(a1 + a1') |> simplify
        @test expr isa SymbolicUtils.BasicSymbolic{SecondQuantizedAlgebra.CNumber}

        expr = average(a1 + a2') |> simplify
        @test expr isa SymbolicUtils.BasicSymbolic{SecondQuantizedAlgebra.CNumber}

        expr = average(a1 + a2) |> simplify
        @test expr isa SymbolicUtils.BasicSymbolic{SecondQuantizedAlgebra.CNumber}
    end
end
