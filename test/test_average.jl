using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "average" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf⊗ha

    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @test isequal(average(2a), 2*average(a))

    @test isequal(average(2*(a+a)), 2*(average(a) + average(a)))
    @test isequal(average(a^2), average(a*a))

    @test isequal(simplify(average(σ)+average(σ)), average(2σ))
    @test iszero(simplify(average(σ)-average(σ)))

    ωc, ωa = cnumbers("ω_c ω_a")
    @test isequal(average(ωc), ωc)
    @test isequal(average(ωc*a), ωc*average(a))
    @test isequal(average(ωc*(a+a')), ωc*average(a) + ωc*average(a'))

    @test iszero(average(a*σ)*average(a) - average(a)*average(a*σ))

    @testset "double average" begin # https://github.com/qojulia/QuantumCumulants.jl/issues/242
        @cnumbers Δ_ g κ η
        hf = FockSpace(:cavity)
        ha1 = NLevelSpace(:atom1, 2)
        ha2 = NLevelSpace(:atom2, 2)
        h = hf ⊗ ha1 ⊗ ha2
        a = Destroy(h, :a)
        s1(i, j) = Transition(h, :s1, i, j, 2)
        s2(i, j) = Transition(h, :s2, i, j, 3)
        H = Δ_*a'*a + g*(a' + a)*(s1(2, 1) + s1(1, 2) + s2(2, 1) + s2(1, 2)) + η*(a' + a)
        J = [a]
        rates = [κ]

        imH = im*H
        op_ = a'a
        rhs_ = commutator(imH, op_)
        rhs_avg = average(rhs_)
        rhs_avg_simplified = SymbolicUtils.simplify(rhs_avg)
        terms = SecondQuantizedAlgebra.undo_average(rhs_avg_simplified)
        for arg in arguments(terms)
            @test arg isa SecondQuantizedAlgebra.QMul
        end
    end
end # testset
