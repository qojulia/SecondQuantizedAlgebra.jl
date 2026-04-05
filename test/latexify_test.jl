using SecondQuantizedAlgebra
using Latexify
using LaTeXStrings
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: simplify, QMul, QAdd, QSym

@testset "latexify" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Operators" begin
        hn = NLevelSpace(:atom, 3, 1)
        hp = PauliSpace(:p)
        hs = SpinSpace(:s, 1 // 2)
        hps = PhaseSpace(:q)

        input = [
            Destroy(h, :a),
            Create(h, :a),
            Transition(hn, :σ, 1, 2),
            Pauli(hp, :σ, 1),
            Pauli(hp, :σ, 2),
            Pauli(hp, :σ, 3),
            Spin(hs, :S, 1),
            Spin(hs, :S, 3),
            Position(hps, :x),
            Momentum(hps, :p),
        ]
        output = [
            L"a",
            L"a^{\dagger}",
            L"{\sigma}^{{12}}",
            L"{\sigma}_{{x}}",
            L"{\sigma}_{{y}}",
            L"{\sigma}_{{z}}",
            L"{S}_{{x}}",
            L"{S}_{{z}}",
            L"\hat{x}",
            L"\hat{p}",
        ]
        for (i, o) in zip(input, output)
            @test latexify(i) == o
        end
    end

    @testset "Transition superscript toggle" begin
        hn = NLevelSpace(:atom, 3, 1)
        σ12 = Transition(hn, :σ, 1, 2)

        transition_superscript(true)
        @test latexify(σ12) == L"{\sigma}^{{12}}"

        transition_superscript(false)
        @test latexify(σ12) == L"{\sigma}_{{12}}"

        transition_superscript(true)  # reset
    end

    @testset "QMul" begin
        input = [
            1 * ad * a,
            3 * ad * a,
            -1 * a,
            QMul(5, QSym[]),
        ]
        output = [
            L"a^{\dagger}a",
            L"3a^{\dagger}a",
            L"-a",
            L"5",
        ]
        for (i, o) in zip(input, output)
            @test latexify(i) == o
        end
    end

    @testset "QAdd" begin
        @test latexify(a + ad) == L"a + a^{\dagger}"
    end

    @testset "Simplify display" begin
        @test latexify(simplify(a * ad)) == L"1 + a^{\dagger}a"
    end

    @testset "Symbolic prefactors" begin
        @variables g
        @test latexify(g * a) == L"ga"
    end

    @testset "Complex prefactors" begin
        result = simplify(Pauli(PauliSpace(:p), :σ, 1) * Pauli(PauliSpace(:p), :σ, 2))
        @test latexify(result) == L"\mathit{i}{\sigma}_{{z}}"
    end

    @testset "Indexed operators" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))

        input = [
            IndexedOperator(Transition(h2, :σ, 1, 2, 2), i),
            IndexedOperator(Destroy(h2, :a, 1), i),
            IndexedOperator(Destroy(h2, :a, 1), i)',
        ]
        output = [
            L"{\sigma}_{i}^{{12}}",
            L"a_{i}",
            L"a_{i}^{\dagger}",
        ]
        for (inp, o) in zip(input, output)
            @test latexify(inp) == o
        end
    end

    @testset "Indexed product" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        @test latexify(IndexedVariable(:g, i) * b' * σ_i) ==
            L"g\left( i \right)b^{\dagger}{\sigma}_{i}^{{12}}"
    end

    @testset "Sum" begin
        @variables N
        h2 = FockSpace(:c) ⊗ NLevelSpace(:atom, 2, 1)
        @qnumbers b::Destroy(h2, 1)
        i = Index(h2, :i, N, NLevelSpace(:atom, 2, 1))
        σ_i = IndexedOperator(Transition(h2, :σ, 1, 2, 2), i)

        H = Σ(IndexedVariable(:g, i) * b' * σ_i, i)
        @test latexify(H) ==
            L"\underset{i}{\overset{N}{\sum}}\left(  g\left( i \right)b^{\dagger}{\sigma}_{i}^{{12}} \right)"
    end

    @testset "MIME text/latex" begin
        @test repr(MIME"text/latex"(), a) == latexify(a)
        @test repr(MIME"text/latex"(), ad) == latexify(ad)
        @test repr(MIME"text/latex"(), a * ad) == latexify(a * ad)
        @test repr(MIME"text/latex"(), a + ad) == latexify(a + ad)
    end
end
