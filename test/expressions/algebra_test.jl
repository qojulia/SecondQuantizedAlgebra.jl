using SecondQuantizedAlgebra
using Symbolics: @variables
using Test
import SecondQuantizedAlgebra: QAdd, QSym, QField, sorted_arguments, CNum,
    _CNUM_ONE, _single_qadd, _to_cnum, simplify

@testset "Algebra (QAdd / QMul construction)" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Basic operator creation" begin
        @qnumbers b::Destroy(h)
        @test hash(b) != hash(b')
        @test isequal(b, b'')
        @test b * b' isa QAdd
        @test b' * b isa QAdd
    end

    @testset "QAdd construction" begin
        @testset "QAdd + QAdd (single-term)" begin
            s = (2 * a) + (3 * ad)
            @test s isa QAdd
            @test length(s) == 2
        end

        @testset "QAdd + QSym" begin
            s1 = (2 * a) + ad
            @test s1 isa QAdd
            @test length(s1) == 2

            s2 = ad + (2 * a)
            @test s2 isa QAdd
        end

        @testset "QAdd + QAdd (merge)" begin
            s = (a + ad) + (3 * a)
            @test s isa QAdd
            @test length(s) == 2  # a and 3a auto-collect to 4a, plus a†
        end

        @testset "QAdd + QAdd" begin
            s = (a + ad) + (a + ad)
            @test s isa QAdd
            @test length(s) == 2  # auto-collects to 2a + 2a†
        end

        @testset "QField + Number" begin
            s1 = a + 5
            @test s1 isa QAdd
            @test length(s1) == 2

            @test 5 + a isa QAdd
            s3 = (a + ad) + 3
            @test s3 isa QAdd
            @test length(s3) == 3
            @test 3 + (a + ad) isa QAdd
        end

        @testset "Subtraction" begin
            s = a - ad
            @test s isa QAdd
            @test length(s) == 2
            @test (a - 3) isa QAdd
        end

        @testset "QAdd * Number" begin
            s = (a + ad) * 3
            @test s isa QAdd
            @test all(t -> prefactor(t) == 3, sorted_arguments(s))
        end

        @testset "QAdd * QSym (distributes)" begin
            s = (a + ad) * a
            @test s isa QAdd
            @test length(s) == 2
            @test all(t -> length(operators(t)) == 2, sorted_arguments(s))
        end

        @testset "QAdd * QAdd (distributes, single-term rhs)" begin
            s = (a + ad) * (2 * a)
            @test s isa QAdd
            @test length(s) == 2
        end

        @testset "QAdd * QAdd (full distribute)" begin
            s = (a + ad) * (a + ad)
            @test s isa QAdd
            @test length(s) == 4
        end

        @testset "Adjoint of QAdd" begin
            s = a + 2 * ad
            sd = s'
            @test sd isa QAdd
            @test length(sd) == 2
        end
    end

    @testset "Multiplication (eager ordering)" begin
        @testset "Eager normal ordering" begin
            # a * a† = a†a + 1 (normal ordered)
            result = a * ad
            @test result[QSym[ad, a]] == _CNUM_ONE
            @test result[QSym[]] == _CNUM_ONE
            @test length(result) == 2

            # a† * a is already normal ordered
            result2 = ad * a
            @test result2[QSym[ad, a]] == _CNUM_ONE
            @test length(result2) == 1
        end

        @testset "Cross-site operators" begin
            h2 = FockSpace(:a) ⊗ FockSpace(:b)
            a1 = Destroy(h2, :a, 1)
            a2 = Destroy(h2, :b, 2)
            m = a1 * a2
            @test m isa QAdd
            @test length(m) == 1
            ops = first(keys(m.arguments)).ops
            @test ops[1].space_index <= ops[2].space_index
        end

        @testset "Division" begin
            m = a / 2
            @test m isa QAdd
            @test m[QSym[a]] == 1 // 2
            @test (a / 2.0) isa QAdd
            @test prefactor(only(sorted_arguments(a / 2.0))) ≈ 0.5
            @test prefactor(only(sorted_arguments(a / 2))) == 1 // 2
            @test ((a + a') / 2.0) isa QAdd
        end

        @testset "Power" begin
            m = a^3
            @test m isa QAdd
            @test length(first(keys(m.arguments)).ops) == 3

            # a^0 = identity (scalar QAdd)
            m0 = a^0
            @test m0 isa QAdd
            @test length(m0) == 1
            @test m0[QSym[]] == _CNUM_ONE
            @test isempty(operators(only(sorted_arguments(m0))))
            @test prefactor(only(sorted_arguments(m0))) == 1
        end

        @testset "Negation" begin
            m = -a
            @test m isa QAdd
            @test m[QSym[a]] == -1
        end

        @testset "Adjoint of QMul" begin
            m = 2 * ad * a
            @test isequal(m', m)  # (2 a†a)† = 2 a†a
        end

        @testset "Adjoint with complex prefactor" begin
            m = (2.0 + 3.0im) * ad * a
            ma = m'
            @test prefactor(only(sorted_arguments(ma))) ≈ conj(2.0 + 3.0im)
        end
    end

    @testset "Symbolic parameters" begin
        @variables g ω
        expr = g * ad * a + ω * ad * a
        @test length(simplify(expr)) == 1
        @test repr(g * a) isa String
        @test repr(g * a + ω * ad) isa String
    end

    @testset "== consistency" begin
        @test _single_qadd(_to_cnum(1), QSym[a]) == _single_qadd(_to_cnum(1), QSym[a])
        @test (a + a') == (a + a')
    end

    @testset "Type stability" begin
        # QSym + QSym / QAdd + ...
        @inferred a + ad
        @inferred (2 * a) + (3 * ad)
        @inferred (2 * a) + ad
        @inferred ad + (2 * a)
        @inferred (a + ad) + (3 * a)
        @inferred (3 * a) + (a + ad)
        @inferred (a + ad) + a
        @inferred a + (a + ad)
        @inferred (a + ad) + (a + ad)
        # QField + Number / Number + QField
        @inferred a + 5
        @inferred 5 + a
        @inferred (2 * a) + 5
        @inferred 5 + (2 * a)
        @inferred (a + ad) + 3
        @inferred 3 + (a + ad)
        # Subtraction
        @inferred a - ad
        @inferred a - 3
        @inferred 3 - a
        # Negation
        @inferred -(a + ad)
        @inferred -a
        # Multiplication
        @inferred a * ad
        @inferred 3 * a
        @inferred a * 2.0
        @inferred (a + ad) * 3
        @inferred 3 * (a + ad)
        @inferred (a + ad) * a
        @inferred a * (a + ad)
        @inferred (a + ad) * (2 * a)
        @inferred (2 * a) * (a + ad)
        @inferred (a + ad) * (a + ad)
        # Division and power
        @inferred a / 2
        @inferred (a + ad) / 2
        @inferred (a + ad) / 2.0
        @inferred a^3
        # Adjoint
        @inferred adjoint(a + 2 * ad)
        # Equality / hashing
        @inferred isequal(a + ad, a + ad)
        @inferred hash(a + ad, UInt(0))
    end

    @testset "sum and reduce" begin
        terms = [a + ad, 2 * a, 3 * ad, a * ad, ad * a]
        manual = terms[1] + terms[2] + terms[3] + terms[4] + terms[5]

        @test sum(terms) == manual
        @test reduce(+, terms) == manual

        args_before = copy(terms[1].arguments)
        sum(terms)
        @test terms[1].arguments == args_before

        @test sum([a + ad]) == a + ad
        @test sum(QAdd[]) == zero(QAdd)
        @test sum(x -> x * a, [a, ad]) == a * a + ad * a

        @test iszero(zero(QAdd))
        @test zero(a + ad) == zero(QAdd)
    end

    @testset "Allocations" begin
        s = a + ad
        m_2a = 2 * a
        m_3ad = 3 * ad

        # Warmup
        a + ad; m_2a + m_3ad; s + m_2a; s + s; s * 3; s * a; s * s

        @test @allocations(a + ad) <= 200
        @test @allocations(m_2a + m_3ad) <= 200
        @test @allocations(s + m_2a) <= 200
        @test @allocations(s + s) <= 200
        @test @allocations(s * 3) <= 500
        @test @allocations(s * a) <= 200
        @test @allocations(s * s) <= 10000
        @test @allocations(-s) <= 500
        @test @allocations(adjoint(s)) <= 200
    end
end
