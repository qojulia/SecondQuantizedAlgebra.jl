using SecondQuantizedAlgebra
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using QuantumOpticsBase
using Test
import SecondQuantizedAlgebra: simplify, QAdd, QSym, QTermDict, _to_cnum, Index, sorted_arguments

@testset "integration" begin
    @testset "Basic operator creation and algebra" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        ad = a'

        @test hash(a) != hash(ad)
        @test isequal(a, ad')

        # Lazy multiplication
        @test a * ad isa QAdd
        @test ad * a isa QAdd
    end

    @testset "Commutation via normal_order" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a * a† normal-ordered → a†a + 1
        m = a * a'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 2
    end

    @testset "Distinct modes on same space don't commute" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        a = Destroy(h1, :a)
        b = Destroy(h2, :b)
        # a and b have same space_index=1 but different names
        m = a * b'
        result = normal_order(m)
        @test length(result) == 1  # no commutation applied
    end

    @testset "Numeric conversion round-trip" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)

        @test to_numeric(a, b) == destroy(b)
        @test to_numeric(a', b) == create(b)
        @test to_numeric(a' * a, b) == create(b) * destroy(b)

        α = 0.1 + 0.2im
        ψ = coherentstate(b, α)
        @test numeric_average(a, ψ) ≈ α
        @test numeric_average(a' * a, ψ) ≈ abs2(α)

        # numeric_average on QAdd
        expr = a + a' * a
        @test numeric_average(expr, ψ) ≈ α + abs2(α)

        # to_numeric with scalar QAdd (empty operators)
        @test to_numeric(QAdd(QTermDict(QSym[] => _to_cnum(3)), Index[], Tuple{Index, Index}[]), b) == 3 * one(b)
    end

    @testset "Multi-space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

        # Operators on different spaces — normal_order doesn't commute them
        m = a * b'
        result = normal_order(m)
        @test length(result) == 1

        # Same-space commutation still works
        m2 = a * a'
        result2 = normal_order(m2)
        @test length(result2) == 2
    end

    @testset "Symbolic parameters" begin
        @variables g ω
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        expr = g * a' * a + ω * a' * a
        result = simplify(expr)
        @test length(result) == 1

        # Printing with symbolic prefactors should not crash
        @test repr(g * a) isa String
        @test repr(g * a + ω * a') isa String
    end

    @testset "Division with non-integer" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test (a / 2.0) isa QAdd
        @test prefactor(only(sorted_arguments(a / 2.0))) ≈ 0.5

        @test (a / 2) isa QAdd
        @test prefactor(only(sorted_arguments(a / 2))) == 1 // 2

        @test ((a + a') / 2.0) isa QAdd
    end

    @testset "Power edge cases" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^0 = identity (scalar QAdd)
        m0 = a^0
        @test m0 isa QAdd
        @test isempty(operators(only(sorted_arguments(m0))))
        @test prefactor(only(sorted_arguments(m0))) == 1
    end

    @testset "Adjoint with complex prefactor" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        m = (2.0 + 3.0im) * a' * a
        ma = m'
        @test prefactor(only(sorted_arguments(ma))) ≈ conj(2.0 + 3.0im)
    end

    @testset "Higher-order normal ordering" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        # a^2 * (a')^2 should produce multiple terms
        m = a * a * a' * a'
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) >= 3
    end

    @testset "== consistency" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)

        @test QAdd(QTermDict(QSym[a] => _to_cnum(1)), Index[], Tuple{Index, Index}[]) == QAdd(QTermDict(QSym[a] => _to_cnum(1)), Index[], Tuple{Index, Index}[])
        @test (a + a') == (a + a')
    end

    @testset "Jaynes-Cummings model" begin
        hc = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2, 1)
        h = hc ⊗ ha
        @qnumbers a::Destroy(h, 1)
        σ12 = Transition(h, :σ, 1, 2, 2)
        σ21 = Transition(h, :σ, 2, 1, 2)

        @variables g_jc
        H_int = g_jc * (a' * σ12 + a * σ21)
        @test H_int isa QAdd

        bf = FockBasis(5)
        bn = NLevelBasis(2)
        bc = bf ⊗ bn
        H_num = to_numeric(H_int, bc)
        @test H_num !== nothing
    end

    @testset "Pauli algebra via normal_order" begin
        h = PauliSpace(:p)
        σx = Pauli(h, :σ, 1)
        σy = Pauli(h, :σ, 2)
        σz = Pauli(h, :σ, 3)

        # σx·σy·σz = iσz·σz = i·I
        m = σx * σy * σz
        result = normal_order(m)
        @test result isa QAdd
        @test length(result) == 1
        @test isempty(operators(only(sorted_arguments(result))))
        @test prefactor(only(sorted_arguments(result))) == im
    end

    @testset "Mixed Fock + Pauli" begin
        h = FockSpace(:c) ⊗ PauliSpace(:p)
        @qnumbers a::Destroy(h, 1)
        σz = Pauli(h, :σ, 3, 2)

        m = a' * a * σz
        @test m isa QAdd
        @test length(operators(only(sorted_arguments(m)))) == 3

        # Normal order: commute a,a† on space 1, leave σz on space 2
        result = normal_order(a * a' * σz)
        @test result isa QAdd
    end

    @testset "Superradiant laser commutators (PR #99 review)" begin
        # Regression test: diagonal split during *(QAdd, QAdd) must substitute
        # sum_idx → ext_idx BEFORE site-sort so same-site composition rules
        # (transition: σ_αβ · σ_βγ = σ_αγ) fire on the correct physical order.
        # See https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/
        # Eqs. (6) and (7).
        hc = FockSpace(:cavity)
        ha = NLevelSpace(:atom, 2)
        h = hc ⊗ ha

        @qnumbers a::Destroy(h, 1)
        σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

        @variables N Δ
        g(idx) = IndexedVariable(:g, idx)

        i = Index(h, :i, N, ha)
        j = Index(h, :j, N, ha)
        k = Index(h, :k, N, ha)

        H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

        # Eq. (6): i [H, a'σⱼ¹²] = -iΔ a'σⱼ¹² + Σᵢ≠ⱼ i g(i) σᵢ²¹σⱼ¹²
        #                          + i g(j) (a'a (σⱼ²² - σⱼ¹¹) + σⱼ²²)
        op3 = a' * σ(1, 2, j)
        result3 = 1im * commutator(H, op3)
        # Verify each expected dict entry is present with correct prefactor.
        d3 = result3.arguments
        # -iΔ * a' * σⱼ¹²
        @test any(d3) do (ops, c)
            ops == [a', σ(1, 2, j)] && isequal(c, _to_cnum(-1im * Δ))
        end
        # i * g(j) * σⱼ²²
        @test any(d3) do (ops, c)
            ops == [σ(2, 2, j)] && isequal(c, _to_cnum(1im * g(j)))
        end
        # i * g(j) * a'a * σⱼ²²
        @test any(d3) do (ops, c)
            ops == [a', a, σ(2, 2, j)] && isequal(c, _to_cnum(1im * g(j)))
        end
        # -i * g(j) * a'a * σⱼ¹¹
        @test any(d3) do (ops, c)
            ops == [a', a, σ(1, 1, j)] && isequal(c, _to_cnum(-1im * g(j)))
        end
        # Σᵢ≠ⱼ i * g(i) * σᵢ²¹σⱼ¹²
        @test any(d3) do (ops, c)
            ops == [σ(2, 1, i), σ(1, 2, j)] && isequal(c, _to_cnum(1im * g(i)))
        end
        @test any(p -> p == (i, j) || p == (j, i), result3.non_equal)

        # Eq. (7): i [H, σⱼ¹²σₖ²¹] = i g(j) a (σⱼ²² - σⱼ¹¹) σₖ²¹
        #                            + i g(k) a' σⱼ¹² (σₖ¹¹ - σₖ²²)
        op4 = σ(1, 2, j) * σ(2, 1, k)
        result4 = 1im * commutator(H, op4)
        d4 = result4.arguments
        @test any(d4) do (ops, c)
            ops == [a, σ(2, 2, j), σ(2, 1, k)] && isequal(c, _to_cnum(1im * g(j)))
        end
        @test any(d4) do (ops, c)
            ops == [a, σ(1, 1, j), σ(2, 1, k)] && isequal(c, _to_cnum(-1im * g(j)))
        end
        @test any(d4) do (ops, c)
            ops == [a', σ(1, 2, j), σ(1, 1, k)] && isequal(c, _to_cnum(1im * g(k)))
        end
        @test any(d4) do (ops, c)
            ops == [a', σ(1, 2, j), σ(2, 2, k)] && isequal(c, _to_cnum(-1im * g(k)))
        end
    end
end
