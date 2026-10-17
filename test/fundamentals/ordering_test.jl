using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: canonical_isless, QAdd, SuperOperator, spre, spost, channels

@testset "Public ordering contract" begin
    h = FockSpace(:fock)
    hn = NLevelSpace(:atom, 3)
    hp = PauliSpace(:pauli)
    hs = SpinSpace(:spin)
    hq = PhaseSpace(:phase)

    precedes(x, y) = canonical_isless(x, y)
    trichotomous(x, y) = precedes(x, y) + precedes(y, x) + isequal(x, y) == 1

    @testset "identity-faithful order" begin
        operators = [
            Destroy(h, :a),
            Create(h, :a),
            Transition(hn, :σ, 1, 2),
            Transition(hn, :σ, 1, 3),
            Pauli(hp, :σ, 1),
            Spin(hs, :S, 2),
            Position(hq, :x),
            Momentum(hq, :p),
        ]
        for left in operators, right in operators
            @test trichotomous(left, right)
        end
    end

    @testset "the order makes sorting reproducible" begin
        a = Destroy(h, :a)
        b = Destroy(h, :b)
        @test sort([b, a]; lt = canonical_isless) == [a, b]
        @test precedes(a, b) ⊻ precedes(b, a)

        i = Index(h, :i, 5, h)
        j = Index(h, :j, 5, h)
        @test trichotomous(IndexedOperator(Destroy(h, :d), i), IndexedOperator(Destroy(h, :d), j))
        @test !isequal(IndexedOperator(Destroy(h, :d), i), IndexedOperator(Destroy(h, :d), j))

        product = FockSpace(:left) ⊗ FockSpace(:right)
        first_site = Create(product, :x, 1)
        second_site = Destroy(product, :x, 2)
        @test precedes(first_site, second_site)
    end

    @testset "sums have a total order" begin
        a = Destroy(h, :a)
        ad = Create(h, :a)
        qa = a * 1
        qad = ad * 1
        for x in (qa, qad, 2 * a, ad * a, qa + qad), y in (qa, qad, 2 * a, ad * a, qa + qad)
            @test trichotomous(x, y)
        end
        @test !precedes(qa, qa)
        reordered = 2 * a + 3 * ad
        @test !precedes(reordered, 3 * ad + 2 * a) && !precedes(3 * ad + 2 * a, reordered)
        @test precedes(qa, ad * a)
    end

    @testset "maps and channels have a total order" begin
        a = Destroy(h, :a)
        maps = [spre(a), spost(a), dissipator(a), 2 * spre(a)]
        for x in maps, y in maps
            @test trichotomous(SuperOperator(x), SuperOperator(y))
        end
        entries = first.(channels(dissipator(a) + dissipator(a' * a)))
        @test trichotomous(entries[1], entries[2])
    end

    @testset "indices take part in the order" begin
        hi = NLevelSpace(:sites, 2)
        @variables N M
        i = Index(hi, :i, N, hi)
        iM = Index(hi, :i, M, hi)
        σ(k) = IndexedOperator(Transition(hi, :σ, 2, 1), k)
        @test trichotomous(σ(i), σ(iM))
        @test trichotomous(Σ(σ(i), i), σ(i) * 1)
        @test trichotomous(Σ(σ(i), i), Σ(σ(iM), iM))
        @test trichotomous(SuperOperator(Σ(spre(σ(i)), i)), SuperOperator(spre(σ(i))))
        sums = [Σ(σ(i), i), σ(i) * 1, σ(iM) * 1, 2 * σ(i)]
        @test sort(sums; lt = canonical_isless) == sort(reverse(sums); lt = canonical_isless)
    end

    @testset "operators and sums compare with each other" begin
        a = Destroy(h, :a)
        mixed = Union{Op, QAdd}[a' * a, a, a * 1, a']
        @test sort(mixed; lt = canonical_isless) == sort(reverse(mixed); lt = canonical_isless)
        @test canonical_isless(a, a * 1) && !canonical_isless(a * 1, a)
        @test @inferred(canonical_isless(mixed[1], mixed[2])) isa Bool
    end

    @testset "it is not Base.isless" begin
        a = Destroy(h, :a)
        b = Destroy(h, :b)
        @test_throws MethodError isless(a, b)
        @test_throws MethodError a < b
        @test_throws MethodError (a + b) < 2a
        @test_throws MethodError sort([b, a])
    end

    @testset "ordering distinguishes type-specific fields" begin
        a = Destroy(h, :a)
        i = Index(h, :i, 5, h)
        j = Index(h, :j, 5, h)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        @test precedes(ai, aj) ⊻ precedes(aj, ai)
        hn2 = NLevelSpace(:atom2, 2)
        σ12 = Transition(hn2, :σ, 1, 2)
        σ21 = Transition(hn2, :σ, 2, 1)
        @test precedes(σ12, σ21) ⊻ precedes(σ21, σ12)
        hp2 = PauliSpace(:p2)
        px = Pauli(hp2, :σ, 1)
        py = Pauli(hp2, :σ, 2)
        @test precedes(px, py) ⊻ precedes(py, px)
        # strict total order trichotomy on a small set
        ops = [a, ai, aj, σ12, px]
        for o1 in ops, o2 in ops
            @test trichotomous(o1, o2)
        end
    end

    @testset "orders by space_index before kind" begin
        h2 = NLevelSpace(:atom, 2) ⊗ FockSpace(:cavity)
        σ = Transition(h2, :σ, 1, 2, 1)
        a2 = Destroy(h2, :a, 2)
        @test precedes(σ, a2)
        @test !precedes(a2, σ)
        @test sort([a2, σ]; lt = canonical_isless) == [σ, a2]
        @test isequal(σ * a2, a2 * σ)
    end
end
