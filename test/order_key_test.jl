using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: order_key, qadd_order_key, term_order_key

@testset "Public ordering contract" begin
    h = FockSpace(:fock)
    hn = NLevelSpace(:atom, 3)
    hp = PauliSpace(:pauli)
    hs = SpinSpace(:spin)
    hq = PhaseSpace(:phase)

    @testset "identity-faithful keys" begin
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
            @test (order_key(left) == order_key(right)) == isequal(left, right)
        end
    end

    @testset "keys make public sorting reproducible" begin
        a = Destroy(h, :a)
        b = Destroy(h, :b)
        @test order_key(a) != order_key(b)
        @test sort([b, a]) == [a, b]
        @test isless(order_key(a), order_key(b)) ⊻ isless(order_key(b), order_key(a))

        i = Index(h, :i, 5, h)
        j = Index(h, :j, 5, h)
        @test order_key(IndexedOperator(Destroy(h, :d), i)) !=
            order_key(IndexedOperator(Destroy(h, :d), j))

        product = FockSpace(:left) ⊗ FockSpace(:right)
        first_site = Create(product, :x, 1)
        second_site = Destroy(product, :x, 2)
        @test isless(order_key(first_site), order_key(second_site))
    end

    @testset "products and sums have total public keys" begin
        a = Destroy(h, :a)
        ad = Create(h, :a)
        term_a = first(collect(a * 1)).first
        term_ad = first(collect(ad * 1)).first
        @test (term_order_key(term_a) == term_order_key(term_ad)) == isequal(term_a, term_ad)

        qa = a * 1
        qad = ad * 1
        @test qadd_order_key(qa) != qadd_order_key(qad)
        @test qadd_order_key(qa) != qadd_order_key(2 * a)
        @test isless(qa, qad) == isless(qadd_order_key(qa), qadd_order_key(qad))
        @test isless(qa, 2 * qa) == isless(qadd_order_key(qa), qadd_order_key(2 * qa))
        @test !isless(qa, qa)

        reordered = 2 * a + 3 * ad
        @test isequal(reordered, 3 * ad + 2 * a)
        @test qadd_order_key(reordered) == qadd_order_key(3 * ad + 2 * a)

        longer = first(collect(ad * a)).first
        @test term_order_key(term_a) < term_order_key(longer)
    end

    @testset "public key is inferable" begin
        a = Destroy(h, :a)
        @inferred order_key(a)
    end

    @testset "ordering distinguishes type-specific fields" begin
        a = Destroy(h, :a)
        i = Index(h, :i, 5, h)
        j = Index(h, :j, 5, h)
        ai = IndexedOperator(a, i)
        aj = IndexedOperator(a, j)
        @test order_key(ai) != order_key(aj)
        hn2 = NLevelSpace(:atom2, 2)
        σ12 = Transition(hn2, :σ, 1, 2)
        σ21 = Transition(hn2, :σ, 2, 1)
        @test order_key(σ12) != order_key(σ21)
        hp2 = PauliSpace(:p2)
        px = Pauli(hp2, :σ, 1)
        py = Pauli(hp2, :σ, 2)
        @test order_key(px) != order_key(py)
        # strict total order trichotomy on a small set
        ops = [a, ai, aj, σ12, px]
        for o1 in ops, o2 in ops
            k1, k2 = order_key(o1), order_key(o2)
            @test (k1 < k2) + (k2 < k1) + (k1 == k2) == 1
        end
    end

    @testset "orders by space_index before kind" begin
        h2 = NLevelSpace(:atom, 2) ⊗ FockSpace(:cavity)
        σ = Transition(h2, :σ, 1, 2, 1)
        a2 = Destroy(h2, :a, 2)
        @test order_key(σ) < order_key(a2)
        @test isless(σ, a2)
        @test !isless(a2, σ)
        @test sort([a2, σ]) == [σ, a2]
        @test isequal(σ * a2, a2 * σ)
    end
end
