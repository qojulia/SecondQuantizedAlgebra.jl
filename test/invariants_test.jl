using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim, get_sum_indices, get_sum_non_equal,
    has_sum_metadata

@testset "Public algebraic properties" begin
    @testset "canonicalization and average round trips" begin
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2, 1)
        a = Destroy(h, :a, 1)
        σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
        @variables N Δ g
        i = Index(h, :i, N, 2)
        j = Index(h, :j, N, 2)
        k = Index(h, :k, N, 2)

        corpus = [
            a * a',
            a' * a + a * a',
            Σ(g * σ(1, 2, i), i),
            Σ(g * σ(1, 2, i), i) * σ(2, 1, j),
            commutator(-Δ * a' * a + Σ(g * a' * σ(1, 2, i), i), σ(1, 2, j)),
            commutator(
                -Δ * a' * a + Σ(g * a' * σ(1, 2, i), i),
                σ(1, 2, j) * σ(2, 1, k)
            ),
            expand_completeness(Σ(σ(1, 1, i), i)),
            substitute(Σ(g * σ(1, 2, i), i), Dict(g => 2)),
            change_index(Σ(g * σ(1, 2, i), i), i, j),
            Rotation(a, Δ),
        ]

        for expression in corpus[1:(end - 1)]
            @test iszero(simplify(expression - normal_order(expression)))
            @test iszero(simplify(undo_average(average(expression)) - expression))
        end
        @test iszero(
            simplify(
                conjugate(a, corpus[end]) - expim(-Δ) * a,
            )
        )
    end

    @testset "cancellation never leaves observable scope" begin
        h = FockSpace(:sites)
        a = Destroy(h, :a)
        @variables N
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        left = Σ(IndexedOperator(a, i), i)
        right = Σ(IndexedOperator(a, j), j)

        @test iszero(left - left)
        @test isempty(get_indices(left - left))
        @test Set(get_indices(left - right)) == Set([i, j])
        @test get_indices((left + right) - left) == [j]
        @test isempty(get_indices(commutator(left, Σ(IndexedOperator(a', j), j))))
    end

    @testset "public algebraic laws" begin
        h = FockSpace(:cavity)
        a = Destroy(h, :a)
        b = Destroy(h, :b)
        expressions = [a, a', a + a', a' * a, (a + a') * (a - a'), a * b']

        for expression in expressions
            @test iszero(commutator(expression, expression))
            @test isequal(adjoint(adjoint(expression)), expression)
        end
        for (A, B, C) in ((a, a', b), (a + a', a, b'), (a' * a, a + b, a' - b'))
            @test iszero(
                simplify(
                    commutator(A, B + C) -
                        commutator(A, B) - commutator(A, C)
                )
            )
        end

        hs = SpinSpace(:spin)
        Sx = Spin(hs, :S, 1)
        Sy = Spin(hs, :S, 2)
        Sz = Spin(hs, :S, 3)
        jacobi = commutator(commutator(Sx, Sy), Sz) +
            commutator(commutator(Sy, Sz), Sx) +
            commutator(commutator(Sz, Sx), Sy)
        @test iszero(simplify(jacobi))
    end

    @testset "indexed average metadata describes its scope" begin
        h = NLevelSpace(:atom, 2)
        @variables N
        i = Index(h, :i, N, h)
        j = Index(h, :j, N, h)
        σ(k) = IndexedOperator(Transition(h, :σ, 1, 2), k)

        averaged = average(Σ(σ(i), i, [j]))
        @test is_indexed_sum(averaged)
        @test has_sum_metadata(averaged)
        @test get_sum_indices(averaged) == [i]
        @test get_sum_non_equal(averaged) == [(i, j)]
        @test iszero(undo_average(averaged) - Σ(σ(i), i, [j]))

        plain = average(Σ(σ(i), i))
        @test get_sum_indices(plain) == [i]
        @test isempty(get_sum_non_equal(plain))
    end
end
