using Test
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics

const sqa = SecondQuantizedAlgebra

@testset "average_sums" begin
    # Common setup for all tests
    N = 2
    ha = NLevelSpace(Symbol(:atom), 2)
    hf = FockSpace(:cavity)
    h = hf ⊗ ha

    ind(i) = Index(h, i, N, ha)

    g(k) = IndexedVariable(:g, k)
    Γij = DoubleIndexedVariable(:Γ, ind(:i), ind(:j))
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    σn(i, j, k) = NumberedOperator(Transition(h, :σ, i, j), k)
    Ω(i, j) = IndexedVariable(:Ω, i, j; identical=false)

    a = Destroy(h, :a)

    @testset "Basic Average Operations" begin
        @test Ω(ind(:i), ind(:i)) == 0
        @test isequal(average(2 * σ(1, 2, ind(:k))), 2 * average(σ(1, 2, ind(:k))))
        @test isequal(
            average(g(ind(:k)) * σ(2, 2, ind(:k))), g(ind(:k)) * average(σ(2, 2, ind(:k)))
        )
        @test isequal(average(g(ind(:k))), g(ind(:k)))
    end

    @testset "Numbered Operators" begin
        @test isequal(
            σn(1, 2, 1) + σn(2, 1, 1),
            NumberedOperator(Transition(h, :σ, 1, 2) + Transition(h, :σ, 2, 1), 1),
        )
    end

    @testset "Index Insertion" begin
        # Basic index insertion tests
        @test isequal(σn(2, 2, 1), insert_index(σ(2, 2, ind(:j)), ind(:j), 1))
        @test isequal(σ(1, 2, ind(:j)), insert_index(σ(1, 2, ind(:j)), ind(:k), 2))
        @test isequal(1, insert_index(1, ind(:k), 1))

        # Variable index insertion
        gamma = insert_index(Γij, ind(:i), 1)
        gamma2 = insert_index(Γij, ind(:j), 2)
        gamma2_ = insert_index(Γij, ind(:i), 1)
        g_ = insert_index(g(ind(:j)), ind(:j), 1)

        @test g_ isa SymbolicUtils.BasicSymbolic
        @test sqa.is_symtype(gamma, sqa.DoubleNumberedVariable)

        gamma_ = insert_index(gamma, ind(:j), 2)
        @test sqa.is_symtype(gamma_, Complex{Real})

        @test !isequal(gamma, gamma_)
        @test !isequal(gamma, g_)
        @test isequal(gamma2_, gamma)
        @test !isequal(gamma, gamma2)

        # Complex expression index insertion
        @test isequal(
            sqa.insert_index(σ(1, 2, ind(:j)) * σn(1, 2, 2), ind(:j), 1),
            sqa.insert_index(
                sqa.insert_index(σ(1, 2, ind(:i)) * σ(1, 2, ind(:j)), ind(:i), 2),
                ind(:j),
                1,
            ),
        )
    end

    @testset "Index Utilities" begin
        sumterm = σ(1, 2, ind(:i)) * σ(2, 1, ind(:j)) * σ(2, 2, ind(:k))
        inds = sqa.get_indices(sumterm)
        @test isequal([ind(:i), ind(:j), ind(:k)], inds)

        @test ind(:i) ∈ sqa.get_indices(g(ind(:i)))

        @cnumbers N_
        ind2(i) = Index(h, i, N_, ha)
    end

    @testset "Sums and Indexed Average Sums" begin
        sum1 = SingleSum(σ(1, 2, ind(:k)), ind(:k))
        sum2 = average(sum1 * σ(1, 2, ind(:l)))

        @test !isequal(σn(2, 2, 1), insert_index(sum2, ind(:j), 1))

        pind = Index(h, :p, 5, ha)
        @test isequal(4 * a, Σ(a, pind, [ind(:i)]))

        @test sqa.indexed_isequal(
            average(Σ(σ(1, 2, ind(:i)), ind(:i))),
            sqa.IndexedAverageSum(average(σ(1, 2, ind(:i))), ind(:i), []),
        )

        @test sqa.indexed_isequal(
            sqa.IndexedAverageSum(g(ind(:i)), ind(:i), []),
            average(Σ(g(ind(:i)), ind(:i), [])),
        )
        @test sqa.is_symtype(
            sqa.IndexedAverageSum(g(ind(:i)), ind(:i), []), sqa.IndexedAverageSum
        )

        @test sqa.IndexedAverageSum(1) == 1
    end

    @testset "Double Sums and Average Operations" begin
        avrgTerm = average(Σ(σ(2, 1, ind(:i)) * σ(1, 2, ind(:j)), ind(:i)))
        @test avrgTerm isa SymbolicUtils.BasicSymbolic && operation(avrgTerm) === +

        ADsum1 = simplify(sqa.IndexedAverageDoubleSum(avrgTerm, ind(:j), [ind(:i)]))
        it_1 = findfirst(x -> sqa.is_symtype(x, IndexedAverageSum), arguments(ADsum1))
        it_2 = findfirst(x -> sqa.is_symtype(x, IndexedAverageDoubleSum), arguments(ADsum1))

        @test ADsum1 isa SymbolicUtils.BasicSymbolic && operation(ADsum1) === +
        @test SymbolicUtils.metadata(arguments(ADsum1)[it_2])[sqa.IndexedAverageDoubleSum] isa
            sqa.IndexedAverageDoubleSum
        @test SymbolicUtils.metadata(arguments(ADsum1)[it_1])[sqa.IndexedAverageSum] isa
            sqa.IndexedAverageSum

        @test sqa.undo_average(arguments(ADsum1)[it_2]) isa sqa.DoubleSum
        @test isequal(
            simplify(
                Σ(Σ(σ(2, 1, ind(:i)) * σ(1, 2, ind(:j)), ind(:i)), ind(:j), [ind(:i)])
            ),
            simplify(sqa.undo_average(ADsum1)),
        )

        # Argument order tests (SymbolicUtils v1.4.0 compatibility)
        @test (
            isequal(
                sqa.sqa_arguments(sqa.sqa_arguments(ADsum1)[it_2]),
                SymbolicUtils.arguments(avrgTerm)[1],
            ) || isequal(
                sqa.sqa_arguments(sqa.sqa_arguments(ADsum1)[it_2]),
                SymbolicUtils.arguments(avrgTerm)[2],
            )
        )

        @test isequal(
            sqa.sqa_arguments(sqa.sqa_arguments(sqa.sqa_arguments(ADsum1)[it_2])),
            SymbolicUtils.arguments(average(σ(2, 1, ind(:i)) * σ(1, 2, ind(:j)))),
        )
    end

    @testset "Special Indexed Averages" begin
        @test sqa.indexed_isequal(
            sqa.SpecialIndexedAverage(average(σ(1, 2, ind(:i))), [(ind(:i), ind(:j))]) +
            sqa.SpecialIndexedAverage(average(σ(2, 1, ind(:j))), [(ind(:i), ind(:j))]),
            sqa.SpecialIndexedAverage(
                average(σ(1, 2, ind(:i))) + average(σ(2, 1, ind(:j))), [(ind(:i), ind(:j))]
            ),
        )

        @test sqa.SpecialIndexedAverage(average(0), [(ind(:i), ind(:j))]) == 0

        specAvrg = sqa.SpecialIndexedAverage(
            average(σ(2, 1, ind(:i)) * σ(1, 2, ind(:j))), [(ind(:i), ind(:j))]
        )

        @test SymbolicUtils.metadata(
            sqa.SpecialIndexedAverage(average(σ(2, 1, ind(:i))), [(ind(:i), ind(:j))])
        )[sqa.SpecialIndexedAverage] isa sqa.SpecialIndexedAverage

        @test isequal(
            sqa.sqa_arguments(specAvrg),
            SymbolicUtils.arguments(average(σ(2, 1, ind(:i)) * σ(1, 2, ind(:j)))),
        )
    end

    @testset "QMul Operations" begin
        @test σ(1, 2, ind(:i)) * σ(2, 1, ind(:j)) * σn(2, 2, 3) isa sqa.QMul
        @test σn(2, 2, 3) * σ(1, 2, ind(:i)) * σ(2, 1, ind(:j)) isa sqa.QMul
    end

    @testset "Utility Functions" begin
        @test isequal("(i≠1)", sqa.SecondQuantizedAlgebra.writeNeqs([(ind(:i), 1)]))
    end

    @testset "Indexed Operator Detection" begin
        # Helper functions for checking if indices occurred in specific terms
        function containsIndexedOps(term::Average)
            arg_ = map(unwrap_const, arguments(term))
            if arg_[1] isa sqa.QMul
                for arg in arg_[1].args_nc
                    if arg isa sqa.IndexedOperator
                        return true
                    end
                end
            else
                return arg_[1] isa sqa.IndexedOperator
            end
            return false
        end
        containsIndex(term::Average, ind::Index) = ind ∈ get_indices(term)

        @test containsIndexedOps(average(a * σ(2, 1, ind(:i)) * σ(1, 2, ind(:j))))
        @test !(containsIndexedOps(average(a' * a)))
    end
end
