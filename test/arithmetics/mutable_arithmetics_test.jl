using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test
import MutableArithmetics as MA
import SecondQuantizedAlgebra:
    QAdd, Index, IndexedOperator, Σ, constraint_pairs, ⊗,
    _QAddBuilder, _build

# Basic correctness vs the explicit `+`-chain is already covered by the
# "sum and reduce" testset in test/expressions/algebra_test.jl (those
# assertions exercise this fast path once the override is loaded). Here we add
# the coverage that testset lacks: coefficient cancellation, indexed-`Σ`
# integration, allocation scaling, inference, and the MA interface.
@testset "mutable_arithmetics" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    ad = a'

    @testset "Coefficient / cancellation correctness" begin
        @variables x y
        # mixed-space + symbolic coefficients
        hf2 = ⊗(FockSpace(:f1), FockSpace(:f2))
        b1 = Destroy(hf2, :b1, 1)
        b2 = Destroy(hf2, :b2, 2)
        terms = [x * b1, y * b2, b1' * b1, 2 * b2' * b2]
        @test sum(terms) == foldl(+, terms)
        @test reduce(+, terms) == foldl(+, terms)

        # polynomial cancellation through the `_addto_key!` gate
        pc = [((x + y)^2) * a, -(x^2 + 2x * y + y^2) * a]
        @test iszero(sum(pc))
        @test sum(pc) == foldl(+, pc)
    end

    @testset "Indexed-Σ integration" begin
        hf = FockSpace(:f)
        af = Destroy(hf, :a)
        i = Index(hf, :i, 10, hf)
        j = Index(hf, :j, 10, hf)
        ai = IndexedOperator(af, i)
        aj = IndexedOperator(af, j)

        s1 = Σ(ai, i)         # indexed, no constraint
        s2 = Σ(aj' * ai, i)   # diagonal split → off-diagonal term carries ne (i≠j)
        vec = [s1, s2]
        chain = s1 + s2

        @test sum(vec) == chain
        @test reduce(+, vec) == chain
        @test sum(vec).indices == chain.indices
        @test Set(constraint_pairs(sum(vec))) == Set(constraint_pairs(chain))
        @test Set(constraint_pairs(reduce(+, vec))) == Set(constraint_pairs(chain))
    end

    @testset "Value semantics & shared-const guard" begin
        terms = [a + ad, 2 * a, 3 * ad]
        args_before = copy(terms[1].arguments)
        sum(terms)
        reduce(+, terms)
        @test terms[1].arguments == args_before

        # the empty sum returns the shared `_ZERO_QADD`; later ops must not touch it
        z = sum(QAdd[])
        a + ad
        ad * a
        @test iszero(z)
    end

    @testset "Allocation scaling" begin
        function build_terms(M)
            hs = ⊗([FockSpace(Symbol(:m, k)) for k in 1:M]...)
            return [
                Float64(k) * (Destroy(hs, Symbol(:a, k), k)' * Destroy(hs, Symbol(:a, k), k))
                    for k in 1:M
            ]
        end
        t8 = build_terms(8)
        t16 = build_terms(16)
        sum(t8); sum(t16)  # warmup
        # linear (not quadratic) growth: doubling M must not triple allocations
        @test @allocations(sum(t16)) < 3 * @allocations(sum(t8))
    end

    @testset "Inference" begin
        terms = [a + ad, 2 * a, 3 * ad]
        @test (@inferred sum(terms)) isa QAdd
        @test (@inferred reduce(+, terms)) isa QAdd
        @test (@inferred _build(_QAddBuilder())) isa QAdd
    end

    @testset "MA interface" begin
        b = _QAddBuilder()
        MA.operate!!(+, b, ad * a)
        MA.operate!!(+, b, ad * a)
        @test _build(b) == 2 * (ad * a)

        # @rewrite over a +-chain equals the chain
        t1 = a + ad
        t2 = 2 * a
        t3 = 3 * ad
        @test (MA.@rewrite t1 + t2 + t3) == t1 + t2 + t3
    end
end
