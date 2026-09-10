using SecondQuantizedAlgebra
using Test

import SecondQuantizedAlgebra: QAdd

@testset "Formal operator Taylor lowering" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    A = a + a'

    c4 = @inferred taylor(cos(A), 0:4)
    s5 = @inferred taylor(sin(A), 0:5)
    e4 = @inferred taylor(expim(A), 0:4)

    @test c4 isa QAdd
    @test isequal(c4, 1 - (1 // 2) * A^2 + (1 // 24) * A^4)
    @test isequal(s5, A - (1 // 6) * A^3 + (1 // 120) * A^5)
    @test isequal(
        e4,
        1 + im * A - (1 // 2) * A^2 - im * (1 // 6) * A^3 + (1 // 24) * A^4,
    )

    @test isequal(taylor(cos(A), 0:0), 1)
    @test isequal(taylor(cos(A), 0:1), 1)
    @test iszero(taylor(sin(A), 0:0))
    @test isequal(taylor(sin(A), 0:1), A)
    @test isequal(taylor(expim(A), 0:0), 1)
    @test isequal(taylor(expim(A), 0:1), 1 + im * A)

    @test isequal(taylor(a * cos(A), 0:2), a * (1 - (1 // 2) * A^2))
    @test isequal(taylor(cos(A) * a, 0:2), (1 - (1 // 2) * A^2) * a)
    @test !isequal(taylor(a * cos(A), 0:2), taylor(cos(A) * a, 0:2))

    wrapped = cos(A) - cos(A) + A
    @test wrapped isa SecondQuantizedAlgebra.QExpr
    @test isequal(taylor(cos(wrapped), 0:2), 1 - (1 // 2) * A^2)

    local_order = 1 - (1 // 2) * A^2
    @test isequal(taylor(cos(A) * cos(A), 0:2), local_order * local_order)

    @test_throws ArgumentError taylor(cos(A), 1:4)
    @test_throws MethodError taylor(cos(A), 0:2:4)
    @test_throws ArgumentError taylor(cos(sin(A)), 0:4)

    i = Index(h, :i, 3, h)
    ai = IndexedOperator(a, i)
    summed = Σ(ai + ai', i)
    @test_throws ArgumentError taylor(cos(summed), 0:2)

    atom = NLevelSpace(:atom, 2)
    σ11 = Transition(atom, :σ, 1, 1)
    high_order = taylor(cos(σ11), 0:22)
    coeff = sum((-1)^k * (big(1) // factorial(big(2k))) for k in 1:11)
    @test isequal(high_order, 1 + coeff * σ11)
end
