using SecondQuantizedAlgebra
using Symbolics: @variables
using LinearAlgebra: det, tr, dot, norm, rmul!, lmul!,
    Symmetric, Hermitian, issymmetric, ishermitian
using Test
import SecondQuantizedAlgebra: Coeff, to_cnum, to_complex, expim

@testset "Coefficient linear algebra" begin
    @testset "determinant expands by minors" begin
        @variables x y z θ

        @test isequal(det(reshape(Coeff[to_cnum(x)], 1, 1)), to_cnum(x))
        @test isequal(det(Coeff[to_cnum(x) to_cnum(y); to_cnum(0) to_cnum(2)]), to_cnum(2x))
        @test isequal(det(Coeff[to_cnum(x) to_cnum(1); to_cnum(2) to_cnum(y)]), to_cnum(x * y - 2))
        @test isequal(
            det(
                Coeff[
                    to_cnum(x) to_cnum(1) to_cnum(0)
                    to_cnum(0) to_cnum(y) to_cnum(2)
                    to_cnum(3) to_cnum(0) to_cnum(z)
                ],
            ),
            to_cnum(x * y * z + 6),
        )
        @test iszero(det(Coeff[to_cnum(x) to_cnum(y); to_cnum(x) to_cnum(y)]))
        @test_throws DimensionMismatch det(Coeff[to_cnum(x) to_cnum(y)])

        # The expansion must agree with the pivoted LU that plain numbers take.
        numeric = ComplexF64[1 2 3 4; 5 1 6 7; 8 9 1 2; 3 4 5 1]
        symbolic = Coeff[to_cnum(numeric[i, j]) for i in axes(numeric, 1), j in axes(numeric, 2)]
        @test to_complex(det(symbolic)) ≈ det(numeric)

        # Four rows exercise the recursion past the hand-written 1x1 and 2x2 bases.
        @test isequal(
            det(
                Coeff[
                    to_cnum(x) to_cnum(1) to_cnum(0) to_cnum(0)
                    to_cnum(0) to_cnum(y) to_cnum(1) to_cnum(0)
                    to_cnum(0) to_cnum(0) to_cnum(x) to_cnum(1)
                    to_cnum(1) to_cnum(0) to_cnum(0) to_cnum(y)
                ],
            ),
            to_cnum(x^2 * y^2 - 1),
        )

        # Lower triangular takes the diagonal-product short circuit.
        @test isequal(
            det(
                Coeff[
                    to_cnum(x) to_cnum(0) to_cnum(0)
                    to_cnum(1) to_cnum(y) to_cnum(0)
                    to_cnum(2) to_cnum(3) to_cnum(z)
                ],
            ),
            to_cnum(x * y * z),
        )

        # Exact phases cancel through the expansion instead of expanding to trigonometry.
        @test isone(det(Coeff[expim(θ) to_cnum(0); to_cnum(0) expim(-θ)]))

        @test isequal(tr(Coeff[to_cnum(x) to_cnum(1); to_cnum(2) to_cnum(y)]), to_cnum(x + y))
    end

    @testset "adjoint and transpose of a coefficient" begin
        @variables x y θ

        @test isequal(transpose(to_cnum(x)), to_cnum(x))
        @test isequal(adjoint(to_cnum(x)), conj(to_cnum(x)))
        @test isequal(to_cnum(2im)', to_cnum(-2im))
        @test isequal(expim(θ)', expim(-θ))

        A = Coeff[to_cnum(x) to_cnum(im); to_cnum(0) to_cnum(y)]
        @test isequal(transpose(A), Coeff[to_cnum(x) to_cnum(0); to_cnum(im) to_cnum(y)])
        @test isequal(A', Coeff[to_cnum(x) to_cnum(0); to_cnum(-im) to_cnum(y)])
        @test isequal((A')', A)
        @test issymmetric(Coeff[to_cnum(x) to_cnum(1); to_cnum(1) to_cnum(y)])
        @test ishermitian(Coeff[to_cnum(x) to_cnum(im); to_cnum(-im) to_cnum(y)])

        # The canonicality condition the `Bogoliubov` docstring states, written as written.
        @variables u v
        S = Coeff[to_cnum(u) to_cnum(v); to_cnum(v) to_cnum(u)]
        J = Coeff[to_cnum(1) to_cnum(0); to_cnum(0) to_cnum(-1)]
        @test isequal(
            S * J * S',
            Coeff[to_cnum(u^2 - v^2) to_cnum(0); to_cnum(0) to_cnum(v^2 - u^2)],
        )
    end

    @testset "inner products and norms" begin
        @variables x y θ

        # `dot` conjugates its first argument.
        @test isequal(dot(to_cnum(im), to_cnum(2)), to_cnum(-2im))
        @test isequal(dot(to_cnum(im), 2), to_cnum(-2im))
        @test isequal(dot(im, to_cnum(2)), to_cnum(-2im))
        @test isequal(dot(Coeff[to_cnum(x), to_cnum(y)], Coeff[to_cnum(x), to_cnum(y)]), to_cnum(x^2 + y^2))

        # `norm` is exact where `abs` is, and refuses where a magnitude has no order.
        @test to_complex(to_cnum(norm(Coeff[to_cnum(3), to_cnum(4)]))) ≈ 5
        @test isone(to_cnum(norm(Coeff[expim(θ)])))
        @test_throws MethodError norm(Coeff[to_cnum(x)])
    end

    @testset "symmetric and hermitian wrappers" begin
        @variables x y

        @test isequal(
            Matrix(Symmetric(Coeff[to_cnum(x) to_cnum(1); to_cnum(1) to_cnum(y)])),
            Coeff[to_cnum(x) to_cnum(1); to_cnum(1) to_cnum(y)],
        )
        @test isequal(
            Matrix(Hermitian(Coeff[to_cnum(x) to_cnum(im); to_cnum(-im) to_cnum(y)])),
            Coeff[to_cnum(x) to_cnum(im); to_cnum(-im) to_cnum(y)],
        )
    end

    @testset "a coefficient broadcasts as a scalar" begin
        @variables x y

        A = Coeff[to_cnum(x) to_cnum(1); to_cnum(2) to_cnum(y)]
        doubled = Coeff[to_cnum(2x) to_cnum(2); to_cnum(4) to_cnum(2y)]

        @test isequal(A .* to_cnum(2), doubled)
        @test isequal(to_cnum(2) .* A, doubled)
        @test isequal(A .+ to_cnum(1), Coeff[to_cnum(x + 1) to_cnum(2); to_cnum(3) to_cnum(y + 1)])
        @test isequal(rmul!(copy(A), to_cnum(2)), doubled)
        @test isequal(lmul!(to_cnum(2), copy(A)), doubled)
    end
end
