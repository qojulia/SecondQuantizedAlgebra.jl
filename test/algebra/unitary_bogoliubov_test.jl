using SecondQuantizedAlgebra
using LinearAlgebra: exp
using QuantumOpticsBase: FockBasis
using Test
using Symbolics: @variables

@testset "Exact bosonic Bogoliubov transformations" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)

    @variables r ϕ θ

    @testset "single-mode squeezing is a Bogoliubov map" begin
        S = [
            cosh(r) exp(im * ϕ) * sinh(r)
            exp(-im * ϕ) * sinh(r) cosh(r)
        ]
        raw = Bogoliubov(a, S)
        named = Squeeze(a, r, ϕ)

        @test raw isa UnitaryTransform
        for op in (a, a')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
            @test iszero(
                simplify(conjugate(op, inv(raw)) - conjugate(op, inv(named))),
            )
        end
        @test iszero(
            simplify(conjugate(a, raw) - cosh(r) * a - exp(im * ϕ) * sinh(r) * a'),
        )
    end

    @testset "passive two-mode mixing uses the same Nambu core" begin
        U = [cos(θ) sin(θ); -sin(θ) cos(θ)]
        V = [0 0; 0 0]
        raw = Bogoliubov((left, right), U, V)
        named = Rotation(left, right, θ)

        @test raw isa UnitaryTransform
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
    end

    @testset "two-mode squeezing uses the same Nambu core" begin
        U = [cosh(r) 0; 0 cosh(r)]
        V = [0 sinh(r); sinh(r) 0]
        raw = Bogoliubov(Op[left, right], U, V)
        named = Squeeze(left, right, r)

        @test raw isa UnitaryTransform
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
    end

    @testset "canonicality is a caller precondition" begin
        @variables u::Number v::Number
        raw = @inferred Bogoliubov(a, [u v; conj(v) conj(u)])
        @test raw isa UnitaryTransform
        @test iszero(simplify(conjugate(a, raw) - u * a - v * a'))

        # Scalar substitution recompiles the affine map. The caller remains responsible for
        # preserving the Bogoliubov contract when substituting arbitrary parameter values.
        resolved = @inferred substitute(raw, Dict(u => 5 // 3, v => 4 // 3))
        @test resolved isa UnitaryTransform
        @test iszero(
            simplify(conjugate(a, resolved) - (5 // 3) * a - (4 // 3) * a'),
        )

        # Canonicality is not a second validation layer. A structurally valid matrix is
        # interpreted under the documented canonicality precondition even if the caller
        # violates it. In that case only the supplied forward map is meaningful.
        assumed = @inferred Bogoliubov(a, [2 0; 0 2])
        @test iszero(simplify(conjugate(a, assumed) - 2 * a))
    end

    @testset "structural validation" begin
        @test_throws ArgumentError Bogoliubov(Op[], Matrix{Int}(undef, 0, 0))
        @test_throws ArgumentError Bogoliubov(a, [1 0 0; 0 1 0])
        @test_throws ArgumentError Bogoliubov((left, left), [1 0; 0 1], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1; 0 0], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1], [0 0 0; 0 0 0])

        identity_map = @inferred Bogoliubov(a, [1 0; 0 1])
        @test iszero(simplify(conjugate(a, identity_map) - a))
        @test iszero(simplify(conjugate(a, inv(identity_map)) - a))

        complex_phase = @inferred Bogoliubov(
            a, ComplexF64[im 0; 0 -im],
        )
        @test iszero(simplify(conjugate(a, complex_phase) - im * a))
    end

    @testset "raw Bogoliubov map has an independent numeric oracle" begin
        dimension = 24
        lowering = zeros(ComplexF64, dimension, dimension)
        for n in 2:dimension
            lowering[n - 1, n] = sqrt(n - 1)
        end

        c = 401 // 399
        s = 40 // 399
        S = [c s; s c]
        raw = Bogoliubov(a, S)
        actual = Matrix(to_numeric(conjugate(a, raw), FockBasis(dimension)).data)
        r0 = asinh(Float64(s))
        unitary = exp((r0 / 2) * ((lowering')^2 - lowering^2))
        expected = unitary' * lowering * unitary
        @test actual[1:6, 1:6] ≈ expected[1:6, 1:6] atol = 1.0e-10
    end
end
