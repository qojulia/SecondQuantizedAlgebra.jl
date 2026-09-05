using SecondQuantizedAlgebra
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim

@testset "Generator-derived exact unitary transformations" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    spin = SpinSpace(:spin)
    Sx = Spin(spin, :S, 1)
    Sy = Spin(spin, :S, 2)
    Sz = Spin(spin, :S, 3)

    pauli = PauliSpace(:pauli)
    σx = Pauli(pauli, :σ, 1)
    σy = Pauli(pauli, :σ, 2)
    σz = Pauli(pauli, :σ, 3)

    @variables θ r ϕ ω t
    @variables α::Number

    function equivalent_on(generators, left_transform, right_transform)
        for generator in generators
            @test iszero(
                simplify(
                    conjugate(generator, left_transform) -
                        conjugate(generator, right_transform),
                ),
            )
        end
    end

    @testset "number rotation" begin
        generic = UnitaryTransform(a' * a, θ)
        named = Rotation(a, θ)
        equivalent_on((a, a'), generic, named)
    end

    @testset "Weyl displacement" begin
        G = im * (α * a' - conj(α) * a)
        generic = UnitaryTransform(G, 1)
        named = Displace(a, α)
        equivalent_on((a, a'), generic, named)
    end

    @testset "single-mode squeeze" begin
        G = (im / 2) * (expim(ϕ) * a'^2 - expim(-ϕ) * a^2)
        generic = UnitaryTransform(G, r)
        named = Squeeze(a, r, ϕ)
        equivalent_on((a, a'), generic, named)
    end

    @testset "beam splitter" begin
        G = im * (left' * right - right' * left)
        generic = UnitaryTransform(G, θ)
        named = BeamSplitter(left, right, θ)
        equivalent_on((left, right, left', right'), generic, named)
    end

    @testset "two-mode squeeze" begin
        G = im * (left' * right' - right * left)
        generic = UnitaryTransform(G, r)
        named = TwoModeSqueeze(left, right, r)
        equivalent_on((left, right, left', right'), generic, named)
    end

    @testset "quadrature rotation and squeeze" begin
        rotation_generator = (x^2 + p^2) / 2
        generic_rotation = UnitaryTransform(rotation_generator, θ)
        named_rotation = Rotation(x, p, θ)
        equivalent_on((x, p), generic_rotation, named_rotation)

        squeeze_generator = (x * p + p * x) / 2
        generic_squeeze = UnitaryTransform(squeeze_generator, r)
        named_squeeze = Squeeze(x, p, r)
        equivalent_on((x, p), generic_squeeze, named_squeeze)
    end

    @testset "scaled exact blocks stay symbolic" begin
        equivalent_on(
            (a, a'),
            UnitaryTransform(2 * a' * a, θ),
            Rotation(a, 2θ),
        )
        equivalent_on(
            (x, p),
            UnitaryTransform(x^2 + p^2, θ),
            Rotation(x, p, 2θ),
        )
    end

    @testset "spin and Pauli rotations" begin
        generic_spin = UnitaryTransform(Sz, θ)
        named_spin = Rotation(Sx, 3, θ)
        equivalent_on((Sx, Sy, Sz), generic_spin, named_spin)

        generic_pauli = UnitaryTransform(σz / 2, θ)
        named_pauli = Rotation(σx, 3, θ)
        equivalent_on((σx, σy, σz), generic_pauli, named_pauli)
    end

    @testset "automatic rotating frame" begin
        H0 = ω * a' * a
        frame = RotatingFrame(H0, t)
        named = Rotation(a, ω * t, t)

        equivalent_on((a, a'), frame, named)
        @test iszero(simplify(gauge_term(frame) - gauge_term(named)))
        @test iszero(simplify(transform(H0, frame)))

        spin_frame = RotatingFrame(Sz, t)
        spin_named = Rotation(Sx, 3, t, t)
        equivalent_on((Sx, Sy, Sz), spin_frame, spin_named)
        @test iszero(simplify(gauge_term(spin_frame) - gauge_term(spin_named)))
    end

    @testset "exact refusal" begin
        @test_throws ArgumentError UnitaryTransform(a, θ)
        @test_throws ArgumentError UnitaryTransform(a - a, θ)
        @test_throws ArgumentError UnitaryTransform(a'^2 * a^2, θ)
        @test_throws ArgumentError UnitaryTransform(x^2, θ)
        @test_throws ArgumentError UnitaryTransform(p^2, θ)
        @test_throws ArgumentError UnitaryTransform(a' * a + a + a', θ)
        @test_throws ArgumentError UnitaryTransform(
            a' * a + (im / 2) * (a'^2 - a^2), θ,
        )

        indexed_space = FockSpace(:indexed)
        i = Index(indexed_space, :i, 3, indexed_space)
        indexed_mode = IndexedOperator(Destroy(indexed_space, :c), i)
        indexed_number = indexed_mode' * indexed_mode
        @test_throws ArgumentError UnitaryTransform(indexed_number, θ)
        @test_throws ArgumentError UnitaryTransform(Σ(indexed_number, i), θ)

        three_modes = FockSpace(:one) ⊗ FockSpace(:two) ⊗ FockSpace(:three)
        one = Destroy(three_modes, :one, 1)
        two = Destroy(three_modes, :two, 2)
        three = Destroy(three_modes, :three, 3)
        chain = one' * two + two' * one + two' * three + three' * two
        @test_throws ArgumentError UnitaryTransform(chain, θ)

        @test_throws ArgumentError RotatingFrame((ω + t) * a' * a, t)
    end
end
