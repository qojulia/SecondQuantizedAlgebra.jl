using Test
using LinearAlgebra: Diagonal, I, norm
using CheckConcreteStructs: all_concrete
using SecondQuantizedAlgebra
using Symbolics: @variables, Num
import SecondQuantizedAlgebra: Coeff, DynamicTime, Op, ParamRelation, QAdd, SiteInfo,
    StaticTime, _CNUM_ONE, _coefficient_matrix, _matrix_unit_rules, _nlevel_gauge,
    _rule_qadd, _substitute_cnum, _to_complex, _validated_transform,
    _zero_qadd, expim

# The unitary API is deliberately tested from its closed forms instead of from a second
# implementation of the removed generic exponentiator. Constructor tables pin one useful
# image from every retained family. Algebraic laws are then tested once on a small,
# noncommuting set, and the independent oracles below use only ordinary Julia matrices.

function _matrix_of(
        q::QAdd, representations::Dict{Op, Matrix{ComplexF64}},
        substitutions::AbstractDict = Dict{Num, Float64}()
    )
    dimension = size(first(values(representations)), 1)
    result = zeros(ComplexF64, dimension, dimension)
    identity_matrix = Matrix{ComplexF64}(I, dimension, dimension)
    for (term, coefficient) in q
        reduced = isempty(substitutions) ? coefficient :
            _substitute_cnum(coefficient, substitutions)
        product = copy(identity_matrix)
        for operator in term.ops
            product *= representations[operator]
        end
        result .+= _to_complex(reduced) .* product
    end
    return result
end

_matrix_of(
    o::Op, representations::Dict{Op, Matrix{ComplexF64}},
    substitutions::AbstractDict = Dict{Num, Float64}()
) =
    _matrix_of(1 * o, representations, substitutions)

function _fock_matrices(dimension::Int)
    a = zeros(ComplexF64, dimension, dimension)
    for n in 1:(dimension - 1)
        a[n, n + 1] = sqrt(n)
    end
    return a, Matrix(a')
end

function _transition_matrices(space::NLevelSpace, name::Symbol, dimension::Int)
    representations = Dict{Op, Matrix{ComplexF64}}()
    for i in 1:dimension, j in 1:dimension
        matrix = zeros(ComplexF64, dimension, dimension)
        matrix[i, j] = 1
        representations[Transition(space, name, i, j)] = matrix
    end
    return representations
end

function _warm_allocated(f)
    f()
    GC.gc()
    return @allocations f()
end

@testset "Exact unitary transformations" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :a, 1)
    right = Destroy(two_modes, :b, 2)

    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    spin = SpinSpace(:spin)
    Sx = Spin(spin, :S, 1)
    Sy = Spin(spin, :S, 2)
    Sz = Spin(spin, :S, 3)

    pauli = PauliSpace(:qubit)
    σx = Pauli(pauli, :σ, 1)
    σy = Pauli(pauli, :σ, 2)
    σz = Pauli(pauli, :σ, 3)

    atom = NLevelSpace(:atom, 2)
    σ11 = Transition(atom, :σ, 1, 1)
    σ12 = Transition(atom, :σ, 1, 2)
    σ21 = Transition(atom, :σ, 2, 1)
    σ22 = Transition(atom, :σ, 2, 2)

    @variables θ φ r ϕ ω Ω η g t s dx dp ωd K
    @variables α::Number
    @variables envelope(t)

    @testset "constructor image table" begin
        cases = [
            (
                name = "Fock displacement",
                transform = Displace(a, α),
                images = [a => a + α, a' => a' + conj(α)],
                inverse_images = [a => a - α, a' => a' - conj(α)],
                compound = a' * a,
                compound_image = (a' + conj(α)) * (a + α),
            ),
            (
                name = "Fock rotation",
                transform = Rotation(a, θ),
                images = [a => expim(-θ) * a, a' => expim(θ) * a'],
                inverse_images = [a => expim(θ) * a, a' => expim(-θ) * a'],
                compound = a + a',
                compound_image = expim(-θ) * a + expim(θ) * a',
            ),
            (
                name = "single-mode squeeze",
                transform = Squeeze(a, r, ϕ),
                images = [
                    a => cosh(r) * a + expim(ϕ) * sinh(r) * a',
                    a' => cosh(r) * a' + expim(-ϕ) * sinh(r) * a,
                ],
                inverse_images = [
                    a => cosh(r) * a - expim(ϕ) * sinh(r) * a',
                    a' => cosh(r) * a' - expim(-ϕ) * sinh(r) * a,
                ],
                compound = a * a,
                compound_image =
                    (cosh(r) * a + expim(ϕ) * sinh(r) * a')^2,
            ),
            (
                name = "beamsplitter",
                transform = Rotation(left, right, θ),
                images = [
                    left => cos(θ) * left + sin(θ) * right,
                    right => cos(θ) * right - sin(θ) * left,
                ],
                inverse_images = [
                    left => cos(θ) * left - sin(θ) * right,
                    right => cos(θ) * right + sin(θ) * left,
                ],
                compound = left' * right,
                compound_image =
                    (cos(θ) * left' + sin(θ) * right') *
                    (cos(θ) * right - sin(θ) * left),
            ),
            (
                name = "two-mode squeeze",
                transform = Squeeze(left, right, r),
                images = [
                    left => cosh(r) * left + sinh(r) * right',
                    right => cosh(r) * right + sinh(r) * left',
                ],
                inverse_images = [
                    left => cosh(r) * left - sinh(r) * right',
                    right => cosh(r) * right - sinh(r) * left',
                ],
                compound = left * right,
                compound_image =
                    (cosh(r) * left + sinh(r) * right') *
                    (cosh(r) * right + sinh(r) * left'),
            ),
            (
                name = "quadrature displacement",
                transform = Displace(x, p, dx, dp),
                images = [x => x + dx, p => p + dp],
                inverse_images = [x => x - dx, p => p - dp],
                compound = x * p,
                compound_image = (x + dx) * (p + dp),
            ),
            (
                name = "quadrature rotation",
                transform = Rotation(x, p, θ),
                images = [
                    x => cos(θ) * x + sin(θ) * p,
                    p => cos(θ) * p - sin(θ) * x,
                ],
                inverse_images = [
                    x => cos(θ) * x - sin(θ) * p,
                    p => cos(θ) * p + sin(θ) * x,
                ],
                compound = x * x + p * p,
                compound_image = x * x + p * p,
            ),
            (
                name = "quadrature squeeze",
                transform = Squeeze(x, p, r),
                images = [x => exp(r) * x, p => inv(exp(r)) * p],
                inverse_images = [x => inv(exp(r)) * x, p => exp(r) * p],
                compound = x * p + p * x,
                compound_image = x * p + p * x,
            ),
            (
                name = "spin rotation",
                transform = Rotation(Sx, 3, θ),
                images = [
                    Sx => cos(θ) * Sx - sin(θ) * Sy,
                    Sy => cos(θ) * Sy + sin(θ) * Sx,
                    Sz => 1 * Sz,
                ],
                inverse_images = [
                    Sx => cos(θ) * Sx + sin(θ) * Sy,
                    Sy => cos(θ) * Sy - sin(θ) * Sx,
                    Sz => 1 * Sz,
                ],
                compound = Sx * Sy,
                compound_image =
                    (cos(θ) * Sx - sin(θ) * Sy) *
                    (cos(θ) * Sy + sin(θ) * Sx),
            ),
            (
                name = "Pauli rotation",
                transform = Rotation(σx, 2, θ),
                images = [
                    σz => cos(θ) * σz - sin(θ) * σx,
                    σx => cos(θ) * σx + sin(θ) * σz,
                    σy => 1 * σy,
                ],
                inverse_images = [
                    σz => cos(θ) * σz + sin(θ) * σx,
                    σx => cos(θ) * σx - sin(θ) * σz,
                    σy => 1 * σy,
                ],
                compound = σx * σz,
                compound_image =
                    (cos(θ) * σx + sin(θ) * σz) *
                    (cos(θ) * σz - sin(θ) * σx),
            ),
            (
                name = "N-level swap",
                transform = Rotation(σ12, [0 1; 1 0]),
                images = [σ11 => 1 * σ22, σ12 => 1 * σ21],
                inverse_images = [σ11 => 1 * σ22, σ12 => 1 * σ21],
                compound = σ11 + 2σ22,
                compound_image = σ22 + 2σ11,
            ),
        ]

        for case in cases
            @testset "$(case.name)" begin
                @test isempty(gauge_term(case.transform))
                for (generator, expected) in case.images
                    @test isequal(simplify(conjugate(generator, case.transform)), simplify(expected))
                end
                inverse = inv(case.transform)
                for (generator, expected) in case.inverse_images
                    @test isequal(simplify(conjugate(generator, inverse)), simplify(expected))
                end
                @test isequal(
                    simplify(conjugate(case.compound, case.transform)),
                    simplify(case.compound_image),
                )
            end
        end
    end

    @testset "all spin axes are complete" begin
        spin_images = (
            (
                Sx => 1 * Sx,
                Sy => cos(θ) * Sy - sin(θ) * Sz,
                Sz => cos(θ) * Sz + sin(θ) * Sy,
            ),
            (
                Sy => 1 * Sy,
                Sz => cos(θ) * Sz - sin(θ) * Sx,
                Sx => cos(θ) * Sx + sin(θ) * Sz,
            ),
            (
                Sz => 1 * Sz,
                Sx => cos(θ) * Sx - sin(θ) * Sy,
                Sy => cos(θ) * Sy + sin(θ) * Sx,
            ),
        )
        for axis in 1:3
            U = Rotation(Sx, axis, θ)
            @test length(generators(U)) == 3
            for image in spin_images[axis]
                @test isequal(conjugate(image.first, U), image.second)
            end
        end
    end

    @testset "timed closed forms and gauges" begin
        moving_displacement = Displace(a, (η + im * ω * t) * t, t)
        @test isequal(conjugate(a, moving_displacement), a + (η + im * ω * t) * t)
        @test isequal(
            gauge_term(moving_displacement),
            -im * (η + 2im * ω * t) * a' +
                im * (η - 2im * ω * t) * a + η * ω * t^2,
        )

        moving_rotation = Rotation(a, ω * t, t)
        @test isequal(conjugate(a, moving_rotation), expim(-ω * t) * a)
        @test isequal(gauge_term(moving_rotation), -ω * a' * a)
        @test isequal(transform(Ω * a' * a, moving_rotation), (Ω - ω) * a' * a)

        moving_squeeze = Squeeze(a, r * t, 0, t)
        @test isequal(
            simplify(gauge_term(moving_squeeze)),
            (im * r * (1 // 2)) * (a * a - a' * a'),
        )
        moving_phase_squeeze = Squeeze(a, r, ω * t, t)
        @test isequal(gauge_term(moving_phase_squeeze), adjoint(gauge_term(moving_phase_squeeze)))

        moving_beamsplitter = Rotation(left, right, θ * t, t)
        expected_beamsplitter_gauge =
            -im * θ * (left' * right - right' * left)
        @test isequal(gauge_term(moving_beamsplitter), expected_beamsplitter_gauge)
        @test isequal(
            gauge_term(Rotation(left', right', θ * t, t)), expected_beamsplitter_gauge,
        )

        moving_two_mode_squeeze = Squeeze(left, right, r * t, t)
        expected_two_mode_gauge = -im * r * (left' * right' - right * left)
        @test isequal(gauge_term(moving_two_mode_squeeze), expected_two_mode_gauge)
        @test isequal(
            gauge_term(Squeeze(left', right', r * t, t)), expected_two_mode_gauge,
        )

        moving_quadratures = Displace(x, p, η * t, ω * t^2, t)
        @test isequal(
            gauge_term(moving_quadratures),
            2ω * t * x - η * p + (η * ω * t^2) * (1 // 2),
        )
        @test isequal(
            gauge_term(Rotation(x, p, θ * t, t)),
            -(θ * (1 // 2)) * (x * x + p * p),
        )
        @test isequal(
            gauge_term(Squeeze(x, p, r * t, t)),
            -(r * (1 // 2)) * (x * p + p * x),
        )

        @test isequal(gauge_term(Rotation(Sx, 3, θ * t, t)), -θ * Sz)
        @test isequal(
            gauge_term(Rotation(σx, 3, θ * t, t)), -(θ * (1 // 2)) * σz,
        )

        phase_matrix = Coeff[expim(-ω * t) 0; 0 expim(ω * t)]
        timed_levels = Rotation(σ12, phase_matrix, t)
        @test isequal(
            conjugate(σ12, timed_levels), expim(ω * t) * expim(ω * t) * σ12,
        )
        @test isequal(gauge_term(timed_levels), -ω * σ11 + ω * σ22)
        @test isequal(
            transform(Ω * (σ11 + σ22), timed_levels),
            (Ω - ω) * σ11 + (Ω + ω) * σ22,
        )
    end

    @testset "automatic bounded Fock displacement" begin
        Hpaper = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')
        static_reference = ω * a' * a + η * (a + a') + g
        static_displacement = Displace(a, static_reference)
        @test static_displacement isa UnitaryTransform{StaticTime}
        @test isequal(conjugate(a, static_displacement), a - η / ω)
        @test transform(static_reference, static_displacement) ==
            ω * a' * a - η^2 / ω + g

        U = Displace(a, Hpaper, t)
        transformed_reference = transform(Hpaper, U)
        operator_terms = [
            (term.ops, coefficient) for
                (term, coefficient) in transformed_reference if !isempty(term.ops)
        ]
        @test operator_terms == [(Op[a', a], SecondQuantizedAlgebra._to_cnum(ω))]

        # Independent numerical comparison with Eq. (2)'s displacement in
        # Phys. Rev. A 110, 042411. Test both sides of the oscillator frequency.
        representations = Dict{Op, Matrix{ComplexF64}}(
            a => zeros(ComplexF64, 2, 2), a' => zeros(ComplexF64, 2, 2),
        )
        for values in (
                Dict(ω => 5.0, ωd => 2.0, Ω => 0.3, t => 0.7),
                Dict(ω => 1.25, ωd => 2.5, Ω => -0.4, t => 0.2),
            )
            actual = _matrix_of(conjugate(a, U) - a, representations, values)[1, 1]
            ωv, ωdv, Ωv, tv = values[ω], values[ωd], values[Ω], values[t]
            expected =
                im * Ωv / (2 * (ωdv - ωv)) * exp(-im * ωdv * tv) -
                im * Ωv / (2 * (ωdv + ωv)) * exp(im * ωdv * tv)
            @test actual ≈ expected atol = 1.0e-13
        end

        constant = Displace(a, static_reference, t)
        @test isequal(conjugate(a, constant), a - η / ω)
        constant_terms = [
            (term.ops, coefficient) for
                (term, coefficient) in transform(static_reference, constant) if
                !isempty(term.ops)
        ]
        @test constant_terms ==
            [(Op[a', a], SecondQuantizedAlgebra._to_cnum(ω))]

        complex_reference = ω * a' * a + α * a' + conj(α) * a
        complex_displacement = Displace(a, complex_reference, t)
        @test isequal(conjugate(a, complex_displacement), a - α / ω)

        multitone_drive = η * cos(ωd * t) + g * sin(2ωd * t)
        multitone_reference = ω * a' * a + multitone_drive * (a + a')
        multitone = Displace(a, multitone_reference, t)
        multitone_terms = [
            (term.ops, coefficient) for
                (term, coefficient) in transform(multitone_reference, multitone) if
                !isempty(term.ops)
        ]
        @test multitone_terms ==
            [(Op[a', a], SecondQuantizedAlgebra._to_cnum(ω))]
        @test isequal(conjugate(a', multitone), conjugate(a, multitone)')

        from_creation = Displace(a', Hpaper, t)
        @test isequal(conjugate(a, from_creation), conjugate(a, U))

        Kerr = (K / 2) * a'^2 * a^2
        full = Hpaper + Kerr
        @test isequal(
            transform(full, U), transform(Hpaper, U) + conjugate(Kerr, U),
        )

        other = Destroy(FockSpace(:automatic_other), :b)
        @test_throws ArgumentError Displace(a, Hpaper + K * a'^2 * a^2, t)
        @test_throws ArgumentError Displace(a, Hpaper + g * other, t)
        @test_throws ArgumentError Displace(a, ω * a' * a + η * a', t)
        @test_throws ArgumentError Displace(
            a, (ω + g * cos(ωd * t)) * a' * a + η * (a + a'), t,
        )
        @test_throws ArgumentError Displace(
            a, ω * a' * a + envelope * (a + a'), t,
        )
        nonlinear_phase = expim(t^2)
        @test_throws ArgumentError Displace(
            a, ω * a' * a + nonlinear_phase * a' + conj(nonlinear_phase) * a, t,
        )
        resonant = expim(-ω * t)
        @test_throws ArgumentError Displace(
            a, ω * a' * a + resonant * a' + conj(resonant) * a, t,
        )
        @test_throws ArgumentError Displace(a, η * (a + a'), t)
        @test_throws ArgumentError Displace(a, η * (a + a'))

        # A possible symbolic resonance remains an exact quotient by design.
        @test Displace(a, Hpaper, t) isa UnitaryTransform{DynamicTime}
    end

    @testset "automatic bounded quadrature displacement" begin
        static_reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 + η * x + dx * p
        static = Displace(x, p, static_reference)
        @test static isa UnitaryTransform{StaticTime}
        determinant = ω * Ω - g^2
        @test isequal(conjugate(x, static), x + (g * dx - Ω * η) / determinant)
        @test isequal(conjugate(p, static), p + (g * η - ω * dx) / determinant)
        static_operator_terms = [
            term.ops for (term, _) in simplify(transform(static_reference, static)) if
                !isempty(term.ops)
        ]
        @test Set(static_operator_terms) == Set((Op[x, x], Op[x, p], Op[p, p]))

        isotropic_reference =
            (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        isotropic = Displace(x, p, isotropic_reference, t)
        @test isotropic isa UnitaryTransform{DynamicTime}
        isotropic_terms = [
            term.ops for (term, _) in transform(isotropic_reference, isotropic) if
                !isempty(term.ops)
        ]
        @test Set(isotropic_terms) == Set((Op[x, x], Op[p, p]))

        representations = Dict{Op, Matrix{ComplexF64}}(
            x => zeros(ComplexF64, 2, 2), p => zeros(ComplexF64, 2, 2),
        )
        for values in (
                Dict(ω => 4.0, ωd => 1.5, η => 0.3, t => 0.7),
                Dict(ω => 1.25, ωd => 2.5, η => -0.4, t => 0.2),
            )
            shifts = ComplexF64[0, 0]
            for sign in (-1, 1)
                frequency = sign * values[ωd]
                matrix = ComplexF64[
                    values[ω] im * frequency
                    -im * frequency values[ω]
                ]
                force = ComplexF64[
                    values[η] * exp(im * frequency * values[t]) / 2, 0,
                ]
                shifts .-= matrix \ force
            end
            actual_x = _matrix_of(
                conjugate(x, isotropic) - x, representations, values,
            )[1, 1]
            actual_p = _matrix_of(
                conjugate(p, isotropic) - p, representations, values,
            )[1, 1]
            @test actual_x ≈ shifts[1] atol = 1.0e-13
            @test actual_p ≈ shifts[2] atol = 1.0e-13
        end

        multitone_reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 +
            (η * cos(ωd * t) + dx * sin(2ωd * t)) * x +
            dp * cos(3ωd * t) * p
        multitone = Displace(x, p, multitone_reference, t)
        multitone_terms = [
            term.ops for (term, _) in transform(multitone_reference, multitone) if
                !isempty(term.ops)
        ]
        @test Set(multitone_terms) == Set((Op[x, x], Op[x, p], Op[p, p]))

        other_phase = PhaseSpace(:automatic_selected) ⊗ PhaseSpace(:automatic_other)
        other_x = Position(other_phase, :other_x, 2)
        nonlinear_phase = expim(t^2)
        @test_throws ArgumentError Displace(p, x, isotropic_reference, t)
        @test_throws ArgumentError Displace(
            x, p, isotropic_reference + other_x, t,
        )
        @test_throws ArgumentError Displace(
            x, p, isotropic_reference + K * x^3, t,
        )
        @test_throws ArgumentError Displace(
            x, p, isotropic_reference + im * x, t,
        )
        @test_throws ArgumentError Displace(
            x, p,
            ((ω + g * cos(ωd * t)) / 2) * x^2 + (Ω / 2) * p^2 + η * x,
            t,
        )
        @test_throws ArgumentError Displace(
            x, p, (ω / 2) * (x^2 + p^2) + envelope * x, t,
        )
        @test_throws ArgumentError Displace(
            x, p,
            (ω / 2) * (x^2 + p^2) + (nonlinear_phase + conj(nonlinear_phase)) * x,
            t,
        )
        @test_throws ArgumentError Displace(
            x, p, (ω / 2) * (x^2 + p^2) + η * cos(ω * t) * x, t,
        )
        @test_throws ArgumentError Displace(x, p, (ω / 2) * x^2 + η * p)
        nonunique = Displace(x, p, (ω / 2) * x^2)
        @test isequal(conjugate(x + p, nonunique), x + p)

        # A possible symbolic determinant resonance remains an exact quotient.
        @test Displace(x, p, isotropic_reference, t) isa
            UnitaryTransform{DynamicTime}
    end

    @testset "cross-cutting algebraic laws" begin
        U = Rotation(a, θ) * Squeeze(a, r, ϕ)
        A = a + 2a'

        homomorphism_transform = Rotation(a, θ)
        @test isequal(
            simplify(conjugate(a * a', homomorphism_transform)),
            simplify(
                conjugate(a, homomorphism_transform) *
                    conjugate(a', homomorphism_transform),
            ),
        )
        @test isequal(conjugate(adjoint(A), U), adjoint(conjugate(A, U)))
        @test isequal(simplify(conjugate(conjugate(A, U), inv(U))), A)
        @test isequal(inv(inv(U)).rules, U.rules)
        @test isequal(adjoint(U).rules, inv(U).rules)

        first = Rotation(a, θ)
        second = Squeeze(a, r, ϕ)
        composed = first * second
        @test isequal(
            simplify(conjugate(A, composed)),
            simplify(conjugate(conjugate(A, first), second)),
        )
        @test !isequal(
            simplify(conjugate(A, first * second)),
            simplify(conjugate(A, second * first)),
        )

        # Composition must retain rules that belong only to either transform. This
        # exercises independent modes, rather than the overlapping Fock rules above.
        independent = Destroy(FockSpace(:independent), :c)
        disjoint = first * Rotation(independent, φ)
        @test length(generators(disjoint)) == 4
        @test isequal(
            simplify(conjugate(a + independent, disjoint)),
            simplify(expim(-θ) * a + expim(-φ) * independent),
        )

        # Pin each composition mode directly; the explicit pairs keep the application
        # order readable and avoid relying on commutativity in the law above.
        composition_pairs = (
            (Rotation(a, θ), Squeeze(a, r)),
            (Rotation(a, θ), Squeeze(a, r * t, 0, t)),
            (Rotation(a, ω * t, t), Squeeze(a, r)),
            (Rotation(a, ω * t, t), Squeeze(a, r * t, 0, t)),
        )
        for (outer, inner) in composition_pairs
            @test isequal(
                simplify(conjugate(A, outer * inner)),
                simplify(conjugate(conjugate(A, outer), inner)),
            )
        end

        timed_first = Rotation(a, ω * t, t)
        timed_second = Displace(a, η * t, t)
        timed_composed = timed_first * timed_second
        H = Ω * a' * a + g * (a + a')
        @test isequal(
            simplify(transform(H, timed_composed)),
            simplify(transform(transform(H, timed_first), timed_second)),
        )
        @test isequal(
            gauge_term(timed_composed),
            conjugate(gauge_term(timed_first), timed_second) + gauge_term(timed_second),
        )
        @test isequal(
            gauge_term(inv(timed_composed)),
            -conjugate(gauge_term(timed_composed), inv(timed_composed)),
        )

        static_H = Ω * a' * a + g * (a + a')
        @test isequal(transform(static_H, U), conjugate(static_H, U))
        @test isequal(transform(static_H, U), adjoint(transform(static_H, U)))
        @test isequal(transform(static_H, timed_composed), adjoint(transform(static_H, timed_composed)))
    end

    @testset "independent direct-matrix oracles" begin
        @testset "low-excitation Fock displacement" begin
            dimension = 24
            lowering, raising = _fock_matrices(dimension)
            representations = Dict{Op, Matrix{ComplexF64}}(a => lowering, a' => raising)
            amplitude = 0.04 + 0.03im
            displacement = exp(amplitude * raising - conj(amplitude) * lowering)
            exact = displacement' * lowering * displacement
            symbolic = _matrix_of(conjugate(a, Displace(a, amplitude)), representations)
            # The finite truncation violates the CCR only at its top state. The first six
            # levels are an exponentially converged oracle for the infinite Fock result.
            @test norm(exact[1:6, 1:6] - symbolic[1:6, 1:6]) < 1.0e-12
        end

        @testset "Pauli rotation" begin
            sx = ComplexF64[0 1; 1 0]
            sy = ComplexF64[0 -im; im 0]
            sz = ComplexF64[1 0; 0 -1]
            representations = Dict{Op, Matrix{ComplexF64}}(σx => sx, σy => sy, σz => sz)
            angle = 0.37
            matrix_unitary = exp(-im * angle * sz / 2)
            oracle = matrix_unitary' * sx * matrix_unitary
            symbolic = _matrix_of(
                conjugate(σx, Rotation(σx, 3, θ)), representations, Dict(θ => angle),
            )
            @test norm(oracle - symbolic) < 1.0e-14
        end

        @testset "finite N-level contraction" begin
            W = Complex{Rational{Int}}[3 // 5 4im // 5; 4im // 5 3 // 5]
            transform_matrix = Rotation(σ12, W)
            representations = _transition_matrices(atom, :σ, 2)
            for operator in (σ11, σ12, σ21, σ22)
                input = representations[operator]
                oracle = ComplexF64.(W)' * input * ComplexF64.(W)
                symbolic = _matrix_of(conjugate(operator, transform_matrix), representations)
                @test norm(oracle - symbolic) < 1.0e-14
            end
        end

        @testset "timed N-level gauge derivative" begin
            frequency = 1.7
            instant = 0.31
            step = 1.0e-6
            Wfun(τ) = Diagonal(ComplexF64[exp(-im * frequency * τ), exp(im * frequency * τ)])
            derivative = (Wfun(instant + step) - Wfun(instant - step)) / (2step)
            oracle = im * derivative' * Wfun(instant)

            symbolic_matrix = Coeff[expim(-ω * t) 0; 0 expim(ω * t)]
            U = Rotation(σ12, symbolic_matrix, t)
            representations = _transition_matrices(atom, :σ, 2)
            symbolic = _matrix_of(
                gauge_term(U), representations, Dict(ω => frequency, t => instant),
            )
            @test norm(oracle - symbolic) < 1.0e-9
        end
    end

    @testset "exact refusal and invalid inputs" begin
        @test_throws ArgumentError Displace(σx, 1)
        @test_throws ArgumentError Rotation(x, θ)
        @test_throws ArgumentError Squeeze(Sx, r)
        @test_throws ArgumentError Rotation(x, x, θ)
        @test_throws ArgumentError Rotation(left, left, θ)
        @test_throws ArgumentError Rotation(Sx, 0, θ)
        @test_throws ArgumentError Rotation(Sx, 4, θ)

        @test_throws MethodError Rotation(a, 1 + im)
        @test_throws MethodError Squeeze(a, 1 + im)
        @test_throws MethodError Displace(x, p, 1 + im, 0)
        @test_throws MethodError Rotation(a, ω * t, 1.0)
        @test_throws MethodError Rotation(a, ω * t, 1 + im)

        @test_throws ArgumentError Rotation(a, ω * t, t) * Rotation(a, ω * s, s)
        @test_throws ArgumentError Rotation(a, ω * t) * Rotation(a, ω * t, t)

        lone_rule = Dict{Op, QAdd}(a => _rule_qadd((_CNUM_ONE, Op[a])))
        @test_throws ArgumentError _validated_transform(
            lone_rule, copy(lone_rule), _zero_qadd(), StaticTime(), ParamRelation[],
        )
        @test_throws MethodError UnitaryTransform(a' * a, θ)

        @test_throws ArgumentError Rotation(σ12, [1 0 0; 0 1 0])
        @test_throws ArgumentError Rotation(σ12, [1 0; 0 2])
        @test_throws ArgumentError Rotation(σ12, [cos(θ) -sin(φ); sin(θ) cos(φ)])
        @test_throws MethodError Rotation(σ12, [0 1; 1 0], t, _zero_qadd())

        collective_space = CollectiveNLevelSpace(:ensemble, 2)
        collective = CollectiveTransition(collective_space, :J, 1, 2)
        @test_throws ArgumentError Rotation(collective, [0 1; 1 0])
        undisplaced = Displace(a, a' * a)
        @test undisplaced isa UnitaryTransform{StaticTime}
        @test iszero(conjugate(a, undisplaced) - a)
        @test_throws MethodError Displace(x, p, a' * a, 0)

        @test !isdefined(SecondQuantizedAlgebra, :RotatingFrame)
        @test !isdefined(SecondQuantizedAlgebra, :DressedFrame)
        @test !isdefined(SecondQuantizedAlgebra, :Bogoliubov)
        @test !isdefined(SecondQuantizedAlgebra, :is_canonical)
        @test !isdefined(SecondQuantizedAlgebra, :constraints)
    end

    @testset "defensive generator coverage" begin
        U = Rotation(a, θ)
        listed = generators(U)
        @test listed == sort(listed)
        @test eltype(listed) === Op
        empty!(listed)
        @test length(generators(U)) == 2
        @test length(U.generators) == 2

        # A transform on a site may leave unrelated sites untouched, but it must never
        # emit a half-transformed expression for a covered generator.
        other = Destroy(FockSpace(:other), :b)
        @test isequal(conjugate(a + other, U), expim(-θ) * a + other)
    end

    @testset "concrete storage and inference" begin
        static = Rotation(a, θ)
        timed = Rotation(a, ω * t, t)
        for U in (static, timed)
            transform_type = typeof(U)
            @test isconcretetype(transform_type)
            @test all_concrete(transform_type; verbose = false)
            @test all(isconcretetype, fieldtypes(transform_type))
            @test !(Any in fieldtypes(transform_type))
            @test all(field -> !(field <: Function), fieldtypes(transform_type))
            @test U.rules isa Dict{Op, QAdd}
            @test U.inverse_rules isa Dict{Op, QAdd}
            @test U.generators isa Vector{Op}
            @test U.sites isa Vector{SiteInfo}
            @test U.relations isa Vector{ParamRelation}
        end
        @test static isa UnitaryTransform{StaticTime}
        @test timed isa UnitaryTransform{DynamicTime}

        @test (@inferred Displace(a, α)) isa UnitaryTransform{StaticTime}
        @test (@inferred Displace(a, η * t, t)) isa UnitaryTransform{DynamicTime}
        automatic_reference = ω * a' * a + η * cos(Ω * t) * (a + a')
        @test (@inferred Displace(a, ω * a' * a + η * (a + a'))) isa
            UnitaryTransform{StaticTime}
        @test (@inferred Displace(a, automatic_reference, t)) isa
            UnitaryTransform{DynamicTime}
        quadrature_reference =
            (ω / 2) * (x^2 + p^2) + η * cos(Ω * t) * x + g * sin(Ω * t) * p
        @test (@inferred Displace(x, p, (ω / 2) * (x^2 + p^2) + η * x)) isa
            UnitaryTransform{StaticTime}
        @test (@inferred Displace(x, p, quadrature_reference, t)) isa
            UnitaryTransform{DynamicTime}
        @test (@inferred Rotation(a, θ)) isa UnitaryTransform{StaticTime}
        @test (@inferred Rotation(a, ω * t, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Squeeze(a, r, ϕ)) isa UnitaryTransform{StaticTime}
        @test (@inferred Squeeze(a, r * t, ϕ, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Rotation(left, right, θ)) isa UnitaryTransform{StaticTime}
        @test (@inferred Rotation(left, right, θ * t, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Squeeze(left, right, r)) isa UnitaryTransform{StaticTime}
        @test (@inferred Squeeze(left, right, r * t, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Displace(x, p, dx, dp)) isa UnitaryTransform{StaticTime}
        @test (@inferred Displace(x, p, dx * t, dp * t, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Rotation(x, p, θ)) isa UnitaryTransform{StaticTime}
        @test (@inferred Squeeze(x, p, r)) isa UnitaryTransform{StaticTime}
        @test (@inferred Rotation(Sx, 3, θ)) isa UnitaryTransform{StaticTime}
        @test (@inferred Rotation(Sx, 3, θ * t, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred Rotation(σ12, [0 1; 1 0])) isa UnitaryTransform{StaticTime}
        @test Rotation(σ12, Diagonal([1, 1])) isa UnitaryTransform{StaticTime}

        level_phases = Coeff[expim(-ω * t) 0; 0 expim(ω * t)]
        @test (@inferred Rotation(σ12, level_phases, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred _coefficient_matrix([0 1; 1 0])) isa Matrix{Coeff}
        coefficients = _coefficient_matrix([0 1; 1 0])
        @test (@inferred _matrix_unit_rules(σ12, coefficients)) isa Dict{Op, QAdd}
        @test (@inferred _nlevel_gauge(σ12, level_phases, t)) isa QAdd

        @test (@inferred conjugate(a, static)) isa QAdd
        @test (@inferred conjugate(a' * a, static)) isa QAdd
        @test (@inferred transform(a' * a, static)) isa QAdd
        @test (@inferred transform(a' * a, timed)) isa QAdd
        @test (@inferred inv(static)) isa UnitaryTransform{StaticTime}
        @test (@inferred inv(timed)) isa UnitaryTransform{DynamicTime}
        @test (@inferred static * Squeeze(a, r)) isa UnitaryTransform{StaticTime}
        @test (@inferred static * Squeeze(a, r * t, 0, t)) isa UnitaryTransform{DynamicTime}
        @test (@inferred timed * Squeeze(a, r)) isa UnitaryTransform{DynamicTime}
        @test (@inferred timed * Squeeze(a, r * t, 0, t)) isa UnitaryTransform{DynamicTime}
    end

    @testset "allocation regression gates" begin
        # These ceilings are the pre-refactor medians plus at most 20% headroom. Warm-up is
        # outside the measured region so the checks describe steady-state API costs.
        static = Rotation(a, θ)
        timed = Rotation(a, ω * t, t)
        static_second = Rotation(a, φ)
        timed_second = Rotation(a, η * t, t)
        observable = a' * a

        @test _warm_allocated(() -> Rotation(a, θ)) <= 144
        @test _warm_allocated(() -> conjugate(a, static)) <= 40
        @test _warm_allocated(() -> conjugate(observable, static)) <= 41
        @test _warm_allocated(() -> transform(observable, timed)) <= 57
        @test _warm_allocated(() -> inv(static)) <= 24
        @test _warm_allocated(() -> static * static_second) <= 112
        @test _warm_allocated(() -> static * timed_second) <= 112
        @test _warm_allocated(() -> timed * static_second) <= 112
        @test _warm_allocated(() -> timed * timed_second) <= 112

        automatic_static = ω * a' * a + η * (a + a')
        automatic_reference = ω * a' * a + η * cos(Ω * t) * (a + a')
        quadrature_static = (ω / 2) * (x^2 + p^2) + η * x + g * p
        quadrature_timed =
            (ω / 2) * (x^2 + p^2) + η * cos(Ω * t) * x + g * sin(Ω * t) * p
        # Julia 1.12 warm medians are 145, 4,556, 232, and 8,211 allocations;
        # every ceiling retains no more than 20% headroom.
        @test _warm_allocated(() -> Displace(a, automatic_static)) <= 174
        @test _warm_allocated(() -> Displace(a, automatic_reference, t)) <= 5100
        @test _warm_allocated(() -> Displace(x, p, quadrature_static)) <= 278
        @test _warm_allocated(() -> Displace(x, p, quadrature_timed, t)) <= 9853
    end
end
