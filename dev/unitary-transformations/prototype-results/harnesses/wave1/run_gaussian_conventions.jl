using Test
using SecondQuantizedAlgebra
using Symbolics

include("gaussian_conventions.jl")
using .GaussianConventionsPrototype

import SecondQuantizedAlgebra: Op, QAdd, _CNUM_ZERO, _add_cnum, _conj_cnum,
    _dt, _iszero_cnum, _mul_cnum, _neg_cnum, _substitute_cnum,
    _to_cnum, _to_complex, to_num

function coefficient_arrays_equal(A, B)
    size(A) == size(B) || return false
    all(isequal(a, b) for (a, b) in zip(A, B))
end

function numeric_coefficients(A, substitutions)
    map(A) do coefficient
        _to_complex(_substitute_cnum(coefficient, substitutions))
    end
end

function nonzero_pattern(A)
    [!_iszero_cnum(A[i, j]) for i in axes(A, 1), j in axes(A, 2)]
end

function warm_measure(f; samples = 21)
    f()
    times = Float64[]
    allocations = Int[]
    bytes = Int[]
    for _ in 1:samples
        GC.gc()
        push!(allocations, @allocations(f()))
        push!(bytes, @allocated(f()))
        push!(times, @elapsed(f()))
    end
    middle = (samples + 1) ÷ 2
    return (
        seconds = sort!(times)[middle],
        allocations = sort!(allocations)[middle],
        bytes = sort!(bytes)[middle],
    )
end

@testset "Gaussian convention spike" begin
    two_fock = FockSpace(:left) ⊗ FockSpace(:right)
    a = Destroy(two_fock, :a, 1)
    b = Destroy(two_fock, :b, 2)
    interleaved = interleaved_order(Op[a, b])
    split = split_order(Op[a, b])

    @test interleaved == Op[a, a', b, b']
    @test split == Op[a, b, a', b']
    @test GaussianConventionsPrototype.adjoint_permutation(interleaved) == [2, 1, 4, 3]
    @test GaussianConventionsPrototype.adjoint_permutation(split) == [3, 4, 1, 2]
    @test to_num.(commutator_matrix(interleaved)) ==
        Complex{Num}[0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0]
    @test to_num.(commutator_matrix(split)) ==
        Complex{Num}[0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]

    @variables ω1 ω2 gR gI λR λI μR μI ηR ηI ξR ξI
    g = gR + im * gI
    λ = λR + im * λI
    μ = μR + im * μI
    η = ηR + im * ηI
    ξ = ξR + im * ξI
    H = ω1 * a' * a + ω2 * b' * b + g * a' * b + conj(g) * b' * a +
        (λ / 2) * a'^2 + (conj(λ) / 2) * a^2 +
        μ * a' * b' + conj(μ) * b * a +
        η * a' + conj(η) * a + ξ * b' + conj(ξ) * b

    direct_i = direct_affine_action(H, interleaved)
    scan_i = scanned_affine_action(H, interleaved)
    direct_s = direct_affine_action(H, split)
    scan_s = scanned_affine_action(H, split)
    @test coefficient_arrays_equal(direct_i[1], scan_i[1])
    @test coefficient_arrays_equal(direct_i[2], scan_i[2])
    @test coefficient_arrays_equal(direct_s[1], scan_s[1])
    @test coefficient_arrays_equal(direct_s[2], scan_s[2])

    # Independent hand-derived first two rows in split Nambu ordering.
    expected_first = _to_cnum.([-im * ω1, -im * g, -im * λ, -im * μ])
    expected_second = _to_cnum.([-im * conj(g), -im * ω2, -im * μ, 0])
    @test all(isequal.(direct_s[1][1, :], expected_first))
    @test all(isequal.(direct_s[1][2, :], expected_second))
    @test isequal(direct_s[2][1], _to_cnum(-im * η))
    @test isequal(direct_s[2][2], _to_cnum(-im * ξ))
    outside = Destroy(FockSpace(:outside), :c)
    @test_throws ErrorException scanned_affine_action(H + outside, interleaved)
    @test_throws ErrorException scanned_affine_action(H + a'^2 * a, interleaved)

    values = Dict(
        ω1 => 1.3, ω2 => 2.1, gR => 0.2, gI => -0.1,
        λR => 0.07, λI => 0.03, μR => -0.11, μI => 0.05,
        ηR => 0.4, ηI => -0.2, ξR => -0.15, ξI => 0.08,
    )
    numeric_action = numeric_coefficients(direct_s[1], values)
    numeric_forcing = numeric_coefficients(direct_s[2], values)
    expected_numeric_first = -im .* ComplexF64[
        1.3, 0.2 - 0.1im, 0.07 + 0.03im, -0.11 + 0.05im,
    ]
    @test numeric_action[1, :] ≈ expected_numeric_first atol = 1.0e-14
    @test numeric_forcing[1] ≈ -im * (0.4 - 0.2im) atol = 1.0e-14

    # Ordering is a permutation only.  Split ordering exposes the particle/hole
    # 2x2 blocks, while interleaving exposes site-local 2x2 blocks.
    permutation = [findfirst(isequal(z), interleaved) for z in split]
    @test numeric_action ≈
        numeric_coefficients(direct_i[1], values)[permutation, permutation]

    # The same quadratic Hamiltonian in q=(x1,p1,x2,p2), using
    # a=(x+im*p)/sqrt(2).  This independently pins every sign, 1/2,
    # sqrt(2), and normal-order scalars.  Besides -(w1+w2)/2 from
    # a' a=(x^2+p^2-1)/2, canonicalizing (lambdaI/2)(xp+px) to x-before-p
    # exposes -im*lambdaI/2 while leaving the Hermitian expression unchanged.
    two_phase = PhaseSpace(:qleft) ⊗ PhaseSpace(:qright)
    x1 = Position(two_phase, :x1, 1)
    p1 = Momentum(two_phase, :p1, 1)
    x2 = Position(two_phase, :x2, 2)
    p2 = Momentum(two_phase, :p2, 2)
    qtwo = Op[x1, p1, x2, p2]
    Hq =
        ((ω1 + λR) / 2) * x1^2 + ((ω1 - λR) / 2) * p1^2 +
        (λI / 2) * (x1 * p1 + p1 * x1) +
        (ω2 / 2) * (x2^2 + p2^2) +
        gR * (x1 * x2 + p1 * p2) + gI * (p1 * x2 - x1 * p2) +
        μR * (x1 * x2 - p1 * p2) + μI * (x1 * p2 + p1 * x2) +
        sqrt(2) * (ηR * x1 + ηI * p1 + ξR * x2 + ξI * p2) -
        (ω1 + ω2) / 2
    Aq, bq = direct_affine_action(Hq, qtwo)
    numeric_Aq = numeric_coefficients(Aq, values)
    numeric_bq = numeric_coefficients(bq, values)
    s2 = inv(sqrt(2.0))
    Tlocal = ComplexF64[s2 s2; -im * s2 im * s2]
    T = zeros(ComplexF64, 4, 4)
    T[1:2, 1:2] .= Tlocal
    T[3:4, 3:4] .= Tlocal
    numeric_Ai = numeric_coefficients(direct_i[1], values)
    numeric_bi = numeric_coefficients(direct_i[2], values)
    @test numeric_Aq ≈ T * numeric_Ai / T atol = 2.0e-14
    @test numeric_bq ≈ T * numeric_bi atol = 2.0e-14
    @test isequal(
        scalar_coefficient(Hq),
        _to_cnum(-(ω1 + ω2) / 2 - im * λI / 2),
    )

    # Both orderings materialize the same rules.  Split order needs no U/V block
    # gather; the final Dict and sorted public generators remain site-local.
    split_from_interleaved = direct_i[1][permutation, permutation]
    @test materialize_rows(direct_i[1], interleaved) ==
        materialize_rows(split_from_interleaved, split)
    @test (@inferred commutator_matrix(interleaved)) isa Matrix{SecondQuantizedAlgebra.Coeff}
    @test (@inferred direct_affine_action(H, interleaved)) isa
        Tuple{Matrix{SecondQuantizedAlgebra.Coeff}, Vector{SecondQuantizedAlgebra.Coeff}}
    @test (@inferred scanned_affine_action(H, interleaved)) isa
        Tuple{Matrix{SecondQuantizedAlgebra.Coeff}, Vector{SecondQuantizedAlgebra.Coeff}}
    @test (@inferred materialize_rows(direct_i[1], interleaved)) isa Dict{Op, QAdd}

    @variables t θ r ϕ u v
    phase = PhaseSpace(:phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)
    qbasis = Op[x, p]

    named = [
        ("Fock displacement", Displace(a, u * t + im * v * t^2, t), Op[a, a']),
        ("Fock rotation", Rotation(a, θ * t, t), Op[a, a']),
        ("single-mode squeeze", Squeeze(a, r * t, ϕ * t, t), Op[a, a']),
        ("beamsplitter", Rotation(a, b, θ * t, t), interleaved),
        ("two-mode squeeze", Squeeze(a, b, r * t, t), interleaved),
        ("quadrature displacement", Displace(x, p, u * t, v * t^2, t), qbasis),
        ("quadrature rotation", Rotation(x, p, θ * t, t), qbasis),
        ("quadrature squeeze", Squeeze(x, p, r * t, t), qbasis),
    ]

    lift_residuals = Dict{String, Any}()
    for (name, U, basis) in named
        @test all(iszero, velocity_residual(U, basis, t))
        reconstructed = reconstructed_gauge(U, basis, t)
        difference = simplify(gauge_term(U) - reconstructed)
        @test iszero(operator_part(difference))
        lift_residuals[name] = to_num(scalar_coefficient(difference))
    end

    # Same affine phase-space rotation, inequivalent unitary lifts:
    # exp(-i θ n) and exp(-i θ (n+1/2)).  The maps M,d and all their
    # derivatives agree under the ladder/quadrature coordinate conversion, but
    # their gauges differ by the central term +dot(θ)/2.
    @test iszero(simplify(lift_residuals["Fock rotation"] - θ / 2))
    @test iszero(simplify(lift_residuals["quadrature rotation"]))
    @test iszero(simplify(lift_residuals["Fock displacement"] - u * v * t^2))
    @test iszero(
        simplify(lift_residuals["quadrature displacement"] - u * v * t^2 / 2),
    )

    # Inversion and composition preserve the velocity equation.  Reconstruction
    # remains unique only modulo the scalar lift.
    composite = Rotation(a, θ * t, t) * Displace(a, (u + im * v) * t, t)
    for U in (composite, inv(composite))
        @test all(iszero, velocity_residual(U, Op[a, a'], t))
        delta = simplify(gauge_term(U) - reconstructed_gauge(U, Op[a, a'], t))
        @test iszero(operator_part(delta))
    end

    println("ordering.interleaved.pattern=", nonzero_pattern(direct_i[1]))
    println("ordering.split.pattern=", nonzero_pattern(direct_s[1]))
    for (name, residual) in sort!(collect(lift_residuals); by = first)
        println("lift.", name, "=", residual)
    end

    mode_count = 33
    many_space = ⊗((FockSpace(Symbol(:mode, i)) for i in 1:mode_count)...)
    lowering = Op[
        Destroy(many_space, Symbol(:a, i), i) for i in 1:mode_count
    ]
    many_i = interleaved_order(lowering)
    many_s = split_order(lowering)
    n = length(many_i)
    matrix_i = fill(_CNUM_ZERO, n, n)
    for i in 1:n
        matrix_i[i, i] = _to_cnum(1)
    end
    # Add a sparse nearest-neighbour ladder map so this is not an identity-only
    # special case.  It need not be canonical: this measures only materialization.
    for i in 1:(n - 2)
        matrix_i[i, i + 2] = _to_cnum(1 // 7)
    end
    ps = [findfirst(isequal(z), many_i) for z in many_s]
    matrix_s = matrix_i[ps, ps]
    measure_i = warm_measure(() -> materialize_rows(matrix_i, many_i))
    measure_s = warm_measure(() -> materialize_rows(matrix_s, many_s))
    measure_direct = warm_measure(() -> direct_affine_action(H, interleaved))
    measure_scan = warm_measure(() -> scanned_affine_action(H, interleaved))
    println("materialize.interleaved=", measure_i)
    println("materialize.split=", measure_s)
    println("action.direct_commutators=", measure_direct)
    println("action.single_scan=", measure_scan)
    println(
        "action.nonzero_coefficients=",
        count(nonzero_pattern(direct_i[1])) + count(!_iszero_cnum, direct_i[2]),
    )
end

