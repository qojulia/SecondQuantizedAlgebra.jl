include("exact_response_solver.jl")

using .ExactResponseSolverPrototype
using SecondQuantizedAlgebra
using LinearAlgebra
using Statistics
using Symbolics
using Test

import SecondQuantizedAlgebra: CNum, _CNUM_ZERO, _to_cnum, _iszero_cnum, to_num

cnum_matrix(x) = Matrix{CNum}(_to_cnum.(x))
cnum_vector(x) = Vector{CNum}(_to_cnum.(x))

function tridiagonal_fixture(diagonal, coupling)
    n = length(diagonal)
    a = fill(_CNUM_ZERO, n, n)
    for i in 1:n
        a[i, i] = _to_cnum(diagonal[i])
        i < n && (a[i, i + 1] = a[i + 1, i] = _to_cnum(coupling))
    end
    return a
end

function block_diagonal(blocks::Matrix{CNum}...)
    n = sum(size(block, 1) for block in blocks)
    result = fill(_CNUM_ZERO, n, n)
    first = 1
    for block in blocks
        last = first + size(block, 1) - 1
        result[first:last, first:last] = block
        first = last + 1
    end
    return result
end

function warm_metrics(f; samples = 7)
    f()
    times = Float64[]
    allocations = Int[]
    for _ in 1:samples
        push!(times, @elapsed f())
        push!(allocations, @allocated f())
    end
    return median(times), median(allocations)
end

function fixture_metrics(name, a, forcing)
    solve = () -> solve_adjugate(a, forcing)
    x, fl, numerator = solve()
    residual, certified = residual_certificate(a, numerator, fl.determinant, forcing)
    elapsed, allocated = warm_metrics(solve)
    certificate_elapsed, certificate_allocated = warm_metrics(
        () -> residual_certificate(a, numerator, fl.determinant, forcing),
    )
    println(
        "METRIC|", name,
        "|n=", size(a, 1),
        "|det_terms=", coefficient_terms(fl.determinant),
        "|det_nodes=", coefficient_nodes(fl.determinant),
        "|det_chars=", coefficient_size(fl.determinant),
        "|adj_chars=", matrix_size(fl.adjugate),
        "|solution_chars=", matrix_size(x),
        "|certificate=", certified,
        "|residual_chars=", matrix_size(residual),
        "|median_us=", round(elapsed * 1e6; digits = 1),
        "|median_alloc=", allocated,
        "|certificate_us=", round(certificate_elapsed * 1e6; digits = 1),
        "|certificate_alloc=", certificate_allocated,
    )
    return x, fl, numerator
end

@variables a b c d e f g w w1 w2 w3 w4 w5 w6 w7 w8

println("SECTION|symbolic-2x2")
a2 = cnum_matrix([a b; c d])
f2 = cnum_vector([e, f])
x2fl, fl2, numerator2fl = fixture_metrics("symbolic_2x2", a2, f2)
x2, det2, numerator2 = solve_2x2(a2, f2)
special_elapsed, special_allocated = warm_metrics(() -> solve_2x2(a2, f2))
@test x2fl == x2
@test fl2.determinant == det2
@test numerator2fl == numerator2
quotient_residual = multiply(a2, x2)
quotient_residual = CNum[quotient_residual[i] - f2[i] for i in eachindex(f2)]
println("CHECK|2x2_specialized_equal=true")
println("CHECK|division_free_certificate=true")
println("CHECK|post_quotient_structural_residual=", all(_iszero_cnum, quotient_residual))
println(
    "METRIC|specialized_symbolic_2x2|median_us=", round(special_elapsed * 1e6; digits = 1),
    "|median_alloc=", special_allocated,
)

println("SECTION|coupled-symbolic")
a4 = tridiagonal_fixture([w1, w2, w3, w4], g)
f4 = cnum_vector([1, 2, 3, 4])
x4, fl4, numerator4 = fixture_metrics("coupled_4x4", a4, f4)

a8 = tridiagonal_fixture([w1, w2, w3, w4, w5, w6, w7, w8], g)
f8 = cnum_vector(collect(1:8))
x8, fl8, numerator8 = fixture_metrics("coupled_8x8", a8, f8)

println("SECTION|block-diagonal")
block_a = block_diagonal(a4, a4)
block_f = vcat(f4, f4)
xblock, flblock, numeratorblock = fixture_metrics("whole_blockdiag_8x8", block_a, block_f)
block_solve = () -> (solve_adjugate(a4, f4), solve_adjugate(a4, f4))
block_elapsed, block_allocated = warm_metrics(block_solve)
block_parts = block_solve()
block_solution = vcat(block_parts[1][1], block_parts[2][1])
block_numerators = vcat(block_parts[1][3], block_parts[2][3])
block_determinants = vcat(
    fill(block_parts[1][2].determinant, 4),
    fill(block_parts[2][2].determinant, 4),
)
rationally_equal = all(eachindex(numeratorblock)) do i
    numeratorblock[i] * block_determinants[i] ==
        block_numerators[i] * flblock.determinant
end
println(
    "METRIC|split_blockdiag_2x4|solution_chars=", matrix_size(block_solution),
    "|median_us=", round(block_elapsed * 1e6; digits = 1),
    "|median_alloc=", block_allocated,
    "|same_solution=", block_solution == xblock,
    "|division_free_rational_equivalence=", rationally_equal,
)

println("SECTION|numeric-oracles")
numeric_cases = [
    ([3.0 1.0; -2.0 4.0], [2.0, -1.0]),
    ([4.0 1.0 0.0 0.5; 1.0 3.0 -0.25 0.0; 0.0 -0.25 2.0 0.75; 0.5 0.0 0.75 5.0], [1.0, 2.0, 3.0, 4.0]),
]
for (i, (an, fn)) in enumerate(numeric_cases)
    xn, _, _ = solve_adjugate(cnum_matrix(an), cnum_vector(fn))
    numeric_x = ComplexF64[getfield(value, :z) for value in xn]
    oracle = an \ fn
    println("ORACLE|numeric_", i, "|error=", norm(numeric_x - oracle))
    @test numeric_x ≈ oracle
end

println("SECTION|singular-semantics")
unique = solve_structural_diagonal(cnum_matrix([2 0; 0 3]), cnum_vector([4, 9]))
nonunique = solve_structural_diagonal(cnum_matrix([2 0; 0 0]), cnum_vector([4, 0]))
@test to_num.(unique) == [2, 3]
@test to_num.(nonunique) == [2, 0]
println("SINGULAR|unique=accepted|solution=", to_num.(unique))
println("SINGULAR|nonunique_disconnected_null=accepted|solution=", to_num.(nonunique))
for (name, aa, ff) in [
        ("inconsistent", [2 0; 0 0], [4, 1]),
        ("timed_resonance", [1 0; 0 0], [0, 1]),
    ]
    try
        solve_structural_diagonal(cnum_matrix(aa), cnum_vector(ff))
        println("SINGULAR|", name, "=unexpected_accept")
    catch error
        println("SINGULAR|", name, "=refused|", sprint(showerror, error))
    end
end
connected_singular = cnum_matrix([1 1; 2 2])
try
    solve_adjugate(connected_singular, cnum_vector([3, 6]))
catch error
    println("SINGULAR|connected_compatible=refused|", sprint(showerror, error))
end
try
    solve_structural_diagonal(connected_singular, cnum_vector([3, 6]))
catch error
    println("SINGULAR|connected_specialized=refused|", sprint(showerror, error))
end

println("SECTION|bareiss-cancellation")
composite = _to_cnum(a * d - b * c)
bareiss_like = composite * _to_cnum(g) / composite
println("BAREISS|expression=", to_num(bareiss_like))
println("BAREISS|cancels_to_g=", bareiss_like == _to_cnum(g))
println("BAREISS|chars=", coefficient_size(bareiss_like))

println("SECTION|inference")
@test @inferred(faddeev_leverrier(a2)) isa FLResult
@test @inferred(solve_adjugate(a2, f2)) isa Tuple{Vector{CNum}, FLResult, Vector{CNum}}
@test @inferred(solve_2x2(a2, f2)) isa Tuple{Vector{CNum}, CNum, Vector{CNum}}
println("INFERENCE|faddeev_leverrier=true|solve_adjugate=true|solve_2x2=true")

println("DONE")

