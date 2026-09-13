function qexpr_unary_batch(A, n)
    out = Vector{SecondQuantizedAlgebra.QExpr}(undef, n)
    @inbounds for i in eachindex(out)
        out[i] = isodd(i) ? cos(A) : sin(A)
    end
    return out
end

function qexpr_nested(A, n)
    result = cos(A)
    for i in 2:n
        result = isodd(i) ? cos(result) : sin(result)
    end
    return result
end

function qexpr_mixed_workflow(a, A)
    cA = cos(A)
    sA = sin(A)
    eA = expim(A)
    sumexpr = 2 * cA + 3 * sA + eA
    product = a * cA * A * sA * a'
    return sumexpr, product
end

function qexpr_collected_sum(cA, n)
    result = zero(cA)
    for i in 1:n
        result += i * cA
    end
    return result
end

function qexpr_mixed_polynomial_sum(a, cA, n)
    result = zero(cA)
    for i in 1:n
        result += i * a + i * cA
    end
    return result
end

function benchmark_operator_functions!(SUITE)
    h = FockSpace(:qexpr)
    a = Destroy(h, :a)
    A = a + a'
    cA = cos(A)
    eA = expim(A)
    cA4 = cA^4

    group = SUITE["Formal operator functions"]
    group["Unary cos construction"] = @benchmarkable cos($A) seconds = 3 evals = 1
    group["64 unary function nodes"] = @benchmarkable qexpr_unary_batch($A, 64) seconds = 3 evals = 1
    group["Nested unary depth 8"] = @benchmarkable qexpr_nested($A, 8) seconds = 3 evals = 1
    group["Mixed formal expression"] = @benchmarkable qexpr_mixed_workflow($a, $A) seconds = 3 evals = 1
    group["Collect 32 identical formal terms"] = @benchmarkable qexpr_collected_sum($cA, 32) seconds = 3 evals = 1
    group["Mixed polynomial/formal sum 32"] = @benchmarkable qexpr_mixed_polynomial_sum($a, $cA, 32) seconds = 3 evals = 1
    group["Formal power 16"] = @benchmarkable ($cA)^16 seconds = 3 evals = 1
    group["cos Taylor order 8"] = @benchmarkable taylor($cA, 0:8) seconds = 3 evals = 1
    group["expim Taylor order 8"] = @benchmarkable taylor($eA, 0:8) seconds = 3 evals = 1
    group["Repeated-node Taylor order 8"] = @benchmarkable taylor($cA4, 0:8) seconds = 3 evals = 1
    return SUITE
end
