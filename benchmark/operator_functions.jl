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
    cA2 = cA^2

    group = SUITE["Formal operator functions"]
    group["Construction"]["Unary cos construction"] = @benchmarkable cos($A) seconds = 3 evals = 1
    group["Construction"]["64 unary function nodes"] = @benchmarkable qexpr_unary_batch($A, 64) seconds = 3 evals = 1
    group["Construction"]["Nested unary depth 8"] = @benchmarkable qexpr_nested($A, 8) seconds = 3 evals = 1
    group["Construction"]["Mixed formal expression"] = @benchmarkable qexpr_mixed_workflow($a, $A) seconds = 3 evals = 1
    group["Arithmetic"]["Collect 32 identical formal terms"] = @benchmarkable qexpr_collected_sum($cA, 32) seconds = 3 evals = 1
    group["Arithmetic"]["Mixed polynomial and formal sum 32"] = @benchmarkable qexpr_mixed_polynomial_sum($a, $cA, 32) seconds = 3 evals = 1
    group["Arithmetic"]["Formal power 16"] = @benchmarkable ($cA)^16 seconds = 3 evals = 1
    group["Taylor expansion"]["cos Taylor order 8"] = @benchmarkable taylor($cA, 0:8) seconds = 3 evals = 1
    group["Taylor expansion"]["expim Taylor order 8"] = @benchmarkable taylor($eA, 0:8) seconds = 3 evals = 1
    group["Taylor expansion"]["Repeated-node Taylor order 4"] = @benchmarkable taylor($cA2, 0:4) seconds = 3 evals = 1
    return SUITE
end
