function nested_commutator(term, n)
    if n == 0
        return term
    else
        return commutator(term, nested_commutator(term, n - 1))
    end
end

function benchmark_commutator!(SUITE)
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)

    term = a'*a
    n=5

    SUITE["Micro Benchmarks"]["Nested commutator"] = @benchmarkable nested_commutator(
        $term, $n
    ) seconds = 20
    return nothing
end
