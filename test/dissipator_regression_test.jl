using Test
using SecondQuantizedAlgebra
using Symbolics: @variables

@testset "PR #99 dissipator regression" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    a = Destroy(h, :a, 1)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)

    @variables N Δ κ γ R ν

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    Jd = adjoint.(J)
    rates = [κ, γ, R, ν]

    function decay(op)
        return simplify(
            expand_completeness(
                0.5 * rates[1] * (2 * Jd[1] * op * J[1] - Jd[1] * J[1] * op - op * Jd[1] * J[1]) +
                    sum(
                    Σ(0.5 * rates[x] * (2 * Jd[x] * op * J[x] - Jd[x] * J[x] * op - op * Jd[x] * J[x]), i)
                        for x in 2:length(J)
                )
            )
        )
    end

    @test isequal(decay(a' * a), simplify(-κ * (a' * a)))
    @test isequal(decay(σ(2, 2, j)), simplify(R - (R + γ) * σ(2, 2, j)))
end
