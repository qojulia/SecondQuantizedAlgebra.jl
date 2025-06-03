using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "average" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

@test isequal(average(2a), 2*average(a))

@test isequal(average(2*(a+a)), 2*(average(a) + average(a)))
@test isequal(average(a^2), average(a*a))

@test isequal(simplify(average(σ)+average(σ)), average(2σ))
@test iszero(simplify(average(σ)-average(σ)))

ωc, ωa = cnumbers("ω_c ω_a")
@test isequal(average(ωc),ωc)
@test isequal(average(ωc*a),ωc*average(a))
@test isequal(average(ωc*(a+a')) , ωc*average(a) + ωc*average(a'))

@test iszero(average(a*σ)*average(a) - average(a)*average(a*σ))

end # testset
