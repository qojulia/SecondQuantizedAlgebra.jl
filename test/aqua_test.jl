using SecondQuantizedAlgebra
using QuantumOpticsBase: expect
using Test

@testset "best practices" begin
    using Aqua
    # TODO: remove `treat_as_own = [expect]` once the piracies in src/numeric.jl
    # are resolved — `expect(::Num, ...)` and `expect(::BasicSymbolic, ...)` are
    # methods on a QuantumOpticsBase function with externally-owned argument types.
    Aqua.test_all(SecondQuantizedAlgebra; piracies = (treat_as_own = [expect],))
end
