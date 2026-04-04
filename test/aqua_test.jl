using SecondQuantizedAlgebra
using Test

@testset "best practices" begin
    using Aqua
    Aqua.test_all(SecondQuantizedAlgebra; ambiguities=false, piracies=false)
end
