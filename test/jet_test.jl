using SecondQuantizedAlgebra
using Test

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        JET.test_package(SecondQuantizedAlgebra; target_modules=(SecondQuantizedAlgebra,))
    end
end
