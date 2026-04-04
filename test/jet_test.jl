using SecondQuantizedAlgebra
using Test

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        rep = report_package(SecondQuantizedAlgebra)
        @show rep
        @test_broken length(JET.get_reports(rep)) == 0
    end
end
