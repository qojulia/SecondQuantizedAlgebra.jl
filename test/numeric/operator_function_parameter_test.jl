using SecondQuantizedAlgebra: FockSpace, Destroy, QuantumOpticsBackend, @qnumbers, @variables, to_numeric
import QuantumOpticsBase as QOB
using Test: @test, @testset, @test_throws

@testset "Formal operator functions require static parameters" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    @variables x y
    b = QOB.FockBasis(3)
    A = a + a'

    for expr in (cos(x * A), x * cos(A), cos(sin(x * A)), cos(A) + x * A)
        @test_throws "supply values for x via `parameter`" to_numeric(
            expr, h, 3; backend = QuantumOpticsBackend(),
        )
    end

    @test_throws "supply values for y via `parameter`" to_numeric(
        y * cos(x * A), h, 3;
        backend = QuantumOpticsBackend(), parameter = Dict(x => 0),
    )
    @test_throws "supply values for y via `parameter`" to_numeric(
        cos(x * A), h, 3;
        backend = QuantumOpticsBackend(), parameter = Dict(x => y),
    )
    @test_throws "nonempty time_parameter is unsupported" to_numeric(
        cos(x * A), h, 3;
        backend = QuantumOpticsBackend(), time_parameter = Dict(x => identity),
    )

    result = to_numeric(
        y * cos(x * A), h, 3;
        backend = QuantumOpticsBackend(), parameter = Dict(x => 0, y => 2),
    )
    @test QOB.dense(result).data == 2 * QOB.dense(one(b)).data
end
