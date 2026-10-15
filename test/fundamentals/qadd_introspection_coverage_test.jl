using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables

@testset "filtered QAdd variable collection" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables x y

    H = (x + im * y) * a + y * a'
    buffer = Set{Any}()
    @test Symbolics.get_variables!(buffer, H, [x, y]) === buffer
    @test buffer == Set(Symbolics.unwrap.([x, y]))
end
