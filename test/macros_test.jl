using SecondQuantizedAlgebra
using Test

@testset "@qnumbers" begin
    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
        @test a.space_index == 1
        @test b.space_index == 2
    end

end
