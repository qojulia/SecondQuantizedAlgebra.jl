using SecondQuantizedAlgebra
using Test

@testset "@qnumbers" begin
    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
        @test acts_on(a) == [1]
        @test acts_on(b) == [2]
        @test is_destroy(a)
        @test is_destroy(b)
    end

end
