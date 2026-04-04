using SecondQuantizedAlgebra
using Test

@testset "@qnumbers" begin
    @testset "Single space" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h)
        @test a isa Destroy
        @test a.name == :a
        @test a.space_index == 1
    end

    @testset "Multiple operators" begin
        h = FockSpace(:c)
        @qnumbers a::Destroy(h) b::Destroy(h)
        @test a.name == :a
        @test b.name == :b
    end

    @testset "Product space" begin
        h = FockSpace(:a) ⊗ FockSpace(:b)
        @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
        @test a.space_index == 1
        @test b.space_index == 2
    end

    @testset "Create operators" begin
        h = FockSpace(:c)
        @qnumbers ad::Create(h)
        @test ad isa Create
        @test ad.name == :ad
    end
end
