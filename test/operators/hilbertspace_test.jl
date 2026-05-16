using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: HilbertSpace
using Test

@testset "hilbert spaces" begin
    @testset "ProductSpace" begin
        h1 = FockSpace(:a)
        h2 = FockSpace(:b)
        h3 = FockSpace(:c)
        h4 = FockSpace(:d)

        h12 = h1 ⊗ h2
        @test h12 isa ProductSpace
        @test h12.spaces == (h1, h2)

        # Associativity
        h123_a = (h1 ⊗ h2) ⊗ h3
        h123_b = h1 ⊗ (h2 ⊗ h3)
        h123_c = h1 ⊗ h2 ⊗ h3
        @test h123_a == h123_b == h123_c
        @test h123_a.spaces == (h1, h2, h3)

        # 4 spaces
        h1234 = h1 ⊗ h2 ⊗ h3 ⊗ h4
        @test h1234.spaces == (h1, h2, h3, h4)
        @test (h1 ⊗ h2) ⊗ (h3 ⊗ h4) == h1234

        # tensor alias
        @test tensor(h1, h2, h3, h4) == h1234

        # isless
        @test isless(h1, h2)
        @test !isless(h3, h2)
    end

    @testset "length" begin
        # Single Hilbert spaces report 1
        @test length(FockSpace(:a)) == 1
        @test length(NLevelSpace(:atom, 3)) == 1
        @test length(PauliSpace(:p)) == 1
        @test length(SpinSpace(:s)) == 1
        @test length(PhaseSpace(:q)) == 1

        # ProductSpace reports number of factor spaces
        @test length(FockSpace(:a) ⊗ FockSpace(:b)) == 2
        @test length(FockSpace(:a) ⊗ NLevelSpace(:atom, 2) ⊗ PauliSpace(:p)) == 3
    end

end
