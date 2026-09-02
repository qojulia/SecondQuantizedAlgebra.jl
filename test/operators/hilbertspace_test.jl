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
        @test repr(h12) == "ℋ(a) ⊗ ℋ(b)"

        # Associativity
        h123_a = (h1 ⊗ h2) ⊗ h3
        h123_b = h1 ⊗ (h2 ⊗ h3)
        h123_c = h1 ⊗ h2 ⊗ h3
        @test h123_a == h123_b == h123_c
        @test repr(h123_a) == "ℋ(a) ⊗ ℋ(b) ⊗ ℋ(c)"

        # 4 spaces
        h1234 = h1 ⊗ h2 ⊗ h3 ⊗ h4
        @test repr(h1234) == "ℋ(a) ⊗ ℋ(b) ⊗ ℋ(c) ⊗ ℋ(d)"
        @test (h1 ⊗ h2) ⊗ (h3 ⊗ h4) == h1234

        # tensor alias
        @test tensor(h1, h2, h3, h4) == h1234

        # isless
        @test isless(h1, h2)
        @test !isless(h3, h2)
    end

    @testset "named space equality" begin
        @test PhaseSpace(:q) == PhaseSpace(:q)
        @test PhaseSpace(:q) != PhaseSpace(:r)
        @test SpinSpace(:s) == SpinSpace(:s)
        @test SpinSpace(:s) != SpinSpace(:t)
        @test hash(PhaseSpace(:q)) == hash(PhaseSpace(:q))
        @test hash(SpinSpace(:s)) == hash(SpinSpace(:s))
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

    @testset "hash + equality (Dict round-trip)" begin
        spaces = HilbertSpace[
            FockSpace(:f),
            NLevelSpace(:atom, 3),
            PauliSpace(:p),
            SpinSpace(:s),
            PhaseSpace(:osc),
        ]
        d = Dict{HilbertSpace, Int}()
        for (i, h) in enumerate(spaces)
            d[h] = i
        end
        for (i, h) in enumerate(spaces)
            @test d[h] == i
        end

        # ProductSpace equality and hashing
        p1 = FockSpace(:a) ⊗ NLevelSpace(:atom, 2)
        p2 = FockSpace(:a) ⊗ NLevelSpace(:atom, 2)
        p3 = FockSpace(:a) ⊗ NLevelSpace(:b, 2)
        @test p1 == p2
        @test hash(p1) == hash(p2)
        @test p1 != p3

        # ProductSpace ordering
        @test isless(p1, p3) || isless(p3, p1)
    end

end
