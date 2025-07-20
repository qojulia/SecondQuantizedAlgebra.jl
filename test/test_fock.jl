using SecondQuantizedAlgebra
using SymbolicUtils
using Test

@testset "fock" begin
    @testset "Basic Operator Creation and Properties" begin
        hf = FockSpace(:c)
        a = Destroy(hf, :a)
        ad = a'
        b = Destroy(hf, :b)

        # Hash uniqueness tests
        @test !isequal(hash(a), hash(ad))
        @test !isequal(hash(a), hash(b))

        # Basic operator relations
        @test isequal(a, ad')
    end

    @testset "Algebraic Operations" begin
        hf = FockSpace(:c)
        a = Destroy(hf, :a)

        # Addition and scalar multiplication
        @test isequal(simplify(a+a), 2*a)
        @test isequal(simplify(a/2 + 0.5*a), a)
    end

    @testset "Commutation Relations" begin
        hf = FockSpace(:c)
        a = Destroy(hf, :a)

        # Canonical commutation relations
        @test isequal(a*a', 1+a'*a)
        @test isequal(simplify(a*a' + 1), 2 + a'*a)

        # More complex commutator identities
        @test isequal(simplify(-1*(a'*a + 1)*a + a), -1*a'*a^2)
        @test isequal(simplify(a'*a*a - a*a'*a), -1*a)
    end

    @testset "Hamiltonian Evolution" begin
        hf = FockSpace(:c)
        a = Destroy(hf, :a)

        # Single mode evolution
        ωc = 0.1313
        H = ωc*a'*a
        da = simplify(1.0im*(H*a - a*H))
        @test isequal(da, (0.0-1.0im)*ωc*a)
    end

    @testset "Symbolic Substitution" begin
        hf = FockSpace(:c)
        a = Destroy(hf, :a)
        @syms x y

        @testset "Substitution with Numbers" begin
            @test iszero(substitute(x*a, Dict(x=>0)))
            @test isequal(substitute(x*a, Dict(x=>1)), a)
            @test iszero(substitute(x*a, Dict(a=>0)))
        end

        @testset "Substitution with Symbols" begin
            @test isequal(substitute(x*a, Dict(x=>y)), y*a)
            @test isequal(substitute(x*a, Dict(a=>y)), x*y)
            @test isequal(substitute(x*(a+a'), Dict(x => y)), y*(a + a'))
        end
    end
    @testset "HilbertsSpace" begin
        h1 = FockSpace(:c1)
        h2 = FockSpace(:c2)
        h3 = FockSpace(:c3)
        h4 = FockSpace(:c4)

        h12 = h1 ⊗ h2
        h23 = h2 ⊗ h3
        h34 = h3 ⊗ h4
        h123 = h12 ⊗ h3
        h1234 = h123 ⊗ h4
        h1234_ = h12 ⊗ h34
        h = h1 ⊗ h2 ⊗ h3 ⊗ h4

        @test h123 == h1 ⊗ h2 ⊗ h3 == h1 ⊗ h23
        @test h == h1234 == h1234_
        @test ⊗(h1) == h1
        @test tensor(h1, h2, h3, h4) == h
        @test isless(h1, h2)
        @test !isless(h3, h2)
        @test isless(h12, h23)
        @test copy(h1) == h1
    end
end
