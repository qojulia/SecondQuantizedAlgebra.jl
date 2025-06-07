using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

@testset "spin" begin
    @testset "Pauli Operators" begin
        @testset "Single Pauli Space" begin
            hs1 = PauliSpace(:Spin1)
            s(axis) = Pauli(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]

            # Basic Pauli algebra
            @test isequal(s(1)*s(2), 1im*s(3))
            @test !isequal(s(1)*s(2), 1im*s(2))
            @test isequal(s(1)*s(3), -1im*s(2))
            @test isequal(s(3)*s(3), 1)
            @test isequal(s(3)*s(1), 1im*s(2))
            @test isequal(s(1)*s(2)*s(3), 1im)
        end

        @testset "Multi-Pauli Space Operations" begin
            hs1 = PauliSpace(:Spin1)
            hs2 = PauliSpace(:Spin2)
            h = hs1 ⊗ hs2

            σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)
            σx(i) = σ(i, 1)
            σy(i) = σ(i, 2)
            σz(i) = σ(i, 3)

            sx(i) = σ(i, :x)
            sy(i) = σ(i, :y)
            sz(i) = σ(i, :z)

            @test sx(1) == σx(1)

            # Commutation relations between different spins
            @test isequal(σx(2)*σx(1), σx(1)*σx(2))
            @test isequal(σy(2)*σx(1), σx(1)*σy(2))
            @test isequal(σz(2)*σz(1), σz(1)*σz(2))
        end

        @testset "Pauli Hamiltonian Construction" begin
            hs1 = PauliSpace(:Spin1)
            hs2 = PauliSpace(:Spin2)
            h = hs1 ⊗ hs2

            σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)
            σx(i) = σ(i, 1)
            σz(i) = σ(i, 3)

            @cnumbers J
            Δi(i) = cnumber(Symbol(:Δ_, i))
            H = Δi(1)*σz(1) + Δi(2)*σz(2) + J*σx(1)*σx(2)

            # Test that Hamiltonian construction works
            @test H isa SecondQuantizedAlgebra.QAdd
        end
    end

    @testset "Collective Spin Operators" begin
        @testset "Single Spin Space" begin
            hcs1 = SpinSpace(:Spin1)
            S(axis) = Spin(hcs1, :S, axis) # axis ∈ [1,2,3] → [x,y,z]

            # Axis equivalences
            @test S(1) == S(:x)
            @test S(2) == S(:Y)
            @test S(:z) == S(:Z)
            @test S(1) ≠ S(2)

            # Basic commutation relations
            @test isequal(simplify(S(:x)*S(:y) - S(:y)*S(:x)), 1im*S(:z))
            @test isequal(simplify(S(:x)*S(:z) - S(:z)*S(:x)), -1im*S(:y))
            @test isequal(simplify(S(:y)*S(:z) - S(:z)*S(:y)), 1im*S(:1))
        end

        @testset "Spin Operator Properties" begin
            hcs1 = SpinSpace(:Spin1)
            S(axis) = Spin(hcs1, :S, axis)

            # Algebraic properties
            @test isequal(S(:x)*S(:y), 1*S(:x)*S(:y))
            @test !isequal(S(:x)*S(:y), S(:x)*S(:z))
            @test isequal(S(1)*S(1), (S(1))^2)
            @test isequal(average(S(1) + S(2)), average(S(2) + S(1)))
            @test isequal(simplify(S(2) + S(2)), 2S(2))

            # Scalar multiplication commutativity
            @test isequal(2*S(1)*S(2)*S(3), S(1)*S(2)*S(3)*2)
        end

        @testset "Multi-Spin Operations" begin
            hcs1 = SpinSpace(:Spin1)
            hcs2 = SpinSpace(:Spin2)
            h = hcs1 ⊗ hcs2

            S(i, axis) = Spin(h, Symbol(:S_, i), axis, i)
            Sx(i) = S(i, 1)
            Sy(i) = S(i, 2)
            Sz(i) = S(i, 3)
            Sm(i) = Sx(i) - 1im*Sy(i)
            Sp(i) = Sx(i) + 1im*Sy(i)

            # Commutation between different spins
            @test isequal(Sx(2)*Sx(1), Sx(1)*Sx(2))
            @test isequal(Sy(2)*Sx(1), Sx(1)*Sy(2))
            @test isequal(Sz(2)*Sz(1), Sz(1)*Sz(2))

            # Note: Commented test case for potential future implementation
            # @test isequal( simplify(Sy(1)Sz(2)Sx(1)), simplify(Sx(1)Sy(1)Sz(2) - 1im*Sz(1)*Sz(2)) )
        end
    end
end
