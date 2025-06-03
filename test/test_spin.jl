using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

@testset "spin" begin
    hs1 = PauliSpace(:Spin1)
    hs2 = PauliSpace(:Spin2)
    h = hs1 ⊗ hs2

    s(axis) = Pauli(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]
    @test isequal(s(1)*s(2), 1im*s(3))
    @test !isequal(s(1)*s(2), 1im*s(2))
    @test isequal(s(1)*s(3), -1im*s(2))
    @test isequal(s(3)*s(3), 1)
    @test isequal(s(3)*s(1), 1im*s(2))
    @test isequal(s(1)*s(2)*s(3), 1im)

    σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)
    σx(i) = σ(i, 1)
    σy(i) = σ(i, 2)
    σz(i) = σ(i, 3)

    sx(i) = σ(i, :x)
    sy(i) = σ(i, :y)
    sz(i) = σ(i, :z)
    sx(1) == σx(1)

    @test isequal(σx(2)*σx(1), σx(1)*σx(2))
    @test isequal(σy(2)*σx(1), σx(1)*σy(2))
    @test isequal(σz(2)*σz(1), σz(1)*σz(2))

    @cnumbers J
    Δi(i) = cnumber(Symbol(:Δ_, i))
    H = Δi(1)*σz(1) + Δi(2)*σz(2) + J*σx(1)*σx(2)

    ### collective spin ###
    hcs1 = SpinSpace(:Spin1)
    hcs2 = SpinSpace(:Spin2)
    h = hcs1 ⊗ hcs2

    S(axis) = Spin(hcs1, :S, axis) # axis ∈ [1,2,3] → [x,y,z]

    S(1)==S(:x)
    S(2)==S(:Y)
    S(:z)==S(:Z)
    S(1)≠S(2)

    @test isequal(simplify(S(:x)*S(:y) - S(:y)*S(:x)), 1im*S(:z))
    @test isequal(simplify(S(:x)*S(:z) - S(:z)*S(:x)), -1im*S(:y))
    @test isequal(simplify(S(:y)*S(:z) - S(:z)*S(:y)), 1im*S(:1))

    @test isequal(S(:x)*S(:y), 1*S(:x)*S(:y))
    @test !isequal(S(:x)*S(:y), S(:x)*S(:z))
    @test isequal(S(1)*S(1), (S(1))^2)
    @test isequal(average(S(1) + S(2)), average(S(2) + S(1)))
    @test isequal(simplify(S(2) + S(2)), 2S(2))

    # error
    @test isequal(2*S(1)*S(2)*S(3), S(1)*S(2)*S(3)*2)

    S(i, axis) = Spin(h, Symbol(:S_, i), axis, i)
    Sx(i) = S(i, 1)
    Sy(i) = S(i, 2)
    Sz(i) = S(i, 3)
    Sm(i) = Sx(i) - 1im*Sy(i)
    Sp(i) = Sx(i) + 1im*Sy(i)

    @test isequal(Sx(2)*Sx(1), Sx(1)*Sx(2))
    @test isequal(Sy(2)*Sx(1), Sx(1)*Sy(2))
    @test isequal(Sz(2)*Sz(1), Sz(1)*Sz(2))
    # @test isequal( simplify(Sy(1)Sz(2)Sx(1)), simplify(Sx(1)Sy(1)Sz(2) - 1im*Sz(1)*Sz(2)) )

end # testset
