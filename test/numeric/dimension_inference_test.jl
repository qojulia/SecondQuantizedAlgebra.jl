using SecondQuantizedAlgebra: NLevelSpace, PauliSpace, FockSpace, PhaseSpace, SpinSpace,
    Transition, Pauli, QuantumOpticsBackend, QuantumToolboxBackend, to_numeric, ⊗, @variables
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
using Test: @test, @testset, @test_throws, @inferred

@testset "Dimensions from the supplied Hilbert space" begin
    atom = NLevelSpace(:atom, 3)
    qubit = PauliSpace(:qubit)
    P = Transition(atom, :P, 2, 2)
    sx = Pauli(qubit, :σ, 1)
    h = atom ⊗ qubit
    Ph = Transition(h, :P, 2, 2, 1)
    @variables λ

    for backend in (QuantumOpticsBackend(), QuantumToolboxBackend())
        @testset "$(typeof(backend))" begin
            for (op, space, dims) in ((P, atom, 3), (sx, qubit, 2), (Ph, h, (3, 2)))
                actual = to_numeric(op, space; backend)
                expected = to_numeric(op, space, dims; backend)
                @test actual ≈ expected
            end
            @test size(to_numeric(Ph, h; backend)) == (6, 6)
            @test to_numeric(1, h; backend) ≈ to_numeric(1, h, (3, 2); backend)
            @test to_numeric(cos(Ph), h; backend) ≈ to_numeric(cos(Ph), h, (3, 2); backend)
            @test to_numeric(λ * Ph, h; backend, parameter = Dict(λ => 0.5)) ≈
                to_numeric(0.5 * Ph, h, (3, 2); backend)
            converted = to_numeric([P, cos(P)], atom; backend)
            expected = [to_numeric(P, atom, 3; backend), to_numeric(cos(P), atom, 3; backend)]
            @test all(isapprox.(converted, expected))
            for space in (FockSpace(:a), PhaseSpace(:q), SpinSpace(:s), atom ⊗ FockSpace(:b))
                @test_throws ArgumentError to_numeric(1, space; backend)
            end
        end
    end

    @test (@inferred to_numeric(P, atom; backend = QuantumOpticsBackend())) isa QOB.Operator
    @test (@inferred to_numeric(P, atom; backend = QuantumToolboxBackend())) isa QTB.QuantumObject
    @test_throws ArgumentError to_numeric(P, atom)

    @testset "QuantumToolbox parameter modes with and without dimensions" begin
        for dims in ((), (3,))
            dynamic = to_numeric(λ * cos(P), atom, dims...; backend = QuantumToolboxBackend())
            numeric = to_numeric(
                λ * cos(P), atom, dims...;
                backend = QuantumToolboxBackend(), default_time_dependence = false
            )
            explicit = to_numeric(
                λ * cos(P), atom, dims...;
                backend = QuantumToolboxBackend(), default_time_dependence = false,
                time_parameter = Dict(λ => (p, t) -> p[λ] * cos(t))
            )
            for t in (0.0, 0.7)
                expected = to_numeric(0.5 * cos(t) * cos(P), atom, 3; backend = QuantumToolboxBackend())
                @test dynamic(Dict(λ => t -> 0.5 * cos(t)), t) ≈ expected
                @test numeric(Dict(λ => 0.5 * cos(t)), t) ≈ expected
                @test explicit(Dict(λ => 0.5), t) ≈ expected
            end
        end
    end
end
