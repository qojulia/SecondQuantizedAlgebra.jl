using SecondQuantizedAlgebra: FockSpace, Destroy, NLevelSpace, Transition, Index,
    IndexedOperator, Σ, QuantumOpticsBackend, QuantumToolboxBackend, @variables, to_numeric
import QuantumOpticsBase as QOB
import QuantumToolbox as QTB
import SciMLOperators
using Test: @test, @testset, @test_throws

@testset "Time-dependent coefficients outside formal functions" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    A = a + a'
    @variables x y z

    for backend in (QuantumOpticsBackend(), QuantumToolboxBackend())
        @testset "$(typeof(backend))" begin
            at = backend isa QuantumOpticsBackend ? ((op, t) -> QOB.dense(QOB.static_operator(op(t)))) : ((op, t) -> op(nothing, t))
            expressions = (
                x * cos(A),
                x * cos(z * A),
                x * cos(A) + y * sin(A),
                (cos(A) + x * a) * (sin(A) + y * a'),
                x * cos(sin(A)),
                cos(A),
            )
            for expr in expressions
                timed = to_numeric(
                    expr, h, 3; backend, parameter = Dict(z => 0.3),
                    time_parameter = Dict(x => t -> 0.5 * cos(t), y => t -> 1 + t),
                )
                @test backend isa QuantumOpticsBackend ? timed isa QOB.TimeDependentSum : timed isa QTB.QuantumObjectEvolution
                for t in (0.0, 0.7, 2.0)
                    expected = to_numeric(
                        expr, h, 3; backend,
                        parameter = Dict(x => 0.5 * cos(t), y => 1 + t, z => 0.3),
                    )
                    @test Matrix(at(timed, t).data) ≈ Matrix(expected.data)
                end
            end

            for expr in (cos(x * A), y * cos(x * A), cos(x * sin(A)))
                @test_throws ArgumentError to_numeric(
                    expr, h, 3; backend, time_parameter = Dict(x => identity, y => identity),
                )
            end
            @test_throws ArgumentError to_numeric(x * cos(A) + y * a, h, 3; backend, time_parameter = Dict(x => identity))
            @test_throws ArgumentError to_numeric(x * cos(A), h, 3; backend, time_parameter = Dict(x => identity), op_type = Matrix)
            lazy = to_numeric(x * cos(A), h, 3; backend, time_parameter = Dict(x => 2), op_type = identity)
            @test Matrix(at(lazy, 0.4).data) ≈ Matrix(to_numeric(2 * cos(A), h, 3; backend).data)
        end
    end

    @testset "Indexed static argument" begin
        atoms = NLevelSpace(:atom, 2)
        i = Index(atoms, :i, 2, atoms)
        S = Σ(IndexedOperator(Transition(atoms, :σ, 2, 2), i), i)
        b = QOB.tensor(QOB.NLevelBasis(2), QOB.NLevelBasis(2))
        timed = to_numeric(x * cos(S), b; time_parameter = Dict(x => t -> 0.5 * cos(t)))
        expected = cos(Matrix(to_numeric(S, b).data))
        @test Matrix(QOB.dense(QOB.static_operator(timed(0.7))).data) ≈ (0.5 * cos(0.7)) * expected
    end

    @testset "QuantumToolbox solver parameters" begin
        backend = QuantumToolboxBackend()
        timed = to_numeric(
            x * cos(A) + y * sin(A), h, 3; backend,
            time_parameter = Dict(x => (p, t) -> p[1] * cos(t), y => t -> sin(t)),
        )
        for p in ([0.2], [0.8]), t in (0.0, 0.7)
            expected = to_numeric(
                x * cos(A) + y * sin(A), h, 3; backend,
                parameter = Dict(x => p[1] * cos(t), y => sin(t)),
            )
            @test Matrix(timed(p, t).data) ≈ Matrix(expected.data)
        end
        @test_throws ArgumentError to_numeric(
            x * cos(A), h, 3; backend = QuantumOpticsBackend(),
            time_parameter = Dict(x => (p, t) -> p[1]),
        )
    end
end
