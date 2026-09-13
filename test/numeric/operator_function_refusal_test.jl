using SecondQuantizedAlgebra
import QuantumOpticsBase as QOB
using Test

@testset "Formal operator numeric refusal — QuantumOptics" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    expr = cos(a + a')
    b = QOB.FockBasis(5)
    ψ = QOB.fockstate(b, 0)
    d = Dict{SecondQuantizedAlgebra.QSym, Any}()

    for thunk in (
            () -> to_numeric(expr, b),
            () -> to_numeric(expr, b, d),
            () -> to_numeric(expr, ψ),
            () -> numeric_average(expr, ψ),
            () -> SecondQuantizedAlgebra.expect(expr, ψ),
        )
        err = try
            thunk()
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("taylor(expr, 0:n)", sprint(showerror, err))
    end

    @test to_numeric(taylor(expr, 0:2), b) !== nothing
end
