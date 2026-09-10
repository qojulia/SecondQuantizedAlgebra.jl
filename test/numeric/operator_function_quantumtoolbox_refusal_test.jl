using SecondQuantizedAlgebra
import QuantumToolbox as QTB
using Test

@testset "Formal operator numeric refusal — QuantumToolbox" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    expr = cos(a + a')
    d = Dict{SecondQuantizedAlgebra.QSym, Any}()

    for thunk in (
            () -> to_numeric(expr, 6),
            () -> to_numeric(expr, 6, d),
            () -> to_numeric(expr, (6,), d),
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

    @test to_numeric(taylor(expr, 0:2), 6) isa QTB.QuantumObject
end
