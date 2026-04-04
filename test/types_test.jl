using SecondQuantizedAlgebra
using Test

@testset "type hierarchy" begin
    @test QSym <: QField
    @test QTerm <: QField
    @test !(QSym <: QTerm)
    @test !(QTerm <: QSym)
    @test isabstracttype(QField)
    @test isabstracttype(QSym)
    @test isabstracttype(QTerm)
end
