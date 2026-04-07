using SecondQuantizedAlgebra
using Test

# JET.report_package is infeasible: methods taking Symbolics types (BasicSymbolic,
# Complex{Num}) send JET deep into Symbolics.jl internals, causing unbounded analysis
# time. Even max_methods=2 doesn't help. Use @test_call on concrete entry points
# instead — JET then sees concrete types (e.g. Destroy, not abstract QSym) and
# avoids the combinatorial subtype × Symbolics blowup.
# @static if isempty(VERSION.prerelease)
#     @testset "Code linting" begin
#         using JET
#         JET.test_package(SecondQuantizedAlgebra;
#             target_modules = (SecondQuantizedAlgebra,))
#     end
# end
