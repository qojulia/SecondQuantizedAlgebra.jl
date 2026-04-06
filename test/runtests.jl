using SecondQuantizedAlgebra
using ParallelTestRunner: ParallelTestRunner

# args = copy(ARGS)
# if !haskey(ENV, "CI") && !("--quickfail" in args)
#     push!(args, "--quickfail")
# end
# ParallelTestRunner.runtests(SecondQuantizedAlgebra, args)
ParallelTestRunner.runtests(SecondQuantizedAlgebra, ARGS)
