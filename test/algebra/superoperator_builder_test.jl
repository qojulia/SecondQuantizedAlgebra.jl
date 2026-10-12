using SecondQuantizedAlgebra
using Symbolics: @variables
using Test

import SecondQuantizedAlgebra: SuperOperator, SuperBuilder, LindbladianBuilder, QAdd,
    accumulate!, build, sprepost, spre, spost
import MutableArithmetics as MA

@testset "Superoperator builders" begin
    h = FockSpace(:c)
    @qnumbers a::Destroy(h)
    @variables g

    maps = [k * sprepost(a^k, a') + g * spre(a' * a) for k in 1:6]
    fold = foldl(+, maps)

    b = SuperBuilder()
    for S in maps
        accumulate!(b, S)
    end
    @test build(b) == fold
    @test sum(maps) == fold
    @test MA.operate(sum, maps) == fold

    b2 = SuperBuilder()
    for k in 1:6
        accumulate!(b2, k, a^k, a')
        accumulate!(b2, g, a' * a, 1)
    end
    @test build(b2) == fold

    Ls = [k * dissipator(a^k) + hamiltonian_action(g * a' * a) for k in 1:4]
    Lfold = foldl(+, Ls)
    lb = LindbladianBuilder{QAdd}()
    for L in Ls
        accumulate!(lb, L)
    end
    @test build(lb) == Lfold
    @test sum(Ls) == Lfold

    sb = SuperBuilder()
    for L in Ls
        accumulate!(sb, L)
    end
    @test build(sb) == SuperOperator(Lfold)
end
