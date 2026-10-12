# Superoperator construction and composition: SecondQuantizedAlgebra against the
# `Liouvillian` of FloquetExpansions.jl, which these types replace.
#
# FloquetExpansions is not registered, so this script is not part of `runbenchmarks.jl`.
# Run it from an environment that develops both packages:
#
#     julia --project=<env> benchmark/superoperators_vs_floquet.jl

using BenchmarkTools
using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: SuperOperator, liouvillian
import FloquetExpansions as FE
using Symbolics: @variables

function chain(M)
    h = ProductSpace(ntuple(k -> FockSpace(Symbol(:c, k)), M))
    a = [Destroy(h, Symbol(:a, k), k) for k in 1:M]
    @variables Δ κ J
    H = sum([Δ * a[k]' * a[k] for k in 1:M]) +
        sum([J * (a[k]' * a[k + 1] + a[k + 1]' * a[k]) for k in 1:(M - 1)])
    return H, a, κ
end

sqa_build(H, a, κ) = liouvillian(H, a; rates = fill(κ, length(a)))
sqa_lowered(H, a, κ) = SuperOperator(sqa_build(H, a, κ))

function fe_build(H, a, κ)
    L = FE.hamiltonian_action(H)
    for x in a
        L = L + κ * FE.dissipator(x)
    end
    return L
end

ms(t) = round(BenchmarkTools.time(t) / 1.0e6; digits = 3)

function main()
    println(rpad("modes", 7), rpad("operation", 34), rpad("SQA [ms]", 12), "FE [ms]")
    for M in (2, 4, 8)
        H, a, κ = chain(M)
        Ls = sqa_build(H, a, κ)
        Ss = SuperOperator(Ls)
        Lf = fe_build(H, a, κ)
        rows = (
            ("build (structured)", @benchmark(sqa_build($H, $a, $κ)), @benchmark(fe_build($H, $a, $κ))),
            ("build + lower to sandwiches", @benchmark(sqa_lowered($H, $a, $κ)), @benchmark(fe_build($H, $a, $κ))),
            ("compose L ∘ L", @benchmark($Ss ∘ $Ss), @benchmark(FE.compose($Lf, $Lf))),
        )
        for (name, ts, tf) in rows
            println(rpad(M, 7), rpad(name, 34), rpad(ms(minimum(ts)), 12), ms(minimum(tf)))
        end
    end
    return nothing
end

main()
