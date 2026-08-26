using SecondQuantizedAlgebra
using Symbolics
using BenchmarkTools
using Statistics
using QuantumOpticsBase
using Test

include(joinpath(@__DIR__, "materializer_prototype.jl"))
using .MaterializerPrototype
const MP = MaterializerPrototype
const SQA = SecondQuantizedAlgebra

@variables α::Number η::Real θ::Real r::Real ϕ::Real
@variables Δa::Real Δb::Real J::Real κ::Real t::Real ω::Real

h = FockSpace(:a) ⊗ FockSpace(:b)
a = Destroy(h, :a, 1)
b = Destroy(h, :b, 2)

atom = NLevelSpace(:atom, 3)
σ = Transition(atom, :σ, 1, 2)
W = Rational{Int}[
    1//3 2//3 2//3
    2//3 1//3 -2//3
    -2//3 2//3 -1//3
]

named = Dict{Symbol, Any}(
    :displacement => Displace(a, α),
    :fock_chain => Rotation(a, θ) * Squeeze(a, r, ϕ),
    :coupled => Rotation(a, b, θ) * Squeeze(a, b, r),
    :nlevel => Rotation(σ, W),
)
pairs = Dict(name => MP.coefficient_pair(U) for (name, U) in named)
generic = Dict(name => MP.materialize_direct(pairs[name], U.relations) for (name, U) in named)
wrapped = Dict(name => MP.materialize_via_wrapper(pairs[name], U.relations) for (name, U) in named)

println("CORRECTNESS")
for name in sort!(collect(keys(named)))
    U = named[name]
    G = generic[name]
    V = wrapped[name]
    println(name, " direct=", MP.same_transform(U, G), " wrapped=", MP.same_transform(U, V))
end

H = Δa * a' * a + Δb * b' * b + J * (a' * b + b' * a) +
    η * (a + a') + κ * (a' * b' + b * a)
Hlevel = sum((i + 2j) * Transition(atom, :σ, i, j) for i in 1:3, j in 1:3)
observables = Dict(
    :displacement => H,
    :fock_chain => H,
    :coupled => H,
    :nlevel => Hlevel,
)
for name in sort!(collect(keys(named)))
    U = named[name]
    G = generic[name]
    expr = observables[name]
    transformed_equal = SQA.transform(expr, U) == SQA.transform(expr, G)
    roundtrip = simplify(SQA.conjugate(SQA.conjugate(expr, G), inv(G))) == simplify(expr)
    println(name, " application=", transformed_equal, " roundtrip=", roundtrip)
end

# The scalar lift is passed separately, consumed into QAdd, and never retained.
operator_gauge = SQA._rule_qadd((SQA._to_cnum(ω), SQA.Op[a', a]))
scalar_gauge = SQA._to_cnum(η)
timed = MP.materialize_direct(
    pairs[:displacement], operator_gauge, scalar_gauge, t,
)
println("scalar-lift type=", typeof(scalar_gauge),
    " stored-only-as-qadd=", timed.gauge isa SQA.QAdd,
    " scalar-present=", any(isempty(term.ops) for (term, _) in timed.gauge))

println("INFERENCE")
inferred_direct = @inferred MP.materialize_direct(
    pairs[:coupled], named[:coupled].relations,
)
inferred_rules = @inferred MP._materialize_rules(pairs[:coupled].forward)
inferred_nested = @inferred MP.nested_conjugate(H, MP.nested(generic[:coupled]))
println("at_inferred=", (
    typeof(inferred_direct), typeof(inferred_rules), typeof(inferred_nested),
))
println("direct=", Base.return_types(MP.materialize_direct,
    Tuple{MP.AffineRulePair, Vector{SQA.ParamRelation}}))
println("rules=", Base.return_types(MP._materialize_rules,
    Tuple{MP.CoefficientRuleMap}))
println("conjugate=", Base.return_types(SQA.conjugate,
    Tuple{SQA.QAdd, SQA.UnitaryTransform{SQA.StaticTime}}))

function stats(f; samples = 21)
    trial = @benchmark $f() samples=21 evals=1 seconds=10
    estimate = median(trial)
    return (time_ns = estimate.time, bytes = estimate.memory, allocs = estimate.allocs)
end

builders = Dict{Symbol, Function}(
    :displacement => () -> Displace(a, α),
    :fock_chain => () -> Rotation(a, θ) * Squeeze(a, r, ϕ),
    :coupled => () -> Rotation(a, b, θ) * Squeeze(a, b, r),
    :nlevel => () -> Rotation(σ, W),
)

println("CONSTRUCTION_BENCHMARKS")
for name in sort!(collect(keys(named)))
    pair = pairs[name]
    relations = named[name].relations
    println(name, " named=", stats(builders[name]),
        " direct=", stats(() -> MP.materialize_direct(pair, relations)),
        " wrapper=", stats(() -> MP.materialize_via_wrapper(pair, relations)))
end

println("HOT_WORKFLOW_BENCHMARKS")
for name in sort!(collect(keys(named)))
    U = generic[name]
    N = MP.nested(U)
    expr = observables[name]
    second = name === :nlevel ? generic[:nlevel] : generic[:coupled]
    Nsecond = MP.nested(second)
    println(name,
        " conjugate_direct=", stats(() -> SQA.conjugate(expr, U)),
        " conjugate_nested=", stats(() -> MP.nested_conjugate(expr, N)),
        " inv_direct=", stats(() -> inv(U)),
        " inv_nested=", stats(() -> MP.nested_inv(N)),
        " display_direct=", stats(() -> sprint(show, U)),
        " display_nested=", stats(() -> MP.nested_render(N)))
    if name === :coupled
        println(name,
            " compose_direct=", stats(() -> U * second),
            " compose_nested=", stats(() -> MP.nested_compose(N, Nsecond)))
    end
end

println("EXPRESSION_SIZES")
for name in sort!(collect(keys(named)))
    U = generic[name]
    terms = sum(length(image) for image in values(U.rules))
    invterms = sum(length(image) for image in values(U.inverse_rules))
    println(name, " generators=", length(U.generators),
        " forward_terms=", terms, " inverse_terms=", invterms,
        " summary_chars=", ncodeunits(sprint(show, U)))
end

numeric_substitutions = Dict(
    Δa => 1.1, Δb => 0.9, J => 0.2, η => 0.1, κ => 0.03,
    θ => 0.4, r => 0.2, ϕ => 0.1, α => 0.15,
)
numeric_basis = FockBasis(3) ⊗ FockBasis(3)
coupled_direct = generic[:coupled]
coupled_nested = MP.nested(coupled_direct)
numeric_direct() = SQA.to_numeric(
    substitute(SQA.conjugate(H, coupled_direct), numeric_substitutions), numeric_basis,
)
numeric_nested() = SQA.to_numeric(
    substitute(MP.nested_conjugate(H, coupled_nested), numeric_substitutions), numeric_basis,
)
println("NUMERIC_CONVERSION")
println("equal=", numeric_direct() == numeric_nested(),
    " direct=", stats(numeric_direct), " nested=", stats(numeric_nested))

nothing

