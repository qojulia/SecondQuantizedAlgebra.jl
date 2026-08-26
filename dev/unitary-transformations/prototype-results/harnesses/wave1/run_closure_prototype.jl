using Test
using SecondQuantizedAlgebra
include(joinpath(@__DIR__, "closure_prototype.jl"))
using .ClosurePrototype: ClosureLimits, discover_closure, infer_operator_context,
    storage_comparison, warm_measure

function record_fixture!(rows, name, space, generator, seeds; limits = ClosureLimits())
    make = () -> discover_closure(space, generator, seeds; limits)
    result = @inferred make()
    measurement = warm_measure(make)
    push!(rows, (
        name = name,
        status = result.status,
        basis = length(result.basis.terms),
        commutators = result.commutators,
        max_degree = result.basis.max_degree,
        offending_degree = result.offending_degree,
        seconds = measurement.seconds,
        bytes = measurement.bytes,
        concrete = isconcretetype(typeof(result)) && isconcretetype(eltype(result.action)),
    ))
    return result
end

rows = NamedTuple[]

hf = FockSpace(:f)
a = Destroy(hf, :a)
number = record_fixture!(rows, "number rotation", hf, a' * a, [a])
@test number.status == :closed
@test length(number.basis.terms) == 2
@test number.basis.max_degree == 1

single_squeeze = record_fixture!(rows, "one-mode squeeze", hf, (a'^2 - a^2) / 2, [a])
@test single_squeeze.status == :closed
@test length(single_squeeze.basis.terms) == 2

h2 = FockSpace(:left) ⊗ FockSpace(:right)
left = Destroy(h2, :a, 1)
right = Destroy(h2, :b, 2)
two_rotation = record_fixture!(
    rows, "two-mode rotation", h2, left' * right + right' * left, [left, right],
)
@test two_rotation.status == :closed
@test length(two_rotation.basis.terms) == 4

two_squeeze = record_fixture!(
    rows, "two-mode squeeze", h2, left' * right' - left * right, [left, right],
)
@test two_squeeze.status == :closed
@test length(two_squeeze.basis.terms) == 4

hp = PauliSpace(:p)
σ = [Pauli(hp, :σ, axis) for axis in 1:3]
pauli = record_fixture!(rows, "Pauli rotation", hp, σ[3], σ)
@test pauli.status == :closed
@test length(pauli.basis.terms) == 3

hs = SpinSpace(:s)
S = [Spin(hs, :S, axis) for axis in 1:3]
spin = record_fixture!(rows, "Spin rotation", hs, S[3], S)
@test spin.status == :closed
@test length(spin.basis.terms) == 3

hn = NLevelSpace(:atom, 3, 2)
σ12 = Transition(hn, :σ, 1, 2)
nlevel_seeds = fundamental_operators(hn)
nlevel = record_fixture!(rows, "N-level mixing", hn, σ12 + σ12', nlevel_seeds)
@test nlevel.status == :closed
@test length(nlevel.basis.terms) == 8
@test nlevel.basis.max_degree == 1

idx = Index(hf, :i, 5, hf)
indexed_a = IndexedOperator(a, idx(2))
indexed = record_fixture!(
    rows, "concrete indexed number rotation", hf,
    indexed_a' * indexed_a, [indexed_a],
)
@test indexed.status == :closed
@test length(indexed.basis.terms) == 2

kerr = record_fixture!(
    rows, "Kerr refusal", hf, a'^2 * a^2, [a];
    limits = ClosureLimits(max_basis = 64, max_degree = 7, max_commutators = 64),
)
@test kerr.status == :max_degree
@test kerr.offending_degree == 9

hfb = FockSpace(:field) ⊗ SpinSpace(:spin)
af = Destroy(hfb, :a, 1)
SX = Spin(hfb, :S, 1, 2)
SY = Spin(hfb, :S, 2, 2)
SZ = Spin(hfb, :S, 3, 2)
spin_boson = record_fixture!(
    rows, "spin-boson refusal", hfb, (af + af') * SZ, [af, SX, SY, SZ];
    limits = ClosureLimits(max_basis = 48, max_degree = 6, max_commutators = 48),
)
@test spin_boson.status in (:max_basis, :max_degree, :max_commutators)

# Independent algebraic oracles, not a second closure implementation.
@test isequal(commutator(a' * a, a), 1 * (-a))
@test isequal(commutator(a' * a, a'), 1 * a')
@test isequal(commutator(σ[3], σ[1]), 1 * (2im * σ[2]))
@test isequal(commutator(S[3], S[1]), 1 * (1im * S[2]))
@test iszero(expand_completeness(
    commutator(σ12 + σ12', σ12) -
        (Transition(hn, :σ, 2, 2) - Transition(hn, :σ, 1, 1)),
))

# Metadata/inference pressure tests.
renamed = SecondQuantizedAlgebra.rename(a, :renamed)
renamed_context = infer_operator_context([renamed, renamed'])
@test string(renamed_context.slots[Int32(1)].family) == "FAMILY_FOCK"

duplicate_context = infer_operator_context([left, right])
@test duplicate_context.contiguous
@test length(duplicate_context.slots) == 2

indexed_context = infer_operator_context([indexed_a])
@test indexed_context.slots[Int32(1)].indexed

nlevel_context = infer_operator_context([Transition(hn, :T, 1, 3)])
@test nlevel_context.slots[Int32(1)].n == 3
@test nlevel_context.slots[Int32(1)].ground == 2

hc = CollectiveNLevelSpace(:collective, 4)
collective_context = infer_operator_context([CollectiveTransition(hc, :S, 1, 2)])
@test !collective_context.slots[Int32(1)].complete
@test collective_context.slots[Int32(1)].n == 0
@test isconcretetype(typeof(collective_context))

storage = storage_comparison(nlevel.basis.terms; repetitions = 200)

println("CLOSURE_ROWS_BEGIN")
for row in rows
    println(row)
end
println("CLOSURE_ROWS_END")
println("STORAGE_COMPARISON=", storage)
println("RESULT_TYPE=", typeof(nlevel))
println("ACTION_TYPE=", typeof(nlevel.action))
println("NLEVEL_CONTEXT=", nlevel_context)
println("COLLECTIVE_CONTEXT=", collective_context)
println("PROTOTYPE_TESTS_PASSED")

