using SecondQuantizedAlgebra
using QuantumOpticsBase
using Test

dat(x) = dense(x).data

struct MockNumericBackend <: SecondQuantizedAlgebra.NumericBackend end
struct MockEagerOperator
    data::Matrix{ComplexF64}
end
mock_identity(n) = ComplexF64[i == j for i in 1:n, j in 1:n]
SecondQuantizedAlgebra.numeric_basis(::MockNumericBackend, ::FockSpace, cutoff) = Int(cutoff) + 1
SecondQuantizedAlgebra.numeric_num_subsystems(::MockNumericBackend, ::Int) = 1
SecondQuantizedAlgebra.numeric_subbasis(::MockNumericBackend, n::Int, slot::Int) =
    slot == 1 ? n : throw(ArgumentError("mock basis has no subsystem slot $slot"))
SecondQuantizedAlgebra.numeric_operator(
    ::MockNumericBackend, ::SecondQuantizedAlgebra.Op, n::Int,
) = mock_identity(n)
SecondQuantizedAlgebra.numeric_embed(::MockNumericBackend, ::Int, ::Int, leaf) = leaf
SecondQuantizedAlgebra.numeric_identity(::MockNumericBackend, n::Int) = mock_identity(n)
function SecondQuantizedAlgebra.numeric_assemble(::MockNumericBackend, n::Int, terms)
    result = zeros(ComplexF64, n, n)
    for (coefficient, factors) in terms
        term = isempty(factors) ? mock_identity(n) : foldl(*, factors)
        result .+= coefficient .* term
    end
    return result
end
SecondQuantizedAlgebra.numeric_materialize(::MockNumericBackend, assembled, ::Nothing) =
    MockEagerOperator(assembled)
SecondQuantizedAlgebra.numeric_materialize(
    ::MockNumericBackend, assembled, ::typeof(identity),
) = assembled

@testset "Numeric backend interface" begin
    @testset "third-party static backend contract" begin
        h = FockSpace(:mock)
        @qnumbers a::Destroy(h)

        eager = to_numeric(2 * a, h, 3; backend = MockNumericBackend())
        lazy = to_numeric(2 * a, h, 3; backend = MockNumericBackend(), op_type = identity)
        @test eager isa MockEagerOperator
        @test eager.data == 2 .* mock_identity(4)
        @test lazy == eager.data
    end

    @testset "public backend hooks and product dims validation" begin
        b = FockBasis(3)
        ψ = fockstate(b, 0)
        @test numeric_backend(ψ) isa QuantumOpticsBackend
        @test numeric_basis(ψ) == b
        @test SecondQuantizedAlgebra.numeric_num_subsystems(QuantumOpticsBackend(), b) == 1

        h = SecondQuantizedAlgebra.var"⊗"(FockSpace(:a), FockSpace(:b))
        a = Destroy(h, :a, 1)
        @test_throws ArgumentError to_numeric(a, h, [2]; backend = QuantumOpticsBackend())
        @test_throws ArgumentError to_numeric(a, h, [2, 3, 99]; backend = QuantumOpticsBackend())
    end

    @testset "HilbertSpace form (uniform entry)" begin
        h = FockSpace(:fock)
        @qnumbers a::Destroy(h)
        b = FockBasis(7)
        ψ = coherentstate(b, 0.2 + 0.1im)

        @test dat(to_numeric(a' * a, h, 7; backend = QuantumOpticsBackend())) ==
            dat(to_numeric(a' * a, b))
        @test numeric_average(a' * a, ψ) ≈
            expect(to_numeric(a' * a, h, 7; backend = QuantumOpticsBackend()), ψ)

        hp = SecondQuantizedAlgebra.var"⊗"(FockSpace(:c), NLevelSpace(:atom, 3, 1))
        ap = Destroy(hp, :a, 1)
        σp = Transition(hp, :σ, 1, 2, 2)
        bp = QuantumOpticsBase.var"⊗"(FockBasis(4), NLevelBasis(3))
        @test dat(to_numeric(ap' * σp, hp, (4, 3); backend = QuantumOpticsBackend())) ==
            dat(to_numeric(ap' * σp, bp))
    end
end
