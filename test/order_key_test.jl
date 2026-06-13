using SecondQuantizedAlgebra
using Test
import SecondQuantizedAlgebra: QSym, Index, IndexedOperator

# `order_key(op)` is the total, identity-faithful ordering key: a comparable tuple that
# orders operators and ties two operators exactly when they are `isequal`.
const order_key = SecondQuantizedAlgebra.order_key

@testset "order_key" begin
    # One product space exercising every concrete QSym subtype. Space indices:
    # c=1 (Fock), atom=2 (NLevel), p=3 (Pauli), s=4 (Spin), q=5 (Phase).
    h = FockSpace(:c) ⊗ NLevelSpace(:atom, 3, 1) ⊗ PauliSpace(:p) ⊗ SpinSpace(:s) ⊗ PhaseSpace(:q)
    N = 4
    ic = Index(h, :i, N, 1)
    jc = Index(h, :j, N, 1)

    a = Destroy(h, :a, 1)
    ad = a'                      # Create, same site
    b = Destroy(h, :b, 1)        # different name
    ai = IndexedOperator(a, ic)  # indexed
    aj = IndexedOperator(a, jc)  # different index name
    σ12 = Transition(h, :σ, 1, 2, 2)
    σ13 = Transition(h, :σ, 1, 3, 2)  # differs only in level j
    px = Pauli(h, :σ, 1, 3)
    py = Pauli(h, :σ, 2, 3)           # differs only in axis
    Sx = Spin(h, :S, 1, 4)
    Sy = Spin(h, :S, 2, 4)            # differs only in axis
    x = Position(h, :x, 5)
    p = Momentum(h, :p, 5)            # differs only in type from a Position-shaped key

    ops = QSym[a, ad, b, ai, aj, σ12, σ13, px, py, Sx, Sy, x, p]

    @testset "ties exactly when isequal" begin
        for o1 in ops, o2 in ops
            @test (order_key(o1) == order_key(o2)) == isequal(o1, o2)
        end
    end

    @testset "distinguishes type-specific fields" begin
        @test order_key(σ12) != order_key(σ13)  # levels
        @test order_key(px) != order_key(py)    # axis
        @test order_key(Sx) != order_key(Sy)    # axis
        @test order_key(x) != order_key(p)      # Position vs Momentum
        @test order_key(ai) != order_key(aj)    # index name
    end

    @testset "every concrete QSym subtype has a method" begin
        for T in (Destroy, Create, Transition, Pauli, Spin, Position, Momentum)
            @test hasmethod(order_key, Tuple{T})
            # a dedicated method, not the generic QSym fallback (which would drop fields)
            @test which(order_key, Tuple{T}).sig !== Tuple{typeof(order_key), QSym}
        end
    end

    @testset "is a strict total order (trichotomy)" begin
        for o1 in ops, o2 in ops
            k1, k2 = order_key(o1), order_key(o2)
            @test (k1 < k2) + (k2 < k1) + (k1 == k2) == 1
        end
    end
end
