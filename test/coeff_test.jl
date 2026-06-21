using Test
using SecondQuantizedAlgebra
using Symbolics: @variables, Num
import SecondQuantizedAlgebra: Coeff, CNum, _to_cnum, _to_complex, to_num, _is_native,
    _conj_cnum, _mul_cnum, _add_cnum, _neg_cnum, _iszero_cnum,
    _CNUM_ONE, _CNUM_ZERO, _CNUM_NEG1, _CNUM_IM, _NUM_ZERO

# Coefficients carry a native `ComplexF64` fast path and a `Complex{Num}` symbolic
# fallback. These tests pin the invariants the rest of the package relies on:
# concrete numbers stay native, genuine symbols stay symbolic, materialization is
# faithful, and equality / hashing are canonical.

@testset "Coeff" begin
    @testset "native classification" begin
        for x in (0, 1, -1, 2, 1.5, 0.7, im, -im, 2 + 3im, 1 // 2, -3 // 4)
            @test _is_native(_to_cnum(x))
        end
        @test CNum === Coeff
    end

    @testset "exactness gate keeps inexact rationals / bignums symbolic" begin
        # 1//3 has no exact ComplexF64, so it must stay on the symbolic path.
        @test !_is_native(_to_cnum(1 // 3))
        @test !_is_native(_to_cnum(big(2)^70))
        # round-trip is still numerically faithful
        @test _to_complex(_to_cnum(1 // 3)) ≈ 1 / 3
    end

    @testset "symbolic fallback" begin
        @variables g
        c = _to_cnum(Complex(Num(g), _NUM_ZERO))
        @test !_is_native(c)
        @test isequal(to_num(c), Complex(Num(g), _NUM_ZERO))
        # a symbolic expression that folds to a number lands back on the native path
        @test _is_native(_to_cnum(Num(g) - Num(g) + 4))
    end

    @testset "to_num round-trip and integer display" begin
        @test isequal(to_num(_to_cnum(2)), Complex(Num(2), Num(0)))
        @test isequal(to_num(_to_cnum(2)), Complex(Num(2), Num(0)))   # integer, not 2.0
        @test string(real(to_num(_to_cnum(2)))) == "2"
        @test string(real(to_num(_to_cnum(0.7)))) == "0.7"
        @test isequal(to_num(_to_cnum(2 + 3im)), Complex(Num(2), Num(3)))
    end

    @testset "_to_complex matches the value" begin
        @test _to_complex(_to_cnum(2)) === ComplexF64(2)
        @test _to_complex(_to_cnum(2 + 3im)) === ComplexF64(2, 3)
        @test _to_complex(_CNUM_ZERO) === ComplexF64(0)
    end

    @testset "conjugation and signed-zero normalization" begin
        # conj(2) produces a -0.0 imaginary part; it must normalize so structurally
        # equal coefficients stay isequal AND hash identically (dict correctness).
        c = _to_cnum(2)
        @test isequal(_conj_cnum(c), c)
        @test hash(_conj_cnum(c)) == hash(c)
        @test isequal(_neg_cnum(_neg_cnum(c)), c)
        @test hash(_neg_cnum(_neg_cnum(c))) == hash(c)
        # complex conjugation flips the sign of the imaginary part
        @test isequal(_conj_cnum(_to_cnum(2 + 3im)), _to_cnum(2 - 3im))
    end

    @testset "native arithmetic stays native and exact" begin
        @test _is_native(_mul_cnum(_to_cnum(2), _to_cnum(3)))
        @test _mul_cnum(_to_cnum(2), _to_cnum(3)) == 6
        @test _add_cnum(_to_cnum(2), _to_cnum(3)) == 5
        @test _mul_cnum(_CNUM_IM, _CNUM_IM) == -1
        @test _iszero_cnum(_add_cnum(_to_cnum(2), _neg_cnum(_to_cnum(2))))
    end

    @testset "mixed comparison with Number / Complex{Num}" begin
        @test _to_cnum(2) == 2
        @test 2 == _to_cnum(2)
        @test _to_cnum(2) == 2 + 0im
        @test _to_cnum(2) == Complex(Num(2), Num(0))
        @test isequal(_to_cnum(2), 2)
        @test !(_to_cnum(2) == 3)
    end

    @testset "scalar arithmetic operators on Coeff" begin
        @test _to_cnum(2) + _to_cnum(3) == 5
        @test _to_cnum(2) + 3 == 5
        @test 3 + _to_cnum(2) == 5
        @test _to_cnum(5) - 2 == 3
        @test -_to_cnum(2) == -2
        @test _to_cnum(2) * 3 == 6
        @test _to_cnum(6) / 2 == 3
    end

    @testset "QAdd equality is canonical across construction paths" begin
        h = FockSpace(:f)
        a = Destroy(h, :a)
        m = 2 * a' * a
        @test isequal(m', m)               # adjoint reproduces the same coefficient
        @test hash(m') == hash(m)
        @test isequal(simplify(m), m)
    end
end
